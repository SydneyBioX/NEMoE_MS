library(glmnet)



glmnet_cv <- function(X, y){
  X <- as.matrix(X)

  temp <- cv.glmnet(x = X, y = y, family = "binomial")

  res <- glmnet(x = X, y = y, family = "binomial", lambda = temp$lambda.min)

  return(res)
}

predict_two_stage <- function(X_new, Z_new, y_new = NULL, twostage, method = "svm", output = "acc"){

  K <- length(twostage$model)
  n <- nrow(X_new)

  center <- twostage$center

  dist_all <- as.matrix(dist(rbind(center, Z_new)))

  dist_sel <- as.matrix(dist_all[(K+1):nrow(dist_all), 1:K])

  r_i <- apply(dist_sel, 1, which.min)

  y_pred <- rep(0, nrow(X_new))

  m <- twostage$model


  if(n == 1){

    if(method == "svm"){

      temp_svm <- attr(predict(m[[r_i]], X_new, type = "response", probability = T), "probabilities")
      y_pred = temp_svm[,"1"]
    }else if(method == "random forest"){
      y_pred = predict(m[[r_i]], X_new, type = "prob")[,"1"]
    }else{
      if("glm" %in% class(m[[r_i]])){
        X_new <- as.data.frame(X_new)
        y_pred = predict(m[[r_i]], X_new, type = "response")
      }else{
        y_pred = predict(m[[r_i]], X_new, type = "response")
      }
    }

  }else{
    for(i in 1:K){

      X_k <- X_new[r_i == i,]
      if(is.vector(X_k)){
        X_k <- matrix(X_k, nrow = 1)
      }

      if(method == "svm"){

        temp_svm <- attr(predict(m[[i]], X_k, type = "response", probability = T), "probabilities")
        y_pred[r_i == i] = temp_svm[,"1"]

      }else if(method == "random forest"){

        y_pred[r_i == i] = predict(m[[i]], X_k, type = "prob")[,2]

      }else if(method == "glmnet"){

        if("glm" %in% class(m[[i]])){
          X_k <- as.data.frame(X_k)
          y_pred[r_i == i] = predict(m[[i]], X_k, type = "response")
        }else{
          y_pred[r_i == i] = predict(m[[i]], X_k, type = "response")
        }
      }
    }
  }

  if(output == "acc"){
    y_hard <- round(y_pred)
    return(sum(y_hard == y_new)/length(y_new))
  }else{
    return(y_pred)
  }
}

two_stage <- function(X, Z, y, K = 2, method = "svm"){

  km_model <- kmeans(Z, K)

  latent <- km_model$cluster

  center <- km_model$centers

  m <- list()

  for(i in 1:K){

    X_k <- X[latent == i,]
    y_k <- y[latent == i]

    if(method == "svm"){
      m[[i]] <- svm(x = X_k, y = as.factor(y_k), kernel = "radial", probability = TRUE)
    }else if(method == "random forest"){
      m[[i]] <- randomForest(x = X_k, y = as.factor(y_k), ntree = 500)
    }else{

      if(class_check(y_k)){
        temp <- try(cv.glmnet(x = X_k, y = as.factor(y_k), family = "binomial"))
        m[[i]] <- glmnet(x = X_k, y = as.factor(y_k), family = "binomial", lambda = temp$lambda.min)
      }else{
        m[[i]] <- glm(as.factor(y_k) ~ X_k, family = "binomial")
      }

    }
  }
  return(list(center = center, model = m))
}


rocpdf<- function(df, method = "RMoE", gd = "actual", model_name = method){

  pred <- df[,method]
  truth <- df[,gd]
  predob<- prediction(pred, truth)
  perf.auc<- performance(predob, measure = 'auc', x.measure = 'cutoff')
  perf<- performance(predob, 'tpr','fpr')
  df<- data.frame(x = attributes(perf)$x.values[[1]],y = attributes(perf)$y.values[[1]], model = model_name)
  return(list(df = df, auc = perf.auc@y.values[[1]]))
}

evaluation <- function(n = 200, p = 50, q = 30, c_e = 2, c_g = 2,
                       eta = 0, rho = 0, K = 2, K0 = 2,  version = "probit",
                       B= 100, num_restart = 0, lambda1 = 0.01,
                       lambda2 = 0.015, verbose = F, init = "kmeans"){

  acc <- matrix(0, nrow = B, ncol = 7)
  Sigma <- (1 - rho)*diag(rep(1,q)) + rho*matrix(1,nrow = q, ncol = q)

  for(b in 1:B){

    data_2 <- genNEMoE(n = 2*n, p = p, q = q, c_e = c_e, c_g = c_g,
                       eta = eta, Sigma = Sigma, K = K, gen_Micro = "dm",
                       p_L = c(), link = version, prev_filt = 0.2)

    mu_X_train <- colMeans(data_2$X_list[[1]][1:n,])
    mu_Z_train <- colMeans(data_2$W[1:n,])
    sd_X_train <- apply(data_2$X_list[[1]], 2, sd)
    sd_Z_train <- apply(data_2$W, 2, sd)

    X_train <- data_2$X_list[[1]][1:n,]
    X_train_in <- scale(X_train)
    Z_train <- data_2$W[1:n,]
    Z_train_in <- scale(Z_train)
    X1_train_in <- cbind(X_train_in, Z_train_in)
    y_train <- data_2$y[1:n]

    X_test <- data_2$X_list[[1]][(n+1):(2*n),]
    X_test_in <- sapply(1:ncol(X_test), function(col){(X_test[,col] - mu_X_train[col])/sd_X_train[col]})
    Z_test <- data_2$W[(n+1):(2*n),]
    Z_test_in <- sapply(1:ncol(Z_test), function(col){(Z_test[,col] - mu_Z_train[col])/sd_Z_train[col]})
    X1_test_in <- cbind(X_test_in, Z_test_in)
    y_test <- data_2$y[(n+1):(2*n)]

    NEMoE_simu <- NEMoE_buildFromList(X_train_in, Z_train_in, y_train, K = K0, lambda1 = lambda1,
                                      lambda2 = lambda2, adapt =T, btr = T,
                                      alpha1 = 0.5, alpha2 = 0.5 , itmax = 1e2,
                                      verbose= verbose, init = init)

    NEMoE_simu <- fitNEMoE(NEMoE_simu, num_restart = 0)

    #        refit_MoE <- refit_cv(X_train_in, Z_train_in, y = y_train, EM_res = res_MoE)

    res_glmnet1 <- try(glmnet_cv(X = X_train_in, y = y_train))

    res_glmnet2 <- try(two_stage(X = X_train_in, Z = Z_train_in, y = y_train, K = 2, method = "glmnet"))

    res_svm1 <- svm(x = X_train_in, y = as.factor(y_train), kernel = "radial", probability = T )

    res_svm2 <- try(two_stage(X = X_train_in, Z = Z_train_in, y = y_train, K = 2, method = "svm"))

    res_rf1 <- randomForest(x = X_train_in, y = as.factor(y_train), ntree = 500)

    res_rf2 <- try(two_stage(X = X_train_in, Z = Z_train_in, y = y_train, K = 2, method = "random forest"))
    pred_MoE <- NEMoE_predict(NEMoE = NEMoE_simu, X_new = list(X_test_in),
                             Z_new = Z_test_in, full = F, transform = F)
#    latent_pred <- cvtLabel(pred_MoE$gating_prob, idx = F)
#    label_pred <- rowSums(latent_pred * pred_MoE$experts_prob[[1]])
    #        acc_refit_MoE <- sum(predict_EM(X_new = X_test_in, Z_new = Z_test_in, EM_res = MoE_refit, output = "hard") == y_test)/n
    acc_MoE <- sum(round(pred_MoE) == y_test)/n
    if( "try-error" %in% class(res_glmnet1)){
      pred_glmnet1 <- rep(which.max(table(y_train)) - 1,n)
    }else{
      pred_glmnet1 <- predict(res_glmnet1, X_test_in, type = "response")
    }
    acc_glmnet1 <- sum(round(pred_glmnet1) == y_test)/ length(y_test)
    pred_svm1 <- as.numeric(as.character(predict(res_svm1, X_test_in, probability = T)))
    acc_svm1 <- sum(round(pred_svm1) == y_test)/length(y_test)
    pred_rf1 <- predict(res_rf1, X_test_in, type = "prob")[,2]
    acc_rf1 <- sum(round(pred_rf1) == y_test)/length(y_test)

    if( "try-error" %in% class(res_glmnet2)){
      pred_glmnet2 <- pred_glmnet1
      acc_glmnet2 <- acc_glmnet1
    }else{
      acc_glmnet2 <- predict_two_stage(X_new = X_test_in, Z_new = Z_test_in, y_new = y_test, twostage = res_glmnet2, method = "glmnet", output = "acc")
    }

    if( "try-error" %in% class(res_svm2)){
      pred_svm2 <- pred_svm1
      acc_svm2 <- acc_svm1
    }else{
      acc_svm2 <- predict_two_stage(X_new = X_test_in, Z_new = Z_test_in, y_new = y_test, twostage = res_svm2, method = "svm", output = "acc")

    }

    if( "try-error" %in% class(res_rf2)){
      pred_rf2 <- pred_rf1
      acc_rf2 <- acc_rf1
    }else{
      acc_rf2 <- predict_two_stage(X_new = X_test_in, Z_new = Z_test_in, y_new = y_test, twostage = res_rf2, method = "random forest", output = "acc")
    }

    acc[b,] <- c(acc_MoE, acc_glmnet1, acc_glmnet2, acc_svm1, acc_svm2, acc_rf1, acc_rf2)

    print(paste0("b=", b,":", paste(acc[b,], collapse = ",")))
  }
  colnames(acc) <- c("NEMoE", "sLR", "sLR II", "svm", "svm II", "rf", "rf II")
  return(acc)
}

evaluation_K <- function(n = 400, p = 50, q = 30, c_e = 1, c_g = 2,
                         eta = 0, rho = 0, K = 3,  version = "probit",
                         B= 100, num_restart = 0, verbose = F, init = "kmeans"){

  acc <- matrix(0, nrow = B, ncol = 7)
  Sigma <- (1 - rho)*diag(rep(1,q)) + rho*matrix(1,nrow = q, ncol = q)

  for(b in 1:B){

    data_2 <- genNEMoE(n = 2*n, p = p, q = q, c_e = c_e, c_g = c_g,
                       eta = eta, Sigma = Sigma, K = K, gen_Micro = "dm",
                       p_L = c(), link = version, prev_filt = 0.2)

    mu_X_train <- colMeans(data_2$X_list[[1]][1:n,])
    mu_Z_train <- colMeans(data_2$W[1:n,])
    sd_X_train <- apply(data_2$X_list[[1]], 2, sd)
    sd_Z_train <- apply(data_2$W, 2, sd)

    X_train <- data_2$X_list[[1]][1:n,]
    X_train_in <- scale(X_train)
    Z_train <- data_2$W[1:n,]
    Z_train_in <- scale(Z_train)
    X1_train_in <- cbind(X_train_in, Z_train_in)
    y_train <- data_2$y[1:n]

    X_test <- data_2$X_list[[1]][(n+1):(2*n),]
    X_test_in <- sapply(1:ncol(X_test), function(col){(X_test[,col] - mu_X_train[col])/sd_X_train[col]})
    Z_test <- data_2$W[(n+1):(2*n),]
    Z_test_in <- sapply(1:ncol(Z_test), function(col){(Z_test[,col] - mu_Z_train[col])/sd_Z_train[col]})
    X1_test_in <- cbind(X_test_in, Z_test_in)
    y_test <- data_2$y[(n+1):(2*n)]

    print("Start glmnet")
    res_glmnet1 <- try(glmnet_cv(X = X_train_in, y = y_train))
    print("Start glmnet2")
    res_glmnet2 <- try(two_stage(X = X_train_in, Z = Z_train_in, y = y_train, K = 2, method = "glmnet"))
    print("Start glmnet3")
    res_glmnet3 <- try(two_stage(X = X_train_in, Z = Z_train_in, y = y_train, K = 3, method = "glmnet"))
    print("Start glmnet4")
    res_glmnet4 <- try(two_stage(X = X_train_in, Z = Z_train_in, y = y_train, K = 4, method = "glmnet"))
    print("Start NEMoE2")
    NEMoE_simu_2 <- NEMoE_buildFromList(X_train_in, Z_train_in, y_train, K = 2, lambda1 = 0.01,
                                      lambda2 = 0.01, adapt =T, btr = T,
                                      alpha1 = 1, alpha2 = 1 ,
                                      verbose= verbose, init = init)

    NEMoE_simu_2 <- fitNEMoE(NEMoE_simu_2, num_restart = 0)
    print("Start NEMoE3")
    NEMoE_simu_3 <- NEMoE_buildFromList(X_train_in, Z_train_in, y_train, K = 3, lambda1 = 0.008,
                                        lambda2 = 0.008, adapt =T, btr = T,
                                        alpha1 = 1, alpha2 = 1 , itmax = 50,
                                        verbose= verbose, init = init)

    NEMoE_simu_3 <- fitNEMoE(NEMoE_simu_3, num_restart = 0)
    print("Start NEMoE4")
    NEMoE_simu_4 <- NEMoE_buildFromList(X_train_in, Z_train_in, y_train, K = 4, lambda1 = 0.006,
                                        lambda2 = 0.006, adapt =T, btr = T,
                                        alpha1 = 1, alpha2 = 1 , itmax = 50,
                                        verbose= verbose, init = init)

    NEMoE_simu_4 <- fitNEMoE(NEMoE_simu_4, num_restart = 0)

#    NEMoE_simu_5 <- NEMoE_buildFromList(X_train_in, Z_train_in, y_train, K = 5, lambda1 = 0.005,
#                                        lambda2 = 0.005, adapt =T, btr = T,
#                                        alpha1 = 1, alpha2 = 1 ,
#                                        verbose= verbose, init = init)

#    NEMoE_simu_5 <- fitNEMoE(NEMoE_simu_5, num_restart = 0)

    pred_MoE_2 <- NEMoE_predict(NEMoE = NEMoE_simu_2, X_new = list(X_test_in),
                              Z_new = Z_test_in, full = F, transform = F)

    acc_MoE_2 <- sum(round(pred_MoE_2) == y_test)/n

    pred_MoE_3 <- NEMoE_predict(NEMoE = NEMoE_simu_3, X_new = list(X_test_in),
                              Z_new = Z_test_in, full = F, transform = F)
    acc_MoE_3 <- sum(round(pred_MoE_3) == y_test)/n

    pred_MoE_4 <- NEMoE_predict(NEMoE = NEMoE_simu_4, X_new = list(X_test_in),
                              Z_new = Z_test_in, full = F, transform = F)
    acc_MoE_4 <- sum(round(pred_MoE_4) == y_test)/n

#    pred_MoE_5 <- NEMoE_predict(NEMoE = NEMoE_simu_5, X_new = list(X_test_in),
#                              Z_new = Z_test_in, full = F, transform = F)
#    acc_MoE_5 <- sum(round(pred_MoE_5) == y_test)/n

    if( "try-error" %in% class(res_glmnet1)){
      pred_glmnet1 <- rep(which.max(table(y_train)) - 1,n)
    }else{
      pred_glmnet1 <- predict(res_glmnet1, X_test_in, type = "response")
    }
    acc_glmnet1 <- sum(round(pred_glmnet1) == y_test)/ length(y_test)

    if( "try-error" %in% class(res_glmnet2)){
      pred_glmnet2 <- pred_glmnet1
      acc_glmnet2 <- acc_glmnet1
    }else{
      acc_glmnet2 <- predict_two_stage(X_new = X_test_in, Z_new = Z_test_in, y_new = y_test, twostage = res_glmnet2, method = "glmnet", output = "acc")
    }

    if( "try-error" %in% class(res_glmnet3)){
      pred_glmnet3 <- pred_glmnet1
      acc_glmnet3 <- acc_glmnet1
    }else{
      acc_glmnet3 <- predict_two_stage(X_new = X_test_in, Z_new = Z_test_in, y_new = y_test, twostage = res_glmnet3, method = "glmnet", output = "acc")
    }

    if( "try-error" %in% class(res_glmnet4)){
      pred_glmnet4 <- pred_glmnet1
      acc_glmnet4 <- acc_glmnet1
    }else{
      acc_glmnet4 <- predict_two_stage(X_new = X_test_in, Z_new = Z_test_in, y_new = y_test, twostage = res_glmnet4, method = "glmnet", output = "acc")
    }

    acc[b,] <- c(acc_glmnet1, acc_MoE_2, acc_glmnet2, acc_MoE_3, acc_glmnet3, acc_MoE_4, acc_glmnet4)

    print(paste0("b=", b,":", paste(acc[b,], collapse = ",")))
  }
  colnames(acc) <- c("sLR", "NEMoE II", "sLR II", "NEMoE III", "sLR III ", "NEMoE IV", "sLR IV")
  return(acc)
}


plt_function <- function(tab_list, label_list, method_list = c("NEMoE", "sLR", "sLR II", "SVM", "SVM II", "RF", "RF II")){

  K <- length(tab_list)

  plt_tab <- data.frame(B = 0, method = "", acc = 0, n = 0, p = 0, q = 0, c_e = 2, c_g = 2, rho = 0, model = "weakly")

  for(i in 1:K){
    colnames(tab_list[[i]]) <- method_list

    tab_temp <- reshape2::melt(tab_list[[i]])

    colnames(tab_temp) <- c("B", "method", "acc")

    tab_temp$n <- label_list[i, 1]
    tab_temp$p <- label_list[i, 2]
    tab_temp$q <- label_list[i, 3]
    tab_temp$c_e <- label_list[i, 4]
    tab_temp$c_g <- label_list[i, 5]
    tab_temp$rho <- label_list[i, 6]
    tab_temp$model <- label_list[i, 7]

    plt_tab <- rbind(plt_tab, tab_temp)
  }
  plt_tab <- plt_tab[-1,]
  return(plt_tab)
}

plt_function_K <- function(tab_list, label_list, method_list = c("sLR", "NEMoE", "NEMoE II", "sLR II", "NEMoE III", "sLR III", "NEMoE IV", "sLR IV")){

  K <- length(tab_list)

  plt_tab <- data.frame(B = 0, method = "", acc = 0, n = 0, p = 0, q = 0, c_e = 2, c_g = 2, rho = 0, model = "weakly")

  for(i in 1:K){
    colnames(tab_list[[i]]) <- method_list

    tab_temp <- reshape2::melt(tab_list[[i]])

    colnames(tab_temp) <- c("B", "method", "acc")

    tab_temp$n <- label_list[i, 1]
    tab_temp$p <- label_list[i, 2]
    tab_temp$q <- label_list[i, 3]
    tab_temp$c_e <- label_list[i, 4]
    tab_temp$c_g <- label_list[i, 5]
    tab_temp$rho <- label_list[i, 6]
    tab_temp$model <- label_list[i, 7]

    plt_tab <- rbind(plt_tab, tab_temp)
  }
  plt_tab <- plt_tab[-1,]
  return(plt_tab)
}

class_check <- function(y, min_obs = 7){

  y <- as.matrix(y)

  n <- nrow(y)
  K <- ncol(y)

  y <- as.numeric(as.character(y))

  y <- matrix(y, nrow = n, ncol = K)

  if(!length(y)){
    return(FALSE)
  }

  if(K == 1){
    y <- cbind(y, 1 - y)
  }

  y_n <- colSums(y)

  if(any(y_n < min_obs)){
    return(FALSE)
  }else{
    return(TRUE)
  }

}



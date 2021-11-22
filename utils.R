library(glmnet)

glmnet_cv <- function(X, y){
  X <- as.matrix(X)

  temp <- suppressWarnings(cv.glmnet(x = X, y = y, family = "binomial"))

  res <- suppressWarnings(glmnet(x = X, y = y, family = "binomial", lambda = temp$lambda.min))

  return(res)
}

predict_two_stage <- function(X_new, Z_new, y_new = NULL, twostage, method = "svm", output = "acc"){

  K <- length(twostage$model)
  n <- nrow(X_new)

  center <- twostage$center

  dist_all <- as.matrix(dist(rbind(center, Z_new)))

  dist_sel <- as.matrix(dist_all[(K+1):nrow(dist_all), 1:K, drop = F])

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

      X_k <- X_new[r_i == i,,drop = F]
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

evaluation_K <- function(n = 400, p = 50, q = 30, c_e = 2, c_g = 2,
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

    res_glmnet1 <- try(glmnet_cv(X = X_train_in, y = y_train))
    res_glmnet2 <- try(two_stage(X = X_train_in, Z = Z_train_in, y = y_train, K = 2, method = "glmnet"))
    res_glmnet3 <- try(two_stage(X = X_train_in, Z = Z_train_in, y = y_train, K = 3, method = "glmnet"))
    res_glmnet4 <- try(two_stage(X = X_train_in, Z = Z_train_in, y = y_train, K = 4, method = "glmnet"))

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

evaluate_level <- function(n = 500, p_L = c(20,50,80,100), eta = 0,
                           rho = 0, K = 2, K0 = 2,  version = "probit",
                           B= 100, num_restart = 0, lambda1 = 0.01,
                           lambda2 = 0.015, verbose = F, init = "glmnet"){

  ari <- matrix(0, nrow = B, ncol = 5)
  Sigma <- (1 - rho)*diag(rep(1,q)) + rho*matrix(1,nrow = q, ncol = q)

  for(b in 1:B){

    data_2 <- genNEMoE(n = 2*n, p = p, q = q, c_e = c_e, c_g = c_g,
                       eta = eta, Sigma = Sigma, K = K, gen_Micro = "mgauss",
                       p_L = c(20,40,80), link = version, prev_filt = 0)
    data_2$X_list <- lapply(data_2$X_list, scale)
    data_2$W <- scale(data_2$W)

    X1_train <- data_2$X_list[[1]][1:n,]
    X2_train <- data_2$X_list[[2]][1:n,]
    X3_train <- data_2$X_list[[3]][1:n,]
    X4_train <- data_2$X_list[[4]][1:n,]
    X_list_train <- list(X1_train, X2_train, X3_train, X4_train)
    Z_train <- data_2$W[1:n,]
    y_train <- data_2$y[1:n]

    X1_test <- data_2$X_list[[1]][(n+1):(2*n),]
    X2_test <- data_2$X_list[[2]][(n+1):(2*n),]
    X3_test <- data_2$X_list[[3]][(n+1):(2*n),]
    X4_test <- data_2$X_list[[4]][(n+1):(2*n),]
    X_list_test <- list(X1_test, X2_test, X3_test, X4_test)
    Z_test <- data_2$W[(n+1):(2*n),]
    y_test <- data_2$y[(n+1):(2*n)]
    logit_test <- 1/(1 + exp(-Z_test %*% data_2$gamma))
    latent_test <- cvtLabel(logit_test, idx = T)

    NEMoE1 <- NEMoE_buildFromList(X1_train, Z_train, y_train, K = 2, lambda1 = 0.01,
                                   lambda2 = 0.01, adapt =T, btr = T,
                                   alpha1 = 0.5, alpha2 = 0.5 ,
                                   verbose= verbose, init = init)

    NEMoE1 <- fitNEMoE(NEMoE1, num_restart = 0)

    pred_logit1 <- 1/(1 + exp(-(cbind(rep(1,n),Z_test) %*% NEMoE1@NEMoE_output$gamma)))
    pred_label1 <- cvtLabel(pred_logit1, idx = T)
    ari1 = mclust::adjustedRandIndex(pred_label1, latent_test)

    NEMoE2 <- NEMoE_buildFromList(X2_train, Z_train, y_train, K = 2, lambda1 = 0.01,
                                  lambda2 = 0.01, adapt =T, btr = T,
                                  alpha1 = 0.5, alpha2 = 0.5 ,
                                  verbose= verbose, init = init)

    NEMoE2 <- fitNEMoE(NEMoE2, num_restart = 0)
    pred_logit2 <- 1/(1 + exp(-(cbind(rep(1,n),Z_test) %*% NEMoE2@NEMoE_output$gamma)))
    pred_label2 <- cvtLabel(pred_logit2, idx = T)
    ari2 = mclust::adjustedRandIndex(pred_label2, latent_test)

    NEMoE3 <- NEMoE_buildFromList(X3_train, Z_train, y_train, K = 2, lambda1 = 0.01,
                                  lambda2 = 0.01, adapt =T, btr = T,
                                  alpha1 = 0.5, alpha2 = 0.5 ,
                                  verbose= verbose, init = init)

    NEMoE3 <- fitNEMoE(NEMoE3, num_restart = 0)
    pred_logit3 <- 1/(1 + exp(-(cbind(rep(1,n),Z_test) %*% NEMoE3@NEMoE_output$gamma)))
    pred_label3 <- cvtLabel(pred_logit3, idx = T)
    ari3 = mclust::adjustedRandIndex(pred_label3, latent_test)

    NEMoE4 <- NEMoE_buildFromList(X4_train, Z_train, y_train, K = 2, lambda1 = 0.01,
                                  lambda2 = 0.01, adapt =T, btr = T,
                                  alpha1 = 0.5, alpha2 = 0.5 ,
                                  verbose= verbose, init = init)

    NEMoE4 <- fitNEMoE(NEMoE4, num_restart = 0)
    pred_logit4 <- 1/(1 + exp(-(cbind(rep(1,n),Z_test) %*% NEMoE4@NEMoE_output$gamma)))
    pred_label4 <- cvtLabel(pred_logit4, idx = T)
    ari4 = mclust::adjustedRandIndex(pred_label4, latent_test)

    NEMoE <- NEMoE_buildFromList(X_list_train, Z_train, y_train, K = 2, lambda1 = 0.01,
                                lambda2 = 0.01, adapt =T, btr = T,
                                alpha1 = 0.5, alpha2 = 0.5 ,
                                verbose= verbose, init = init)

    NEMoE <- fitNEMoE(NEMoE, num_restart = 0)
    pred_logit5 <- 1/(1 + exp(-(cbind(rep(1,n),Z_test) %*% NEMoE@NEMoE_output$gamma)))
    pred_label5 <- cvtLabel(pred_logit5, idx = T)
    ari5 = mclust::adjustedRandIndex(pred_label5, latent_test)

    ari[b,] <- c(ari1, ari2, ari3, ari4, ari5)
    print(paste0("b=", b,":", paste(ari[b,], collapse = ",")))

  }
  return(ari)

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

filter_comp <- function(X, thresh_func = var, thresh = 1e-4){

  X <- as.matrix(X)

  X_stat <- apply(X, 2, thresh_func)

  X <- X[, (X_stat > thresh)]

  return(X)
}

trans_comp <- function(X, eps = 1e-4, method = "asin", scale = T){

  X <- as.matrix(X)

  if(method == "clr"){

    X[X <= eps] <- eps

    X_trans <- t(scale(log(t(X))))

  }else if(method == "asin"){

    X_trans <- asin(sqrt(X))
  }

  if(scale){
    X_trans <- scale(X_trans)
  }

  return(X_trans)
}

evaluation_real <- function(X = Micro_scale, Z = Nutri_scale, y = Y_scale){

  X1 <- X_list[[1]]
  X2 <- X_list[[2]]
  X3 <- X_list[[3]]
  X4 <- X_list[[4]]
  X5 <- X_list[[5]]

  n <- nrow(X1)
  p1 <- ncol(X1)
  p2 <- ncol(X2)
  p3 <- ncol(X3)
  p4 <- ncol(X4)
  p5 <- ncol(X5)
  q <- ncol(Z)

  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  X3 <- as.matrix(X3)
  X4 <- as.matrix(X4)
  X5 <- as.matrix(X5)
  Z <- as.matrix(Z)

  pred <- array(0, dim = c(n, 16, 5))

  for(i in 1:n){

    X1_train <- X1[-i, ]
    X2_train <- X2[-i, ]
    X3_train <- X3[-i, ]
    X4_train <- X4[-i, ]
    X5_train <- X5[-i, ]
    X_train_list <- list(X1_train, X2_train, X3_train, X4_train, X5_train)

    Z_train <- Z[-i, ]
    y_train <- y[-i]

    X1_test <- t(as.matrix(X1[i,]))
    X2_test <- t(as.matrix(X2[i,]))
    X3_test <- t(as.matrix(X3[i,]))
    X4_test <- t(as.matrix(X4[i,]))
    X5_test <- t(as.matrix(X5[i,]))
    X_test_list <- list(X1_test, X2_test, X3_test, X4_test, X5_test)

    Z_test <- t(as.matrix(Z[i,]))
    y_test <- y[i]

    X1_comb_train <- cbind(X1_train, Z_train)
    X1_comb_test <- cbind(X1_test, Z_test)
    X2_comb_train <- cbind(X2_train, Z_train)
    X2_comb_test <- cbind(X2_test, Z_test)
    X3_comb_train <- cbind(X3_train, Z_train)
    X3_comb_test <- cbind(X3_test, Z_test)
    X4_comb_train <- cbind(X4_train, Z_train)
    X4_comb_test <- cbind(X4_test, Z_test)
    X5_comb_train <- cbind(X5_train, Z_train)
    X5_comb_test <- cbind(X5_test, Z_test)

    NEMoE_phy_NM <- fitNEMoE(NEMoE_buildFromList(X1_train, Z_train, y_train, lambda1 = 0.005, lambda2 = 0.01, verbose = F, init = "glmnet"))
    NEMoE_ord_NM <- fitNEMoE(NEMoE_buildFromList(X2_train, Z_train, y_train, lambda1 = 0.02, lambda2 = 0.01, verbose = F, init = "glmnet"))
    NEMoE_fam_NM <- fitNEMoE(NEMoE_buildFromList(X3_train, Z_train, y_train, lambda1 = 0.03, lambda2 = 0.02, verbose = F, init = "glmnet"))
    NEMoE_gen_NM <- fitNEMoE(NEMoE_buildFromList(X4_train, Z_train, y_train, lambda1 = 0.034, lambda2 = 0.02, verbose = F, init = "glmnet"))
    NEMoE_asv_NM <- fitNEMoE(NEMoE_buildFromList(X5_train, Z_train, y_train, lambda1 = 0.04, lambda2 = 0.02, verbose = F, init = "glmnet"))


    pred_MoE_NM1 <- NEMoE_predict(NEMoE_phy_NM, list(X1_test), Z_test, full = F)
    pred_MoE_NM2 <- NEMoE_predict(NEMoE_ord_NM, list(X2_test), Z_test, full = F)
    pred_MoE_NM3 <- NEMoE_predict(NEMoE_fam_NM, list(X3_test), Z_test, full = F)
    pred_MoE_NM4 <- NEMoE_predict(NEMoE_gen_NM, list(X4_test), Z_test, full = F)
    pred_MoE_NM5 <- NEMoE_predict(NEMoE_asv_NM, list(X5_test), Z_test, full = F)

    pred[i,1,1:5] <- c(pred_MoE_NM1, pred_MoE_NM2, pred_MoE_NM3, pred_MoE_NM4, pred_MoE_NM5)

    NEMoE_phy_MM <- fitNEMoE(NEMoE_buildFromList(X1_train, X1_train,  y_train, lambda1 = 0.01, lambda2 = 0.01, verbose = F, init = "glmnet"))
    NEMoE_ord_MM <- fitNEMoE(NEMoE_buildFromList(X2_train, X2_train,  y_train, lambda1 = 0.02, lambda2 = 0.01, verbose = F, init = "glmnet"))
    NEMoE_fam_MM <- fitNEMoE(NEMoE_buildFromList(X3_train, X3_train,  y_train, lambda1 = 0.03, lambda2 = 0.02, verbose = F, init = "glmnet"))
    NEMoE_gen_MM <- fitNEMoE(NEMoE_buildFromList(X4_train, X4_train,  y_train, lambda1 = 0.034, lambda2 = 0.02, verbose = F, init = "glmnet"))
    NEMoE_asv_MM <- fitNEMoE(NEMoE_buildFromList(X5_train, X5_train,  y_train, lambda1 = 0.04, lambda2 = 0.02, verbose = F, init = "glmnet"))

    pred_MoE_MM1 <- NEMoE_predict(NEMoE_phy_MM, X_new = list(X1_test), Z_new = X1_test, full = F)
    pred_MoE_MM2 <- NEMoE_predict(NEMoE_ord_MM, X_new = list(X2_test), Z_new = X2_test, full = F)
    pred_MoE_MM3 <- NEMoE_predict(NEMoE_fam_MM, X_new = list(X3_test), Z_new = X3_test, full = F)
    pred_MoE_MM4 <- NEMoE_predict(NEMoE_gen_MM, X_new = list(X4_test), Z_new = X4_test, full = F)
    pred_MoE_MM5 <- NEMoE_predict(NEMoE_asv_MM, X_new = list(X5_test), Z_new = X5_test, full = F)

    pred[i,2,1:5] <- c(pred_MoE_MM1, pred_MoE_MM2, pred_MoE_MM3, pred_MoE_MM4, pred_MoE_MM5)

    NEMoE_NN <- fitNEMoE(NEMoE_buildFromList(Z_train, Z_train, Response = y_train, lambda1 = 0.024, lambda2 = 0.025, verbose = F))

    pred_MoE_NN <- NEMoE_predict(NEMoE_NN, X_new = list(Z_test), Z_new = Z_test, full = F)
    pred[i,3,1:5] <- c(pred_MoE_NN, pred_MoE_NN, pred_MoE_NN, pred_MoE_NN, pred_MoE_NN)

    NEMoE_phy_MN <- fitNEMoE(NEMoE_buildFromList(Z_train, X1_train, y_train, lambda1 = 0.005, lambda2 = 0.01, verbose = F, init = "glmnet"))
    NEMoE_ord_MN <- fitNEMoE(NEMoE_buildFromList(Z_train, X2_train, y_train, lambda1 = 0.02, lambda2 = 0.01, verbose = F, init = "glmnet"))
    NEMoE_fam_MN <- fitNEMoE(NEMoE_buildFromList(Z_train, X3_train, y_train, lambda1 = 0.03, lambda2 = 0.02, verbose = F, init = "glmnet"))
    NEMoE_gen_MN <- fitNEMoE(NEMoE_buildFromList(Z_train, X4_train, y_train, lambda1 = 0.034, lambda2 = 0.02, verbose = F, init = "glmnet"))
    NEMoE_asv_MN <- fitNEMoE(NEMoE_buildFromList(Z_train, X5_train, y_train, lambda1 = 0.04, lambda2 = 0.02, verbose = F, init = "glmnet"))

    pred_MoE_MN1 <- NEMoE_predict(NEMoE_phy_MN, list(Z_test), X1_test, full = F)
    pred_MoE_MN2 <- NEMoE_predict(NEMoE_ord_MN, list(Z_test), X2_test, full = F)
    pred_MoE_MN3 <- NEMoE_predict(NEMoE_fam_MN, list(Z_test), X3_test, full = F)
    pred_MoE_MN4 <- NEMoE_predict(NEMoE_gen_MN, list(Z_test), X4_test, full = F)
    pred_MoE_MN5 <- NEMoE_predict(NEMoE_asv_MN, list(Z_test), X5_test, full = F)
    pred[i,4,1:5] <- c(pred_MoE_MN1, pred_MoE_MN2, pred_MoE_MN3, pred_MoE_MN4, pred_MoE_MN5)

    NEMoE_phy_ext <- fitNEMoE(NEMoE_buildFromList(X1_comb_train, X1_comb_train, y_train, lambda1 = 0.005, lambda2 = 0.01, verbose = F, init = "glmnet"))
    NEMoE_ord_ext <- fitNEMoE(NEMoE_buildFromList(X2_comb_train, X2_comb_train, y_train, lambda1 = 0.02, lambda2 = 0.01, verbose = F, init = "glmnet"))
    NEMoE_fam_ext <- fitNEMoE(NEMoE_buildFromList(X3_comb_train, X3_comb_train, y_train, lambda1 = 0.03, lambda2 = 0.02, verbose = F, init = "glmnet"))
    NEMoE_gen_ext <- fitNEMoE(NEMoE_buildFromList(X4_comb_train, X4_comb_train, y_train, lambda1 = 0.034, lambda2 = 0.02, verbose = F, init = "glmnet"))
    NEMoE_asv_ext <- fitNEMoE(NEMoE_buildFromList(X5_comb_train, X5_comb_train, y_train, lambda1 = 0.04, lambda2 = 0.02, verbose = F, init = "glmnet"))

    pred_MoE_ext1 <- NEMoE_predict(NEMoE_phy_ext, list(X1_comb_test), X1_comb_test, full = F)
    pred_MoE_ext2 <- NEMoE_predict(NEMoE_ord_ext, list(X2_comb_test), X2_comb_test, full = F)
    pred_MoE_ext3 <- NEMoE_predict(NEMoE_fam_ext, list(X3_comb_test), X3_comb_test, full = F)
    pred_MoE_ext4 <- NEMoE_predict(NEMoE_gen_ext, list(X4_comb_test), X4_comb_test, full = F)
    pred_MoE_ext5 <- NEMoE_predict(NEMoE_asv_ext, list(X5_comb_test), X5_comb_test, full = F)
    pred[i,5,1:5] <- c(pred_MoE_ext1, pred_MoE_ext2, pred_MoE_ext3, pred_MoE_ext4, pred_MoE_ext5)


    NEMoE_phy_NM_3 <- fitNEMoE(NEMoE_buildFromList(X1_train, Z_train, y_train, lambda1 = 0.008, lambda2 = 0.01, verbose = F, K = 3, alpha2 = 1, init = "glmnet"), restart_it = 0)
    NEMoE_ord_NM_3 <- fitNEMoE(NEMoE_buildFromList(X2_train, Z_train, y_train, lambda1 = 0.02, lambda2 = 0.01, verbose = F, K = 3, alpha2 = 1, init = "glmnet"), restart_it = 0)
    NEMoE_fam_NM_3 <- fitNEMoE(NEMoE_buildFromList(X3_train, Z_train, y_train, lambda1 = 0.03, lambda2 = 0.02, verbose = F, K = 3, alpha2 = 1, init = "glmnet"), restart_it = 0)
    NEMoE_gen_NM_3 <- fitNEMoE(NEMoE_buildFromList(X4_train, Z_train, y_train, lambda1 = 0.034, lambda2 = 0.02, verbose = F, K = 3, alpha2 = 1, init = "glmnet"), restart_it = 0)
    NEMoE_asv_NM_3 <- fitNEMoE(NEMoE_buildFromList(X5_train, Z_train, y_train, lambda1 = 0.04, lambda2 = 0.02, verbose = F, K = 3, alpha2 = 1, init = "glmnet"), restart_it = 0)

    pred_MoE_NM3_1 <- NEMoE_predict(NEMoE_phy_NM_3, list(X1_test), Z_test, full = F)
    pred_MoE_NM3_2 <- NEMoE_predict(NEMoE_ord_NM_3, list(X2_test), Z_test, full = F)
    pred_MoE_NM3_3 <- NEMoE_predict(NEMoE_fam_NM_3, list(X3_test), Z_test, full = F)
    pred_MoE_NM3_4 <- NEMoE_predict(NEMoE_gen_NM_3, list(X4_test), Z_test, full = F)
    pred_MoE_NM3_5 <- NEMoE_predict(NEMoE_asv_NM_3, list(X5_test), Z_test, full = F)
    pred[i,6,1:5] <- c(pred_MoE_NM3_1, pred_MoE_NM3_2, pred_MoE_NM3_3, pred_MoE_NM3_4, pred_MoE_NM3_5)

    NEMoE_phy_NM_4 <- fitNEMoE(NEMoE_buildFromList(X1_train, Z_train, y_train, lambda1 = 0.01, lambda2 = 0.01, verbose = F, K = 4, alpha2 = 1, init = "glmnet"), restart_it = 0)
    NEMoE_ord_NM_4 <- fitNEMoE(NEMoE_buildFromList(X2_train, Z_train, y_train, lambda1 = 0.02, lambda2 = 0.01, verbose = F, K = 4, alpha2 = 1, init = "glmnet"), restart_it = 0)
    NEMoE_fam_NM_4 <- fitNEMoE(NEMoE_buildFromList(X3_train, Z_train, y_train, lambda1 = 0.03, lambda2 = 0.02, verbose = F, K = 4, alpha2 = 1, init = "glmnet"), restart_it = 0)
    NEMoE_gen_NM_4 <- fitNEMoE(NEMoE_buildFromList(X4_train, Z_train, y_train, lambda1 = 0.034, lambda2 = 0.02, verbose = F, K = 4, alpha2 = 1, init = "glmnet"), restart_it = 0)
    NEMoE_asv_NM_4 <- fitNEMoE(NEMoE_buildFromList(X5_train, Z_train, y_train, lambda1 = 0.04, lambda2 = 0.02, verbose = F, K = 4, alpha2 = 1, init = "glmnet"), restart_it = 0)

    pred_MoE_NM4_1 <- NEMoE_predict(NEMoE_phy_NM_4, list(X1_test), Z_test, full = F)
    pred_MoE_NM4_2 <- NEMoE_predict(NEMoE_ord_NM_4, list(X2_test), Z_test, full = F)
    pred_MoE_NM4_3 <- NEMoE_predict(NEMoE_fam_NM_4, list(X3_test), Z_test, full = F)
    pred_MoE_NM4_4 <- NEMoE_predict(NEMoE_gen_NM_4, list(X4_test), Z_test, full = F)
    pred_MoE_NM4_5 <- NEMoE_predict(NEMoE_asv_NM_4, list(X5_test), Z_test, full = F)
    pred[i,7,1:5] <- c(pred_MoE_NM4_1, pred_MoE_NM4_2, pred_MoE_NM4_3, pred_MoE_NM4_4, pred_MoE_NM4_5)

    res_glmnet1 <- glmnet_cv(X = X1_train, y = y_train)
    res_glmnet2 <- glmnet_cv(X = X2_train, y = y_train)
    res_glmnet3 <- glmnet_cv(X = X3_train, y = y_train)
    res_glmnet4 <- glmnet_cv(X = X4_train, y = y_train)
    res_glmnet5 <- glmnet_cv(X = X5_train, y = y_train)

    pred_glmnet1 <- as.matrix(predict(res_glmnet1, newx = X1_test, type = "response"))
    pred_glmnet2 <- as.matrix(predict(res_glmnet2, newx = X2_test, type = "response"))
    pred_glmnet3 <- as.matrix(predict(res_glmnet3, newx = X3_test, type = "response"))
    pred_glmnet4 <- as.matrix(predict(res_glmnet4, newx = X4_test, type = "response"))
    pred_glmnet5 <- as.matrix(predict(res_glmnet5, newx = X5_test, type = "response"))
    pred[i,8,1:5] <- c(pred_glmnet1, pred_glmnet2, pred_glmnet3, pred_glmnet4, pred_glmnet5)

    res_glmnet2_1 <- try(two_stage(X = X1_train, Z = Z_train, y = y_train, K = 2, method = "glmnet"))
    res_glmnet2_2 <- try(two_stage(X = X2_train, Z = Z_train, y = y_train, K = 2, method = "glmnet"))
    res_glmnet2_3 <- try(two_stage(X = X3_train, Z = Z_train, y = y_train, K = 2, method = "glmnet"))
    res_glmnet2_4 <- try(two_stage(X = X4_train, Z = Z_train, y = y_train, K = 2, method = "glmnet"))
    res_glmnet2_5 <- try(two_stage(X = X5_train, Z = Z_train, y = y_train, K = 2, method = "glmnet"))
    pred_glmnet2_1 <- as.matrix(predict_two_stage(X_new = X1_test, Z_new = Z_test, y_new = y_test, twostage = res_glmnet2_1, method = "glmnet", output = "reponse"))
    pred_glmnet2_2 <- as.matrix(predict_two_stage(X_new = X2_test, Z_new = Z_test, y_new = y_test, twostage = res_glmnet2_2, method = "glmnet", output = "reponse"))
    pred_glmnet2_3 <- as.matrix(predict_two_stage(X_new = X3_test, Z_new = Z_test, y_new = y_test, twostage = res_glmnet2_3, method = "glmnet", output = "reponse"))
    pred_glmnet2_4 <- as.matrix(predict_two_stage(X_new = X4_test, Z_new = Z_test, y_new = y_test, twostage = res_glmnet2_4, method = "glmnet", output = "reponse"))
    pred_glmnet2_5 <- as.matrix(predict_two_stage(X_new = X5_test, Z_new = Z_test, y_new = y_test, twostage = res_glmnet2_5, method = "glmnet", output = "reponse"))
    pred[i,9,1:5] <- c(pred_glmnet2_1, pred_glmnet2_2, pred_glmnet2_3, pred_glmnet2_4, pred_glmnet2_5)

    res_svm1 <- svm(x = X1_train, y = as.factor(y_train), kernel = "radial", probability = T )
    res_svm2 <- svm(x = X2_train, y = as.factor(y_train), kernel = "radial", probability = T )
    res_svm3 <- svm(x = X3_train, y = as.factor(y_train), kernel = "radial", probability = T )
    res_svm4 <- svm(x = X4_train, y = as.factor(y_train), kernel = "radial", probability = T )
    res_svm5 <- svm(x = X5_train, y = as.factor(y_train), kernel = "radial", probability = T )

    pred_svm1 <- as.matrix(attr(predict(res_svm1, X1_test, probability = T), "probabilities")[,"1"])
    pred_svm2 <- as.matrix(attr(predict(res_svm2, X2_test, probability = T), "probabilities")[,"1"])
    pred_svm3 <- as.matrix(attr(predict(res_svm3, X3_test, probability = T), "probabilities")[,"1"])
    pred_svm4 <- as.matrix(attr(predict(res_svm4, X4_test, probability = T), "probabilities")[,"1"])
    pred_svm5 <- as.matrix(attr(predict(res_svm5, X5_test, probability = T), "probabilities")[,"1"])
    pred[i,10,1:5] <- c(pred_svm1, pred_svm2, pred_svm3, pred_svm4, pred_svm5)


    res_svm2_1 <- try(two_stage(X = X1_train, Z = Z_train, y = y_train, K = 2, method = "svm"))
    res_svm2_2 <- try(two_stage(X = X2_train, Z = Z_train, y = y_train, K = 2, method = "svm"))
    res_svm2_3 <- try(two_stage(X = X3_train, Z = Z_train, y = y_train, K = 2, method = "svm"))
    res_svm2_4 <- try(two_stage(X = X4_train, Z = Z_train, y = y_train, K = 2, method = "svm"))
    res_svm2_5 <- try(two_stage(X = X5_train, Z = Z_train, y = y_train, K = 2, method = "svm"))

    pred_svm2_1 <- as.matrix(predict_two_stage(X_new = X1_test, Z_new = Z_test, y_new = y_test, twostage = res_svm2_1, method = "svm", output = "reponse"))
    pred_svm2_2 <- as.matrix(predict_two_stage(X_new = X2_test, Z_new = Z_test, y_new = y_test, twostage = res_svm2_2, method = "svm", output = "reponse"))
    pred_svm2_3 <- as.matrix(predict_two_stage(X_new = X3_test, Z_new = Z_test, y_new = y_test, twostage = res_svm2_3, method = "svm", output = "reponse"))
    pred_svm2_4 <- as.matrix(predict_two_stage(X_new = X4_test, Z_new = Z_test, y_new = y_test, twostage = res_svm2_4, method = "svm", output = "reponse"))
    pred_svm2_5 <- as.matrix(predict_two_stage(X_new = X5_test, Z_new = Z_test, y_new = y_test, twostage = res_svm2_5, method = "svm", output = "reponse"))
    pred[i,11,1:5] <- c(pred_svm2_1, pred_svm2_2, pred_svm2_3, pred_svm2_4, pred_svm2_5)

    res_rf1 <- randomForest(x = X1_train, y = as.factor(y_train), ntree = 500)
    res_rf2 <- randomForest(x = X2_train, y = as.factor(y_train), ntree = 500)
    res_rf3 <- randomForest(x = X3_train, y = as.factor(y_train), ntree = 500)
    res_rf4 <- randomForest(x = X4_train, y = as.factor(y_train), ntree = 500)
    res_rf5 <- randomForest(x = X5_train, y = as.factor(y_train), ntree = 500)

    pred_rf1 <- as.matrix(predict(res_rf1, X1_test, type = "prob")[,"1"])
    pred_rf2 <- as.matrix(predict(res_rf2, X2_test, type = "prob")[,"1"])
    pred_rf3 <- as.matrix(predict(res_rf3, X3_test, type = "prob")[,"1"])
    pred_rf4 <- as.matrix(predict(res_rf4, X4_test, type = "prob")[,"1"])
    pred_rf5 <- as.matrix(predict(res_rf5, X5_test, type = "prob")[,"1"])
    pred[i,12,1:5] <- c(pred_rf1, pred_rf2, pred_rf3, pred_rf4, pred_rf5)

    res_rf2_1 <- try(two_stage(X = X1_train, Z = Z_train, y = y_train, K = 2, method = "random forest"))
    res_rf2_2 <- try(two_stage(X = X2_train, Z = Z_train, y = y_train, K = 2, method = "random forest"))
    res_rf2_3 <- try(two_stage(X = X3_train, Z = Z_train, y = y_train, K = 2, method = "random forest"))
    res_rf2_4 <- try(two_stage(X = X4_train, Z = Z_train, y = y_train, K = 2, method = "random forest"))
    res_rf2_5 <- try(two_stage(X = X5_train, Z = Z_train, y = y_train, K = 2, method = "random forest"))
    pred_rf2_1 <- as.matrix(predict_two_stage(X_new = X1_test, Z_new = Z_test, y_new = y_test, twostage = res_rf2_1, method = "random forest", output = "reponse"))
    pred_rf2_2 <- as.matrix(predict_two_stage(X_new = X2_test, Z_new = Z_test, y_new = y_test, twostage = res_rf2_2, method = "random forest", output = "reponse"))
    pred_rf2_3 <- as.matrix(predict_two_stage(X_new = X3_test, Z_new = Z_test, y_new = y_test, twostage = res_rf2_3, method = "random forest", output = "reponse"))
    pred_rf2_4 <- as.matrix(predict_two_stage(X_new = X4_test, Z_new = Z_test, y_new = y_test, twostage = res_rf2_4, method = "random forest", output = "reponse"))
    pred_rf2_5 <- as.matrix(predict_two_stage(X_new = X5_test, Z_new = Z_test, y_new = y_test, twostage = res_rf2_5, method = "random forest", output = "reponse"))
    pred[i,13,1:5] <- c(pred_rf2_1, pred_rf2_2, pred_rf2_3, pred_rf2_4, pred_rf2_5)

    NEMoE_l <- fitNEMoE(NEMoE_buildFromList(Microbiome = X_train_list, Nutrition = Z_train, Response = y_train,
                                 lambda1 = c(0.005, 0.012, 0.013, 0.023, 0.025),
                                 lambda2 = 0.02, alpha1 = 0.5, alpha2 = 0.5,
                                 cvParams = createCVList(g1 = 10, shrink = 0.4,
                                                         track = F), itmax = 1e3))

    pred_res_l <- NEMoE_predict(X_new = X_test_list, Z_new = Z_test, NEMoE = NEMoE_l, full = T)$output
    pred[i,14,1:5] <- c(pred_res_l[1], pred_res_l[2], pred_res_l[3], pred_res_l[4], pred_res_l[4])

    pred_res_l <- NEMoE_predict(X_new = X_test_list, Z_new = Z_test, NEMoE = NEMoE, full = T)$output
    pred[i,15,1:5] <- c(pred_res_l[1], pred_res_l[2], pred_res_l[3], pred_res_l[4], pred_res_l[4])
    pred[i,16,1:5] <- rep(y_test,5)

    print(paste0("i(Family):",i, ";", pred[i,,3]))
    print(paste0("i(Genus):",i, ";", pred[i,,4]))

  }

  return(pred)
}


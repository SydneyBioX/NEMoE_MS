## Preprocessing

library(phyloseq)
## Load ASV table
load("E:\\data\\Precision Nutrition\\data\\PD\\PD_silva\\ourdata_silva_with_tree.RData")
otu_tab0 = as.matrix(res$otu_tab)
tax_tab0 = res$taxTab
tree = res$tree

## Split sample names

sampleNames0 <- c()
for(i in 1:nrow(otu_tab0)){
  name_temp = rownames(otu_tab0)[i]
  sampleNames0[i] = strsplit(name_temp, "_")[[1]][1]
}
rownames(otu_tab0) = sampleNames0

sample_sp = strsplit(sampleNames0,"-t-")
sample_info0 = data.frame(id = 0,t = 0)
for (i in 1:length(sample_sp)){
  if(length(sample_sp[[i]]) == 2){
    id = sample_sp[[i]][1]
    id = sub(pattern = "-",replacement = "/",x = id)
    sample_time = sample_sp[[i]][2]
    if (sample_time == "minus2"){
      sample_time = -2
    }else if(is.na(as.numeric(sample_time))){
      sample_time = 0
    }else{
      sample_time = as.numeric(sample_time)
    }

  }else{
    id = sample_sp[[i]][1]
    sample_time = NA
  }
  sample_info0 = rbind(sample_info0,c(id,sample_time))
}
sample_info0 = sample_info0[-1,]
rownames(sample_info0) = sampleNames0
rownames(otu_tab0) = sampleNames0

## Load meta data

Clinical1 <- readxl::read_xlsx("D:\\Rcode\\microbiome\\old\\stage3\\Project Database - All Time Intervals.xlsx",sheet = 1)
var_name = c("PD ID", "Age","Sex ID","Cohort Identifier ID","Height (cm)","Weight (kg)","BMI","Ethnicity ID", "Marital Status ID","Education ID","Employment ID","Income ID","Health Insurance ID","Support Services ID","Allied Health in Last Yr ID","Physio ID", "PD Duration","PD Phenotype ID","Late Onset >60 yrs ID","Genetic ID", "Rehabilitation ID", "Dyskinesia ID","On-Off Fluctuation ID","Wearing off ID","Anosmia ID", "Neuroleptic Use ID","PD FHx ID","Pesticide ID", "Current Smoker ID","Prior Smoker ID","Never Smoked ID" ,"No. Cups (day)","EtOH Consumption ID","EtOH < 1 weekly ID","EtOH 1-6 days week ID","EtOH Daily ID", "Diabetes ID","DM Duration (years)","HbA1c...111", "Medication PD Rx ID","L-Dopa ID","Duodopa ID","DBS ID","Apomorphine ID","Dopamine Agonist ID","MAOB ID","Anticholinergic ID","COMT ID","Amantidine ID", "Physical Functioning","Role Functioning_Physical","Energy_Fatigue","Emotional Well-being","Social Functioning","General Health","Health Change","Physical Component Summary"  ,"Pain","Chronic Pain ID","Pain Score","Walk 1 km ID","Climb 1 flight stairs ID","IPAQ MET-min/week","IPAQ Sitting hr/day","IPAQ Category Score ID", "Beck's Depression Inventory","Depressed ID","MOCA Visuospacial","MOCA Naming","MOCA Attention","MOCA Language","MOCA Abstraction","MOCA Delayed Recall","MOCA Orientation","MOCA total score","Mild Cognitive Impairment ID","PD Dementia ID", "Constipation ID...74","Constipation ID...170","Cleveland Constipation Score","Strinaing ID","Hard Stools ID","Incomplete Evacucation ID","Manual Manoeuvres ID","Spontaneous bowel mvts / week ID","Rome IV Total Score","LDQ Score","Most troublesome symptom ID","Bristol Score", "NMSS Cardiovascular","NMSS Sleep and fatigue","NMSS Mood and cognition","NMSS Perceptual problems","NMSS Attention and Memory","NMSS Gastrointestinal","NMSS Urinary","NMSS Sexual","NMSS Miscellaneous","NMSS Total Score", "PDQ Mobility","PDQ ADL","PDQ Emotions","PDQ Stigma","PDQ Social","PDQ Cognitions","PDQ Communication","PDQ Body Pain","PDQ Summary Index", "UPDRS Speech","UPDRS Facial Expression","UPDRS Rigidity","UPDRS Finger Tapping","UPDRS Hand Movements","UPDRS Pronation Supination","UPDRS Toe Tapping","UPDRS Leg Agility","UPDRS Rising from Chair","UPDRS Gait","UPDRS Freezing","UPDRS Postural Stability","UPDRS Posture","UPDRS Body Bradykinesia","UPDRS Postural Hand Tremor","UPDRS Kinetic Hand Tremor","UPDRS Rest Tremor","UPDRS Consistancy of Rest Tremor","UPDRS Total Score","Hoehn and Yahr Stage","ESR","CRP","Total Cholesterol", "LDL", "HDL", "Glucose", "Trigs", "HbA1c...262", "Albumin","ICD ID","RBD ID","OT ID","Speech Pathologist ID","Dietician ID","Pack YHx","Type of Tabaco ID","Last Abx use (months)","Telemedicine Interest ID","Distance to clinic (km)","Total LED (mg)","Head Trauma ID","Vegetarian Diet ID","Tot_Energy_with_Dietary_Fibre_kJ_per_day","Tot_Energy_without_Dietary_Fibre_kJ_per_day","Tot_Moisture","Tot_Protein_g_per_day","Tot_Fat_g_per_day","Tot_Fibre_g_per_day","Tot_Alcohol_g_per_day","Tot_of_Total_Sugars_g_per_day","Tot_Free_Sugars_g_per_day","Tot_Added_Sugars_g_per_day","Tot_Avail_Carbs_g_per_day","Tot_Calcium_mg_per_day","Tot_Iron_mg_per_day","Tot_Magn_mg_per_day","Tot_Potas_mg_per_day","Tot_Sodium_mg_per_day","Tot_Zinc_mg_per_day","Tot_Retinol_ug_per_day","Tot_Beta_car_ug_per_day","Tot_Vit_A_RE_ug_per_day","Tot_Thiamin_mg_per_day","Tot_Ribo_mg_per_day","Tot_B12_ug_per_day","Tot_Dietary_Folate_ug_per_day","Tot_Vit_C_mg_per_day")
Clinical1 <- Clinical1[!is.na(Clinical1$`Neurogenetics ID`),]

## Zymo and replicates

sample_info_micro <- sample_info0
sample_info_micro$type <- "normal"

sampleNames_micro <- rownames(sample_info_micro)

flag1 <- substr(sampleNames_micro,1,4)
flag2 <- substr(sampleNames_micro,12,14)

sample_info_micro$type[flag1 == "Zymo"] = "Zymo"
sample_info_micro$t[flag1 == "Zymo"] = 0
sample_info_micro$type[flag2 == "rpt"] = "rpt"

## Split time point

Clinical1$t = Clinical1$Intervals
Clinical1$t[Clinical1$t == 1] = -2
Clinical1$t[Clinical1$t == 2] = 0
Clinical1$t[Clinical1$t == 3] = 2
Clinical1$t[Clinical1$t == 6] = 12
Clinical1$t[Clinical1$t == 5] = 6
idx <- c()
for(i in 1:nrow(sample_info_micro)){
  if(sample_info_micro$type[i] == "Zymo"){
    idx[i] <- NA
  }else{
    temp1 <- which(Clinical1$`Neurogenetics ID` == sample_info_micro[i,1])
    if(length(temp1) == 0){
      idx[i] <- NA
    }else if(length(temp1) == 1){
      idx[i] = temp1
    }else{
      temp2 <- which(Clinical1$t[temp1] == sample_info_micro[i,2])
      if(length(temp2)){
        idx[i] = temp1[temp2]
      }
      else{
        idx[i] = temp1[1]
      }
    }
  }

}
sample_info = cbind(sample_info0, Clinical1[idx,var_name])
colnames(sample_info)[3] = "PD"
colnames(sample_info)[6] = "Treatment"

sampleNames <- rownames(otu_tab0)
meta_data <-as.data.frame(sample_info)
rownames(meta_data) = sampleNames

## Create phyloseq object

ps <- phyloseq(otu_table(otu_tab0, taxa_are_rows = F), sample_data(meta_data), tax_table(tax_tab0))
keep_sample <- !grepl("Zymo|rpt",sample_names(ps))
ps <- subset_samples(ps, keep_sample)

## Filter out ASV with low prevalence
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

#colnames(prevdf) = c("Prevalence","Total Abundance")

prevdf_phylum = plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

colnames(prevdf_phylum) = c("Phylum","Mean of Prevalence", "Sum of Prevalence")

filterPhyla = c("Planctomycetota ","Elusimicrobiota","Campylobacterota","Cyanobacteria", "Fusobacteriota", "Thermoplasmatota")

# Filter entries with unidentified Phylum.
ps1 = subset_taxa(ps, !Phylum %in% filterPhyla)

prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))

#ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
#  geom_hline(yintercept = 0.1, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
#  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
#  facet_wrap(~Phylum) + theme(legend.position="none") + theme_bw()

prevalenceThreshold = floor(0.05 * nsamples(ps))
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps_filt = prune_taxa(keepTaxa, ps)
temp_dim = dim(ps_filt@otu_table)

ps0 <- subset_samples(ps_filt, (t == 0))

Z_nutri <- ps0@sam_data[,c(153:177)]
Z_nutri$Pratio <- unlist(Z_nutri[,"Tot_Protein_g_per_day"] * 17 / Z_nutri[,"Tot_Energy_without_Dietary_Fibre_kJ_per_day"])
Z_nutri$Fratio <- unlist(Z_nutri[,"Tot_Fat_g_per_day"] * 37 / Z_nutri[,"Tot_Energy_without_Dietary_Fibre_kJ_per_day"])
Z_nutri$Cratio <- unlist(Z_nutri[,"Tot_Avail_Carbs_g_per_day"] * 17 / Z_nutri[,"Tot_Energy_without_Dietary_Fibre_kJ_per_day"])
Z_nutri$PC <- unlist(Z_nutri[,"Tot_Protein_g_per_day"] / Z_nutri[,"Tot_Avail_Carbs_g_per_day"])

Z_nutri1 <- na.omit(Z_nutri)
ps1 <- subset_samples(ps0, sample_names(ps0) %in% rownames(Z_nutri1))

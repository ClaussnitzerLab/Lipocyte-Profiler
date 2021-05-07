######LP profile  of tail ends of polygenic risk#####
##set working directory###
#setwd("/Users/sophiestrobel/Desktop/LP/data raw/LP data import")
####input data you want to analyse####
#combined_data<-read.csv(file = 'LP_data.csv')
#combined_data<-combined_data[,-1]

######subset data set#####
#all_data_nonFFA<-subset(combined_data, combined_data$FFA == "0")
all_data_FFA<-subset(combined_data, combined_data$FFA == "1")
# 
# data_nonFFA_day0_vc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "vc" & all_data_nonFFA$day == "0")
# data_nonFFA_day3_vc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "vc" & all_data_nonFFA$day == "3")
# data_nonFFA_day8_vc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "vc" & all_data_nonFFA$day == "8")
# data_nonFFA_day14_vc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "vc" & all_data_nonFFA$day == "14")
# 
# data_nonFFA_day0_sc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "sc" & all_data_nonFFA$day == "0")
# data_nonFFA_day3_sc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "sc" & all_data_nonFFA$day == "3")
# data_nonFFA_day8_sc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "sc" & all_data_nonFFA$day == "8")
# data_nonFFA_day14_sc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "sc" & all_data_nonFFA$day == "14")
# 
# data_FFA_day0_vc<-subset(all_data_FFA, all_data_FFA$cellType == "vc" & all_data_FFA$day == "0")
# data_FFA_day3_vc<-subset(all_data_FFA, all_data_FFA$cellType == "vc" & all_data_FFA$day == "3")
#data_FFA_day8_vc<-subset(all_data_FFA, all_data_FFA$cellType == "vc" & all_data_FFA$day == "8")
#data_FFA_day14_vc<-subset(all_data_FFA, all_data_FFA$cellType == "vc" & all_data_FFA$day == "14")
# 
# data_FFA_day0_sc<-subset(all_data_FFA, all_data_FFA$cellType == "sc" & all_data_FFA$day == "0")
# data_FFA_day3_sc<-subset(all_data_FFA, all_data_FFA$cellType == "sc" & all_data_FFA$day == "3")
# data_FFA_day8_sc<-subset(all_data_FFA, all_data_FFA$cellType == "sc" & all_data_FFA$day == "8")
# data_FFA_day14_sc<-subset(all_data_FFA, all_data_FFA$cellType == "sc" & all_data_FFA$day == "14")

##### use same n numbers for all comparison, adj. n to n of D8####
# AP_subset<-data_FFA_day8_vc
# patient<-AP_subset$patientID
#############input data you want to analyse##########
AP_subset<-data_FFA_day14_vc
######adj n number###
#AP_subset<-subset(AP_subset, AP_subset$patientID %in% patient)
#####subset in meta and LP feature data######
meta<-AP_subset[, 1:36]
input<-AP_subset[, 37:3038]
#####filter data, remove blocklisted features, SmallBODIPY objects and features with 0´s or NA´s across all samples#####
#####see methods LP manuscripte Laber and Strobel et al.#####
variables<-input

variables[,] <- sapply(sapply(variables[,], as.numeric), as.numeric) 
back<-as.data.frame(t(variables))

var.selection.1<-back[!grepl( "_Costes_" , rownames(back) ) ,  ]
var.selection.2<-var.selection.1[!grepl( "_Manders_" , rownames(var.selection.1)) ,  ]
var.selection.3<-var.selection.2[!grepl( "_RWC_" , rownames(var.selection.2)) ,  ]
var.selection.4<-var.selection.3[!grepl( "SmallBODIPY" , rownames(var.selection.3) ) ,  ]

mean<-rowMeans(var.selection.4)
var<-cbind(mean, var.selection.4)

var.1<-subset(var, var$mean != "0")
var.1a<- var.1[complete.cases(var.1), ]

variables<-as.data.frame(t(var.1a[,2:36]))

#####merge meta and LP data###
df<-cbind(meta, variables)

####categories PRS percentiles####
####high risk > 75% -> 2
####low risk < 25% -> 1
library(dplyr)

PRSm <- df %>%
  dplyr::mutate(PRSHOMA = ifelse(
    df$HOMAPERC>= 75,"2",
    ifelse (df$HOMAPERC<= 25, "0", "1")))

PRSmm <- PRSm %>%
  dplyr::mutate(PRSWHR = ifelse(
    PRSm$PERCWHRc   >= 75	, "2",
    ifelse (PRSm$PERCWHRc <= 25	, "0", "1")))

PRSmmm <- PRSmm %>%
  dplyr::mutate(PRSWHRf = ifelse(
    PRSmm$PERCWHRf   >= 75	, "2",
    ifelse (PRSmm$PERCWHRf <= 25	, "0", "1")))

PRSmmmm <- PRSmmm %>%
  dplyr::mutate(PRST2D = ifelse(
    PRSmmm$T2DPERC   >= 75	, "2",
    ifelse (PRSmmm$T2DPERC <= 25	, "0", "1")))

####consider just feamles for WHRadjBMI PRS#####
#PRSmmmm_f<-subset(PRSmmmm, PRSmmmm$sex == "2")
#########################################################
##########t-test##########
#########case = top########
#########control = bottom #########

case<-subset(PRSmmmm, PRSmmmm$PRSHOMA == "2")

control<-subset(PRSmmmm, PRSmmmm$PRSHOMA  == "0")

############build matrix for anova
data_anov<-rbind(case, control)

#####as numeric LP features data, change range depending on input#####
data_anov[,37:2723] <- sapply(sapply(data_anov[,37:2723], as.character),as.numeric) 

#######ANOVA######
colvector<-as.character(colnames(variables))

library(sjstats)

lm_model<-list()

for(i in colvector){
  
  lm <- lm(data_anov[,i] ~ PRSHOMA  
                          + age
                          + BMI
                          + sex     #add sex as variable if your dataset includes both feamle and male
                          + batch, data = data_anov)
  
  lm_model[[i]]<-anova(lm)
}

#####extract p-values #####
pvalue<-list()

for(i in colvector){
  
  pvalue[[i]]<-lm_model[[i]][["Pr(>F)"]]
  
}
###structure pvalue dataframe#####
pvalue_df<-do.call(cbind, pvalue)
pvalue_df_1<-as.data.frame(t(pvalue_df))
pvalue_df_1<-pvalue_df_1[, -6] #adj depending on number of variables
colnames(pvalue_df_1)<-c("PRSHOMA.pvalue", 
                         "age.pvalue", 
                         "BMI.pvalue", 
                         "sex.pvalue", #add sex as variable if your dataset includes both feamle and male
                         "batch.pvalue")

###calculate q-values (FDR)#####
library(qvalue)

p<-pvalue_df_1$PRSHOMA.pvalue

qvalues<-qvalue(p)
hist(qvalues)

q<-qvalues[["qvalues"]]

##combine q-values and p-values#####
pvalue_data<-cbind(pvalue_df_1, q)

###extract F value#####
Fvalue<-list()
for(i in colvector){
  
  Fvalue[[i]]<-lm_model[[i]][["F value"]]
  }

Fvalue_df<-do.call(cbind, Fvalue)
Fvalue_df_1<-as.data.frame(t(Fvalue_df))
Fvalue_df_1<-Fvalue_df_1[, -6] #adj depending on number of variables
colnames(Fvalue_df_1)<-c("PRSHOMA.Fvalue", 
                         "age.F", 
                         "BMI.F", 
                         "sex.F" ,#add sex as variable if your dataset includes both feamle and male
                         "batch.F")

####extract effect size (eta_sq)####
stats<-list()

for(i in colvector){
  
  stats[[i]]<-anova_stats(lm_model[[i]])
  
}

anova_stats(lm_model[[3]])

effect<-sapply(stats,function(x) x[7])
effect_df<-do.call(cbind, effect)
effect_df_1<-as.data.frame(t(effect_df))
effect_df_1<-effect_df_1[, -6] #adj depending on number of variables
colnames(effect_df_1)<-c("PRSHOMA.eta_sq",   
                         "age.eta_sq", 
                         "BMI.eta_sq", 
                        "sex.eta_sq", #add sex as variable if your dataset includes both feamle and male
                         "batch.eta_sq")
features<-colvector

####t-test to calculate t-statistic####
#####as numeric and as.matrix case and control dataframe####
mcluster<-control[,37:2723]

mcluster[,] <- sapply(sapply(mcluster[,], as.character),as.numeric) 

matrixcontrol<-as.matrix(mcluster)
#######numeric matrix case 
mcluster<-case[,37:2723]

mcluster[,] <- sapply(sapply(mcluster[,], as.character),as.numeric) 

matrixcase<-as.matrix(mcluster)

###t-test###
outtest<- list ()

for(i in colvector){
  
  x<-matrixcase[,i]
  y<-matrixcontrol[,i]
  
  outtest[[i]]<-t.test(x,y)
}

pvalue_ttest <- data.frame(matrix(unlist(outtest), nrow=length(outtest), byrow=T))

rownames(pvalue_ttest)<-colvector
colnames(pvalue_ttest)<-c("t-test", "df", "pvalue", "confin1", "confin2", "meanx", "meany", "diffmean", "stderr", "alternative", "method", "dataname")

###calculate qvalue, FDR####
library(qvalue)
p<- as.numeric(as.character(pvalue_ttest$pvalue))

qvalues<-qvalue(p)
hist(qvalues)

q<-qvalues[["qvalues"]]

####combine q-value, p-value and t-statistics####
pFDR<-cbind(p, q, pvalue_ttest)

ttest_pFDR<-as.data.frame(t(pFDR))
#####structure dataframes####
ttest_pFDR<-rbind(ttest_pFDR["t-test",],ttest_pFDR["p",],ttest_pFDR["q",] )

ttest_pFDR<-as.data.frame(t(ttest_pFDR))
colnames(ttest_pFDR)<-c("t-test", "p","qvalue")

ttest_pFDR[,] <- sapply(sapply(ttest_pFDR[,], as.character),as.numeric) 

###combine ANOVA an ttest in one output file####

volcano<-cbind(pvalue_data, effect_df_1, Fvalue_df_1 , ttest_pFDR, features )


####define color categories###
volcano <- volcano %>%
  dplyr::mutate(feature_color = 
                  ifelse(grepl("BODIPY", rownames(volcano)),'Lipid', 
                         ifelse(grepl("AGP", rownames(volcano)),'AGP', 
                                ifelse(grepl("Mito", rownames(volcano)),'Mito',
                                       ifelse(grepl("DNA", rownames(volcano)),'DNA',
                                              'other')))))

volcano$feature_color<-factor(volcano$feature_color,
                              levels = c("Mito","AGP","Lipid","DNA", "other"))
####remove missing values###
volcano<- volcano[complete.cases(volcano), ]

####set significance level###
volcano <- volcano %>%
  dplyr::mutate(feature_sig = 
                  ifelse(volcano$q <= 0.05 & volcano$PRSHOMA.pvalue <= 0.05,'5%FDR','not'))


######plot####
library(ggplot2)
plot <- ggplot(volcano, aes(x= `t-test` , y=-log10(PRSHOMA.pvalue))) +
  geom_point(aes(color=feature_color,
                 size = feature_sig,
                 alpha = feature_sig ))+
  
              xlab("t-statistics") +
              ylab("-log10 p-value") +
              ylim(0,5)+
              xlim(-6,6)+
              scale_size_manual(name = "", 
                    values = c("5%FDR" = 1.5, "not" = 1))+
              scale_alpha_manual(name = "", 
                     values = c("5%FDR" = 0.75, "not" = 0.1))+
              scale_color_manual(name = "", 
                     values = c("Mito" = "#f56464" ,
                                "AGP" = "#f2cf41",
                                "Lipid" = "#4dac26",
                                "DNA" = "#65a6db",
                                "other" = "#959ca3"))+
              theme_bw() +
              theme(strip.background = element_rect(colour = "black",
                                        fill = "#fdfff4"))


pdf("PRS HOMA-IR visceral D14.pdf")
plot
dev.off()

write.table(pvalue_data, file = "PRS-HOMA profile visceral D14.csv", sep = ",", col.names = NA,
            qmethod = "double")

###########
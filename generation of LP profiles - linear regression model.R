#########linear regression of LP profile and PRS scores####
####set working directory#####
setwd("/Users/sophiestrobel/Desktop/LP/data raw/LP data import")

####load process specifc PRS scores cluster 1-3 ####
# cluster<-read.csv(file = "Miriam cluster_combined5_MOBB.csv")
# # 
# names(cluster)[names(cluster) == "FID"] <- "patientID"

###load LP data set#####
combined_data<-read.csv(file = "LP_data.csv.csv")


###merge data sets#####
#combined_cluster <- merge(x = cluster, y = combined_data, by = "patientID",  all.y= T, sort = FALSE)

####subset data set####
all_data_nonFFA<-subset(combined_data, combined_data$FFA == "0")
#all_data_FFA<-subset(combined_data, combined_data$FFA == "1")

# data_nonFFA_day0_vc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "vc" & all_data_nonFFA$day == "0")
# data_nonFFA_day3_vc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "vc" & all_data_nonFFA$day == "3")
# data_nonFFA_day8_vc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "vc" & all_data_nonFFA$day == "8")
# data_nonFFA_day14_vc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "vc" & all_data_nonFFA$day == "14")

# data_nonFFA_day0_sc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "sc" & all_data_nonFFA$day == "0")
# data_nonFFA_day3_sc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "sc" & all_data_nonFFA$day == "3")
data_nonFFA_day8_sc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "sc" & all_data_nonFFA$day == "8")
data_nonFFA_day14_sc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "sc" & all_data_nonFFA$day == "14")

# data_FFA_day0_vc<-subset(all_data_FFA, all_data_FFA$cellType == "vc" & all_data_FFA$day == "0")
# data_FFA_day3_vc<-subset(all_data_FFA, all_data_FFA$cellType == "vc" & all_data_FFA$day == "3")
# data_FFA_day8_vc<-subset(all_data_FFA, all_data_FFA$cellType == "vc" & all_data_FFA$day == "8")
#data_FFA_day14_vc<-subset(all_data_FFA, all_data_FFA$cellType == "vc" & all_data_FFA$day == "14")

# data_FFA_day0_sc<-subset(all_data_FFA, all_data_FFA$cellType == "sc" & all_data_FFA$day == "0")
# data_FFA_day3_sc<-subset(all_data_FFA, all_data_FFA$cellType == "sc" & all_data_FFA$day == "3")
# data_FFA_day8_sc<-subset(all_data_FFA, all_data_FFA$cellType == "sc" & all_data_FFA$day == "8")
# data_FFA_day14_sc<-subset(all_data_FFA, all_data_FFA$cellType == "sc" & all_data_FFA$day == "14")

##### use same n numbers for all comparison, adj. n to n of D8####
# AP_subset<-data_nonFFA_day8_sc
# patient<-AP_subset$patientID
#############input data you want to analyse##########
AP_subset<-data_nonFFA_day14_sc
######adj n number###
AP_subset<-subset(AP_subset, AP_subset$patientID %in% patient)
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
var.selection.5<-var.selection.3[!grepl( "SmallBODIPY" , rownames(var.selection.3) ) ,  ]

mean<-rowMeans(var.selection.5)
var<-cbind(mean, var.selection.5)

var.1<-subset(var, var$mean != "0")
var.1a<- var.1[complete.cases(var.1), ]

variables<-as.data.frame(t(var.1a[,2:25]))

#####merge meta and LP data###
df<-cbind(meta, variables)

.#####lm model####
####define LP features you want to test####
####default: all LP features ####
colvector<-as.character(colnames(variables))

library(sjstats)

lm_model<-list()
lm_anov<-list()

for(i in colvector){
  
  lm_model[[i]] <- summary.lm(lm(df[,i] ~ score.cluster4 + age+ BMI+ batch +sex , data = df))
  lm_anov[[i]]<-anova(lm(df[,i] ~ score.cluster4 + age+ BMI+ batch +sex, data = df))
  }
####extract p-values of linear regression####
pvalue.l<-list()

for(i in colvector){
  
  pvalue.l[[i]]<-lm_model[[i]][["coefficients"]]
  
}
###structure pvalue dataframe####
coef<-as.data.frame(sapply(pvalue.l,function(x) x[2]))
colnames(coef)<-"estimate"
pvalue_df<-sapply(pvalue.l,function(x) x[,4])
pvalue_df_1<-as.data.frame(t(pvalue_df))
colnames(pvalue_df_1)<-c("intersect", "score.pvalue", "age.pvalue", "BMI.pvalue", "batch.pvalue", "sex.pvalue")

###calculate q-values (FDR) of linear regression####
library(qvalue)

p<-pvalue_df_1$score.pvalue

qvalues<-qvalue(p)
hist(qvalues)

q.l<-qvalues[["qvalues"]]

###extract p-values of ANOVA####
pvalue.a<-list()

for(i in colvector){
  
  pvalue.a[[i]]<-lm_anov[[i]][["Pr(>F)"]]
  
}
#####structure pvalue dataframe########
pvalue_df.a<-do.call(cbind, pvalue.a)
pvalue_df.a_1<-as.data.frame(t(pvalue_df.a))
pvalue_df.a_1<-pvalue_df.a_1[, -6]
colnames(pvalue_df.a_1)<-c("score.anova.pvalue", "age.anova.pvalue", "BMI.anova.pvalue", "batch.anova.pvalue","sex.anova.pvalue")

####calculate q-value for ANOVA p-value###
p<-pvalue_df.a_1$score.anova.pvalue

qvalues<-qvalue(p)
hist(qvalues)

q.a<-qvalues[["qvalues"]]

###merge dataframes####
features<-colvector
pvalue_data<-cbind(pvalue_df_1, q.l, pvalue_df.a_1, q.a,coef, features )
####save output####
write.table(pvalue_data, file = "PRS lipodystrophy profile subq D14.csv", sep = ",", col.names = NA,
            qmethod = "double")

####plot results as volcano plot###

volcano<-pvalue_data

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
                  ifelse(volcano$q.a <= 0.1 & volcano$score.anova.pvalue <= 0.05,'5%FDR','not'))


######plot####
library(ggplot2)
plot <- ggplot(volcano, aes(x= estimate , y=-log10(score.anova.pvalue))) +
          geom_point(aes(color=feature_color,
                          size = feature_sig,
                          alpha = feature_sig ))+
  
          xlab("estimate") +
          ylab("-log10 p-value") +
          ylim(0,5)+
          xlim(-100,100)+
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


pdf("PRS lipodystrophy profile subq D14 volcano plot.pdf")
plot
dev.off()
##############
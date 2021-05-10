###Lipocyte Profiler on AMSCS treated with isoproteronol and rosi data in AMSCs#####
################################################
#####set working directory#####
setwd("/Users/sophiestrobel/Desktop/LP/LP data import")
####load raw data files#####
stim_plate<-read.csv(file = 'features_well_median_plate_norm_treatment_only.csv')
stim_plate.1<-unite(stim_plate, newcol, c("Metadata_Plate", "Metadata_Well"), remove=FALSE)
########plates with treatment####
treatment_plates<-c("BR00112230", "BR00112232", "BR00112233", "BR00112236", "BR00112238", "BR00112240",  "BR00112242")

#####meta data###
meta_data<-read.csv(file = "MOBB_Phenotype.csv")

layouts<-read.csv(file = 'layouts.csv')
layout_select<-subset(layouts, layouts$plate %in% treatment_plates)
layout_select.1<-unite(layout_select, newcol, c("plate", "well"), remove=FALSE)

treat_layout<-read.csv(file = 'treatment_layout.csv')
treat_layout.1<-unite(treat_layout[, 1:4], newcol, c("Metadata_Plate", "Metadata_Well"), remove=FALSE)

####merge data files#####
merged_stim <- merge(x = layout_select.1, y = treat_layout.1, by = "newcol", sort = FALSE)
merged_stim.1 <- merge(x = merged_stim, y = stim_plate.1, by = "newcol", sort = FALSE)

merged_stimulus<-merged_stim.1[, -c(7,8,9, 11,12)]
names(merged_stimulus)[names(merged_stimulus) == "FFA.x"] <- "FFA"
merged_stimulus.1<-unite(merged_stimulus, Group.1, c("patientID", "treatment", "FFA", "cellType"), remove=FALSE)
merged_stimulus.1<-merged_stimulus.1[, - 1]
merged_stimulus.1[,8:3009] <- sapply(sapply(merged_stimulus.1[,8:3009], as.character),as.numeric) 

#######change ER into BODIPY####
colnames(merged_stimulus.1) <- gsub("_ER", "_BODIPY", colnames(merged_stimulus.1))
###############################
#####include control samples of LP_data.csv###
###load data and subset###
#combined_data<-read.csv(file = 'LP_data.csv')
#combined_data<-combined_data[,-1]
D14<-subset(combined_data, combined_data$day == "14")
treatment<-rep("DMSO", 129)

D14.1<-cbind(treatment, D14)
D14.2<-D14.1[, c(37,32,27,3,33,1,38:3039)]

D14.3<-unite(D14.2, Group.1, c("patientID", "treatment", "FFA", "cellType"), remove=FALSE)
#subset to matching samples with treated samples###
group<-merged_stimulus.1$Group.1

D14.4<-subset(D14.3, D14.3$Group.1 %in% group)
names(D14.4)[names(D14.4) == "Metadata_Plate"] <- "plate"

#####merge data files###
df.1<-rbind(D14.4, merged_stimulus.1)


#####subset data files#####
iso_ffa_vc<-subset(df.1, df.1$cellType == "vc" & df.1$FFA == "1" & df.1$treatment !="rosi" )
#iso_nonffa_vc<-subset(merged_stimulus.1, merged_stimulus.1$cellType == "vc" & merged_stimulus.1$FFA == "0" & merged_stimulus.1$treatment !="rosi" )
#iso_ffa_sc<-subset(merged_stimulus.1, merged_stimulus.1$cellType == "sc" & merged_stimulus.1$FFA == "1" & merged_stimulus.1$treatment !="rosi" )
#iso_nonffa_sc<-subset(merged_stimulus.1, merged_stimulus.1$cellType == "sc" & merged_stimulus.1$FFA == "0" & merged_stimulus.1$treatment !="rosi" )

#####aggregate data you want to normalise####
data<-iso_ffa_vc

agg<-aggregate(data[,8:3009], by = list(data$Group.1), FUN= "median")

agg_meta<-data[!duplicated(data$Group.1), ]
agg_meta<-agg_meta[, 1:7]

######merge aggregated data and meta data####
df<-merge(x = agg_meta, y =agg, by = "Group.1", sort = FALSE)

#####subset in meta and feature data######
meta<-df[, 1:7]
input<-df[, 8:3009]

#####filter data, remove blocklisted features, SmallBODIPY objects and features with 0´s or NA´s across all samples#####
#####see methods LP manuscripte Laber and Strobel et al.#####
input[,] <- sapply(sapply(input[,], as.character),as.numeric) 

back<-as.data.frame(t(input))
var.selection.1<-back[!grepl( "_Costes_" , rownames(back) ) ,  ]
var.selection.2<-var.selection.1[!grepl( "_Manders_" , rownames(var.selection.1)) ,  ]
var.selection.3<-var.selection.2[!grepl( "_RWC_" , rownames(var.selection.2)) ,  ]
var.selection.4<-var.selection.3[!grepl( "SmallBODIPY" , rownames(var.selection.3) ) ,  ]
#
mean<-rowMeans(var.selection.4)
var<-cbind(mean, var.selection.4)

var.1<-subset(var, var$mean != "0")
var.1a<- var.1[complete.cases(var.1), ]

variables<-as.data.frame(t(var.1a[,2:11]))

######merge meta und data####
df.2<-cbind(meta, variables)
######LP ANOVA###########
case<-subset(df.2, df.2$treatment== "iso")
patient<-case$patientID

control<-subset(df.2, df.2$treatment == "DMSO")
####subset control that individual match ###
control<-subset(control, control$patientID %in% patient)

######build dataframe####
data_anov<-rbind(case, control)

#########as numeric features data, change range depending on input####
data_anov[,8:2725] <- sapply(sapply(data_anov[,8:2725], as.character),as.numeric) 
############define variables you want to test
colvector<-as.character(colnames(variables))

library(sjstats)

lm_model<-list()

for(i in colvector){
  
  lm <- lm(data_anov[,i] ~ treatment 
                        + patientID
                        + plate, 
                          data = data_anov)
  
  lm_model[[i]]<-anova(lm)
}

######extract p-values####### 
pvalue<-list()

for(i in colvector){
  
  pvalue[[i]]<-lm_model[[i]][["Pr(>F)"]]
  
}

##########structure pvalue dataframe
pvalue_df<-do.call(cbind, pvalue)
pvalue_df_1<-as.data.frame(t(pvalue_df))
pvalue_df_1<-pvalue_df_1[, -4]
colnames(pvalue_df_1)<-c("stim.pvalue",
                         "patientID.pvalue" ,
                         "plate.pvalue"
                         )

#########calculate q-values (FDR)#####
library(qvalue)

p<-pvalue_df_1$stim.pvalue

qvalues<-qvalue(p)
hist(qvalues)

q<-qvalues[["qvalues"]]

#######combine q-values and p-values#####
pvalue_data<-cbind(pvalue_df_1, q)

##############extract F value#####
Fvalue<-list()
for(i in colvector){
  
  Fvalue[[i]]<-lm_model[[i]][["F value"]]
}

Fvalue_df<-do.call(cbind, Fvalue)
Fvalue_df_1<-as.data.frame(t(Fvalue_df))
Fvalue_df_1<-Fvalue_df_1[, -4]
colnames(Fvalue_df_1)<-c("stim.Fvalue", 
                         "patientID.F" ,
                         "plate.F"
                         )

#######extract effect size (eta_sq)#####
stats<-list()

for(i in colvector){
  
  stats[[i]]<-anova_stats(lm_model[[i]])
  
}

anova_stats(lm_model[[3]])

effect<-sapply(stats,function(x) x[7])
effect_df<-do.call(cbind, effect)
effect_df_1<-as.data.frame(t(effect_df))
effect_df_1<-effect_df_1[, - 4]
colnames(effect_df_1)<-c("stim.eta_sq",  
                         "patientID.eta_sq" ,
                         "plate.eta_sq"
                         )
#
features<-colvector

#####t-test to calculate t-statistic####
#####as numeric and as.matrix case and control dataframe####
mcluster<-control[,8:2725]

mcluster[,] <- sapply(sapply(mcluster[,], as.character),as.numeric) 

matrixcontrol<-as.matrix(mcluster)
#######numeric matrix case 
mcluster<-case[,8:2725]

mcluster[,] <- sapply(sapply(mcluster[,], as.character),as.numeric) 

matrixcase<-as.matrix(mcluster)

##t-test##
outtest<- list ()

for(i in colvector){
  
  x<-matrixcase[,i]
  y<-matrixcontrol[,i]
  
  outtest[[i]]<-t.test(x,y)
}

pvalue_ttest <- data.frame(matrix(unlist(outtest), nrow=length(outtest), byrow=T))

rownames(pvalue_ttest)<-colvector
colnames(pvalue_ttest)<-c("t-test", "df", "pvalue", "confin1", "confin2", "meanx", "meany", "diffmean", "stderr", "alternative", "method", "dataname")


###calculate qvalue, FDR#####
library(qvalue)
p<- as.numeric(as.character(pvalue_ttest$pvalue))

qvalues<-qvalue(p)
hist(qvalues)

q<-qvalues[["qvalues"]]

###combine q-value, p-value and t-statistics#####
pFDR<-cbind(p, q, pvalue_ttest)

ttest_pFDR<-as.data.frame(t(pFDR))
##############structure
ttest_pFDR<-rbind(ttest_pFDR["t-test",],ttest_pFDR["p",],ttest_pFDR["q",] )

ttest_pFDR<-as.data.frame(t(ttest_pFDR))
colnames(ttest_pFDR)<-c("t-test", "p","qvalue")

ttest_pFDR[,] <- sapply(sapply(ttest_pFDR[,], as.character),as.numeric) 

####combine anova an ttest in one output file####

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
                  ifelse(volcano$q <= 0.05 & volcano$stim.pvalue <= 0.05,'5%FDR','not'))


######plot####
library(ggplot2)
plot <- ggplot(volcano, aes(x= `t-test` , y=-log10(stim.pvalue))) +
  geom_point(aes(color=feature_color,
                 size = feature_sig,
                 alpha = feature_sig ))+
  
  xlab("t-statistics") +
  ylab("-log10 p-value") +
  ylim(0,8)+
  xlim(-4,4)+
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


pdf("profile isoproterenol visceral +ffa D14.pdf")
plot
dev.off()

write.table(pvalue_data, file = "profile BCL2-KD data.csv", sep = ",", col.names = NA, qmethod = "double")
#########################################


###Lipocyte Profiler on BCL2-KD data in AMSCs#####
#scripte for loading and merging of data inputs#
################################################
#####set working directory#####
setwd("/Users/sophiestrobel/Desktop/LP/data raw/LP data import")
####load raw data files#####
BCL2<-read.csv(file = "BCL2KD.csv")
BCL2<-BCL2[,-1]

#####structure dataframe#####
names(BCL2)[names(BCL2) == "Metadata_metadata_day"] <- "Day"
names(BCL2)[names(BCL2) == "Metadata_metadata_sample"] <- "sample"

######subset data of siBCL2 and siNT###
BCL2KD<-BCL2[grepl( "si" , BCL2$sample ) ,  ]

#######structure data#####
BCL2.1<-BCL2KD %>% separate(sample, c("sitype", "ID"), " ", extra = "merge")

#######silenced samples during differentiation#####
BCL2.diff<-BCL2.1[grepl( "5 " , BCL2.1$ID ) ,  ]

BCL2.1diff<-BCL2.diff %>% separate(ID, c("KD", "ID"), " ", extra = "merge")
BCL2.diff<-BCL2.1diff[,c(  -3, -4, -5)]
#######silenced samples before differentiation####
BCL2.pre<-BCL2.1[!grepl( "5 " , BCL2.1$ID ) ,  ]
BCL2.pre<-BCL2.pre[,c(  -3, -4, -5)]
KD<-rep("pre", 160)
BCL2.pre<-cbind(BCL2.pre[,1:3], KD, BCL2.pre[,4:3041])

########combine data####
BCL2_combined<-rbind(BCL2.pre,BCL2.diff)

########create new variable for position in the plate#####
library(stringr)
well<-str_split_fixed(BCL2_combined$Metadata_Well, "[[:digit:]]", 2)
colnames(well)<-c("row", "column")

BCL2_combined_data<-cbind(well, BCL2_combined)

# ##########subset data#####
# #########AMSCs ID´s######498 522 523 578 582#####
########into diff = silenced samples during differentiation and pre = silenced samples before differentiation####
BCL2_pre<-subset(BCL2_combined,BCL2_combined$KD == "pre")
BCL2_diff<-subset(BCL2_combined,BCL2_combined$KD == "5")

########subset into days of differentiation#####
day0<-subset(BCL2_pre, BCL2_pre$Day == "0")
day3<-subset(BCL2_pre, BCL2_pre$Day == "3")
day8<-subset(BCL2_pre, BCL2_pre$Day == "8")
day14<-subset(BCL2_pre, BCL2_pre$Day == "14")

#######decide which subset you want to analyse#####
AP_subset<-day14

##############
aggregated<-aggregate(AP_subset[,7:3042], by = list(AP_subset$ID,AP_subset$sitype ), FUN= "median")

meta<-aggregated[, 1:2]
input<-aggregated[, 3:3038]
########ER to BODIPY#####
colnames(input) <- gsub("_ER", "_BODIPY", colnames(input))

#####filter data, remove blocklisted features, SmallBODIPY objects and features with 0´s or NA´s across all samples#####
#####see methods LP manuscripte Laber and Strobel et al.#####
variables<-input

variables[,] <- sapply(sapply(variables[,], as.numeric), as.numeric) 
back<-as.data.frame(t(variables))
# 
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
df<-cbind(meta, variables)

####define case = BCL2-KD and control
case<-subset(df, df$Group.2 == "siBCL2")

control<-subset(df, df$Group.2  == "siNT")
############build matrix for anova######
data_anov<-rbind(case, control)

#########as numeric features data, change range depending on input######
data_anov[,3:2738] <- sapply(sapply(data_anov[,3:2738], as.character),as.numeric) 

#########BCL2 inhibition is pro-apototic######
#######adj with cel density in the well #####
y<-as.data.frame(t(variables))
input.1<-y[!grepl( "Cells_Neighbors_PercentTouching_Adjacent" , rownames(y)) ,  ]

########define LP features you want to test#####
colvector<-rownames(input.1)

########ANOVA#####
library(sjstats)

lm_model<-list()

for(i in colvector){
  
  lm <- lm(data_anov[,i] ~ Group.2 
                        + Cells_Neighbors_PercentTouching_Adjacent, 
                          data = data_anov)
        lm_model[[i]]<-anova(lm)
 }
##########extract p-values######## 
pvalue<-list()

for(i in colvector){
  
  pvalue[[i]]<-lm_model[[i]][["Pr(>F)"]]
  
}

##########structure pvalue dataframe######
pvalue_df<-do.call(cbind, pvalue)
pvalue_df_1<-as.data.frame(t(pvalue_df))
pvalue_df_1<-pvalue_df_1[, -3]
colnames(pvalue_df_1)<-c("KD.pvalue",
                         "density.pvalue") 

##################calculate q-values (FDR)######
library(qvalue)

p<-pvalue_df_1$KD.pvalue

qvalues<-qvalue(p)
hist(qvalues)

q<-qvalues[["qvalues"]]

############combine q-values and p-values########
pvalue_data<-cbind(pvalue_df_1, q)

##############extract F value######
Fvalue<-list()
for(i in colvector){
  
  Fvalue[[i]]<-lm_model[[i]][["F value"]]
}

Fvalue_df<-do.call(cbind, Fvalue)
Fvalue_df_1<-as.data.frame(t(Fvalue_df))
Fvalue_df_1<-Fvalue_df_1[, -3]
colnames(Fvalue_df_1)<-c("KD.F", 
                         "density.F") 

#############extract effect size eta_sq
stats<-list()

for(i in colvector){
  
  stats[[i]]<-anova_stats(lm_model[[i]])
  
}

effect<-sapply(stats,function(x) x[7])
effect_df<-do.call(cbind, effect)
effect_df_1<-as.data.frame(t(effect_df))
effect_df_1<-effect_df_1[, -3]
colnames(effect_df_1)<-c("KD.eta",
                         "density.eta")

features<-colvector

######t-test to calculate t-statistic#######
#####as numeric and as.matrix case and control dataframe#####
mcluster<-control[, 3:2738]

mcluster[,] <- sapply(sapply(mcluster[,], as.character),as.numeric) 

matrixcontrol<-as.matrix(mcluster)
#######numeric matrix case##### 
mcluster<-case[,3:2738]

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
####calculate qvalue, FDR#####
library(qvalue)
p<- as.numeric(as.character(pvalue_ttest$pvalue))

qvalues<-qvalue(p)
hist(qvalues)

q<-qvalues[["qvalues"]]
######### combine q-value, p-value and t-statistics########
pFDR<-cbind(p, q, pvalue_ttest)

ttest_pFDR<-as.data.frame(t(pFDR))
#######structure#########
ttest_pFDR<-rbind(ttest_pFDR["t-test",],ttest_pFDR["p",],ttest_pFDR["q",] )

ttest_pFDR<-as.data.frame(t(ttest_pFDR))
colnames(ttest_pFDR)<-c("t-test", "p","qvalue")

ttest_pFDR[,] <- sapply(sapply(ttest_pFDR[,], as.character),as.numeric) 

##########combine anova an ttest in one output file#######
volcano<-cbind(pvalue_data,effect_df_1,Fvalue_df_1, ttest_pFDR, features)


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
                  ifelse(volcano$q <= 0.01 & volcano$KD.pvalue <= 0.05,'1%FDR','not'))


######plot####
library(ggplot2)
plot <- ggplot(volcano, aes(x= `t-test` , y=-log10(KD.pvalue))) +
            geom_point(aes(color=feature_color,
                 size = feature_sig,
                 alpha = feature_sig ))+
  
            xlab("t-statistics") +
            ylab("-log10 p-value") +
            ylim(0,10)+
            xlim(-25,25)+
            scale_size_manual(name = "", 
                    values = c("1%FDR" = 1.5, "not" = 1))+
            scale_alpha_manual(name = "", 
                     values = c("1%FDR" = 0.75, "not" = 0.1))+
            scale_color_manual(name = "", 
                     values = c("Mito" = "#f56464" ,
                                "AGP" = "#f2cf41",
                                "Lipid" = "#4dac26",
                                "DNA" = "#65a6db",
                                "other" = "#959ca3"))+
            theme_bw() +
            theme(strip.background = element_rect(colour = "black",
                                        fill = "#fdfff4"))


pdf("profile BCL2-KD D14.pdf")
plot
dev.off()

write.table(pvalue_data, file = "profile BCL2-KD data.csv", sep = ",", col.names = NA, qmethod = "double")

###############################

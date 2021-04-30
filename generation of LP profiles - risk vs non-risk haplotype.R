#####Generation of LP profile of risk and non-risk haplotype####
#####set working directory###
setwd("/Users/sophiestrobel/Desktop/LP/data raw/LP data import")

####load genotype data ####
#SNPs<-read.csv(file = 'singlelocus.csv')

#merged_locus<-merge(x = SNPs, y = combined_data, by = "patientID", all.y = T, sort = FALSE)

###load LP data set#####
combined_data<-read.csv(file = "LP_data.csv.csv")


######subset depending on which anylsis you want to run
#####susbet depot
all_data_nonFFA<-subset(combined_data, combined_data$FFA == "0")
all_data_FFA<-subset(combined_data, combined_data$FFA == "1")

# data_nonFFA_day0_vc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "vc" & all_data_nonFFA$day == "0")
# data_nonFFA_day3_vc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "vc" & all_data_nonFFA$day == "3")
# data_nonFFA_day8_vc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "vc" & all_data_nonFFA$day == "8")
# data_nonFFA_day14_vc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "vc" & all_data_nonFFA$day == "14")
# 
# data_FFA_day0_vc<-subset(all_data_FFA, all_data_FFA$cellType == "vc" & all_data_FFA$day == "0")
# data_FFA_day3_vc<-subset(all_data_FFA, all_data_FFA$cellType == "vc" & all_data_FFA$day == "3")
# data_FFA_day8_vc<-subset(all_data_FFA, all_data_FFA$cellType == "vc" & all_data_FFA$day == "8")
# data_FFA_day14_vc<-subset(all_data_FFA, all_data_FFA$cellType == "vc" & all_data_FFA$day == "14")
# 
# data_nonFFA_day0_sc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "sc" & all_data_nonFFA$day == "0")
# data_nonFFA_day3_sc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "sc" & all_data_nonFFA$day == "3")
# data_nonFFA_day8_sc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "sc" & all_data_nonFFA$day == "8")
# data_nonFFA_day14_sc<-subset(all_data_nonFFA, all_data_nonFFA$cellType == "sc" & all_data_nonFFA$day == "14")

#############input data you want to analyse##########
AP_subset<-data_nonFFA_day8_sc
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


####multi-way analysis of variance (ANOVA)#####
####case = risk allele
####control = non-risk allele

case<-subset(df, df$rs12454712 == "TT")

control<-subset(df, df$rs12454712  == "CC")

#####merge case and control###
data_anov<-rbind(case, control)

####transform LP features to numeric and change range depending on input####
data_anov[,37:2749] <- sapply(sapply(data_anov[,37:2749], as.character),as.numeric) 

###define variables you want to test###
colvector<-as.character(colnames(variables))

#####ANOVA#####
library(sjstats)

lm_model<-list()

for(i in colvector){
  
  lm <- lm(data_anov[,i] ~ rs12454712 +BMI +age +batch + sex , data = data_anov)
  
  lm_model[[i]]<-anova(lm)
  
  }
####extract p-values of ANOVA###
pvalue<-list()

for(i in colvector){
  
  pvalue[[i]]<-lm_model[[i]][["Pr(>F)"]]
  
}
###structure pvalue dataframe####
pvalue_df<-do.call(cbind, pvalue)
pvalue_df_1<-as.data.frame(t(pvalue_df))
pvalue_df_1<-pvalue_df_1[, -6]
colnames(pvalue_df_1)<-c("SNP.pvalue","BMI.pvalue", "age.pvalue", "sex.pvalue", "batch.pvalue")

####calculate q-values (FDR)#####
library(qvalue)

p<-pvalue_df_1$SNP.pvalue

qvalues<-qvalue(p)
hist(qvalues)

q<-qvalues[["qvalues"]]

############combine q-values and p-values
pvalue_data<-cbind(pvalue_df_1, q)

#####extract F value#####
Fvalue<-list()
for(i in colvector){
  
  Fvalue[[i]]<-lm_model[[i]][["F value"]]
}

Fvalue_df<-do.call(cbind, Fvalue)
Fvalue_df_1<-as.data.frame(t(Fvalue_df))
Fvalue_df_1<-Fvalue_df_1[, -6]
colnames(Fvalue_df_1)<-c("SNP.Fvalue","BMI.Fvalue", "age.F", "sex.F", "batch.F")

####extract effect size eta_sq#####
stats<-list()

for(i in colvector){
  
  stats[[i]]<-anova_stats(lm_model[[i]])
  
}

anova_stats(lm_model[[3]])

effect<-sapply(stats,function(x) x[7])
effect_df<-do.call(cbind, effect)
effect_df_1<-as.data.frame(t(effect_df))
effect_df_1<-effect_df_1[, -6]
colnames(effect_df_1)<-c("SNP.eta_sq","BMI.eta_sq", "age.eta_sq", "sex.eta_sq", "batch.eta_sq")

features<-colvector

####t-test to calculate t-statistic######
####as numeric and as.matrix case and control dataframe#####
mcluster<-control[,37:2749]

mcluster[,] <- sapply(sapply(mcluster[,], as.character),as.numeric) 

matrixcontrol<-as.matrix(mcluster)
#######numeric matrix case 
mcluster<-case[,37:2749]

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

####calculate qvalue (FDR)####
library(qvalue)
p<- as.numeric(as.character(pvalue_ttest$pvalue))

qvalues<-qvalue(p)
hist(qvalues)

q<-qvalues[["qvalues"]]

######### combine q-value, p-value and t-statistics######
pFDR<-cbind(p, q, pvalue_ttest)

ttest_pFDR<-as.data.frame(t(pFDR))
##############structure
ttest_pFDR<-rbind(ttest_pFDR["t-test",],ttest_pFDR["p",],ttest_pFDR["q",] )

ttest_pFDR<-as.data.frame(t(ttest_pFDR))
colnames(ttest_pFDR)<-c("t-test", "p","qvalue")

ttest_pFDR[,] <- sapply(sapply(ttest_pFDR[,], as.character),as.numeric) 

#####combine anova and ttest in one output file######
volcano<-cbind(pvalue_data,effect_df_1,Fvalue_df_1, ttest_pFDR, features)

###########plot  volcano#####
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
                  ifelse(volcano$q <= 0.05 & volcano$SNP.pvalue <= 0.05,'5%FDR','not'))


######plot####
library(ggplot2)
plot <- ggplot(volcano, aes(x= `t-test` , y=-log10(SNP.pvalue))) +
              geom_point(aes(color=feature_color,
                 size = feature_sig,
                 alpha = feature_sig ))+
  
              xlab("estimate") +
              ylab("-log10 p-value") +
              ylim(0,4)+
              xlim(-5,5)+
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


pdf("risk vs non-risk profile subq D8 volcano plot.pdf")
plot
dev.off()
######
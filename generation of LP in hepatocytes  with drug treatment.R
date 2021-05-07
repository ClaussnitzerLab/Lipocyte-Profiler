#########Lipocyte Profiler om hepatocytes######
##set working directory###
setwd("/Users/sophiestrobel/Desktop/LP/Hepatocytes 2021")
####load input data###
#####4% oxygen#####
plate1<-"BR00101180_normalized.gct"

data_low_oxygen<-read.table(plate1, skip = 2, header = TRUE, sep = "\t")
df_low_oxygen<-as.data.frame(t(data_low_oxygen))
colnames(df_low_oxygen)<-as.character(data_low_oxygen$id)
df_low_oxygen<-df_low_oxygen[-c(1,2), ]                   
#####21% oxygen#####
plate2<-"BR00101181_normalized.gct"

data_normal_oxygen<-read.table(plate2, skip = 2, header = TRUE, sep = "\t")
df_normal_oxygen<-as.data.frame(t(data_normal_oxygen))
colnames(df_normal_oxygen)<-as.character(data_normal_oxygen$id)
df_normal_oxygen<-df_normal_oxygen[-c(1,2), ]  
########merge data plates####
data_hep<-rbind(df_low_oxygen, df_normal_oxygen)

########define new variables column and well position in the well####
library(stringr)
well<-str_split_fixed(data_hep$Well, "[[:digit:]]", 2)
colnames(well)<-c("row", "column")

data_hep<-cbind(well, data_hep)
data_hep<-data_hep[, -c(3:8)]
######split data into meta and LP data###
meta<-data_hep[, 1:4]
input<-data_hep[, 5:3040]

########LP data numeric####
input[,] <- sapply(sapply(input[,], as.character),as.numeric) 
########ER to BODIPY#####
colnames(input) <- gsub("_ER", "_BODIPY", colnames(input))
#####filter data, remove blocklisted features, SmallBODIPY objects and features with 0´s or NA´s across all samples#####
#####see methods LP manuscripte Laber and Strobel et al.#####
back<-as.data.frame(t(input))
var.selection.1<-back[!grepl( "_Costes_" , rownames(back) ) ,  ]
var.selection.2<-var.selection.1[!grepl( "_Manders_" , rownames(var.selection.1)) ,  ]
var.selection.3<-var.selection.2[!grepl( "_RWC_" , rownames(var.selection.2)) ,  ]
var.selection.4<-var.selection.3[!grepl( "SmallBODIPY" , rownames(var.selection.3) ) ,  ]

#########remove LP features with 0´s or NA´s across data###
mean<-rowMeans(var.selection.4)
var<-cbind(mean, var.selection.4)

var.1<-subset(var, var$mean != "0")
var.1a<- var.1[complete.cases(var.1), ]
########adj n number####
variables<-as.data.frame(t(var.1a[,2:145]))

#######merge data####
df<-cbind(meta, variables)

######choose condition to analyse########
normal<-subset(df, df$condition_O2== "21")

### choose drug to analyse####
###split into treatment group and control group####
#1uM Rotenone
#0.5mM DMOG
#5mM metformin
#0.3mM OA

case<-subset(normal, normal$treatment == "0.3mM OA" )

control<-subset(normal, normal$treatment == "CP CTRL")

############build matrix for anova####
data_anov<-rbind(case, control)

#########as numeric features data, change range depending on input#####
data_anov[,5:2741] <- sapply(sapply(data_anov[,5:2741], as.character),as.numeric) 
#####define variables you want to test######
colvector<-as.character(colnames(variables))

######ANOVA#####
library(sjstats)

lm_model<-list()

for(i in colvector){
  
  lm <- lm(data_anov[,i] ~ treatment
                          + row, 
                           data = data_anov)
  
  lm_model[[i]]<-anova(lm)
}

####extract p-values##### 
pvalue<-list()

for(i in colvector){
  
  pvalue[[i]]<-lm_model[[i]][["Pr(>F)"]]
  
}
####structure pvalue dataframe
pvalue_df<-do.call(cbind, pvalue)
pvalue_df_1<-as.data.frame(t(pvalue_df))
pvalue_df_1<-pvalue_df_1[, -3]
colnames(pvalue_df_1)<-c("treat.pvalue","pos.pvalue")
##################calculate q-values (FDR)
library(qvalue)

p<-pvalue_df_1$treat.pvalue

qvalues<-qvalue(p)
hist(qvalues)

q<-qvalues[["qvalues"]]

##########combine q-values and p-values
pvalue_data<-cbind(pvalue_df_1, q)

##############extract F value
Fvalue<-list()
for(i in colvector){
  
  Fvalue[[i]]<-lm_model[[i]][["F value"]]
}

Fvalue_df<-do.call(cbind, Fvalue)
Fvalue_df_1<-as.data.frame(t(Fvalue_df))
Fvalue_df_1<-Fvalue_df_1[, -3]
colnames(Fvalue_df_1)<-c("treat.Fvalue", "pos.F")
#############extract effect size eta_sq
stats<-list()

for(i in colvector){
  
  stats[[i]]<-anova_stats(lm_model[[i]])
  
}

effect<-sapply(stats,function(x) x[7])
effect_df<-do.call(cbind, effect)
effect_df_1<-as.data.frame(t(effect_df))
effect_df_1<-effect_df_1[, -3]
colnames(effect_df_1)<-c("treat.eta_sq",  "pos.eta_sq")

features<-colvector

######t-test to calculate t-statistic#####
#####as numeric and as.matrix case and control dataframe
mcluster<-control[,5:2741]

mcluster[,] <- sapply(sapply(mcluster[,], as.character),as.numeric) 

matrixcontrol<-as.matrix(mcluster)
#######numeric matrix case 
mcluster<-case[,5:2741]

mcluster[,] <- sapply(sapply(mcluster[,], as.character),as.numeric) 

matrixcase<-as.matrix(mcluster)

##t-test###
outtest<- list ()

for(i in colvector){
  
  x<-matrixcase[,i]
  y<-matrixcontrol[,i]
  
  outtest[[i]]<-t.test(x,y)
}

pvalue_ttest <- data.frame(matrix(unlist(outtest), nrow=length(outtest), byrow=T))

rownames(pvalue_ttest)<-colvector
colnames(pvalue_ttest)<-c("t-test", "df", "pvalue", "confin1", "confin2", "meanx", "meany", "diffmean", "stderr", "alternative", "method", "dataname")


####calculate qvalue, FDR##########
library(qvalue)
p<- as.numeric(as.character(pvalue_ttest$pvalue))

qvalues<-qvalue(p)
hist(qvalues)

q<-qvalues[["qvalues"]]

######### combine q-value, p-value and t-statistics
pFDR<-cbind(p, q, pvalue_ttest)

ttest_pFDR<-as.data.frame(t(pFDR))
##############structure
ttest_pFDR<-rbind(ttest_pFDR["t-test",],ttest_pFDR["p",],ttest_pFDR["q",] )

ttest_pFDR<-as.data.frame(t(ttest_pFDR))
colnames(ttest_pFDR)<-c("t-test", "p","qvalue")

ttest_pFDR[,] <- sapply(sapply(ttest_pFDR[,], as.character),as.numeric) 

##########combine anova an ttest in one output file
volcano<-cbind(pvalue_data, effect_df_1, Fvalue_df_1 , ttest_pFDR, features )

###########plot  volcano####
#######assign  categories to channel####
volcano <- volcano %>%
  dplyr::mutate(feature_color = 
                  ifelse(grepl("BODIPY", rownames(volcano)),'Lipid', 
                         ifelse(grepl("AGP", rownames(volcano)),'AGP', 
                                ifelse(grepl("Mito", rownames(volcano)),'Mito',
                                       ifelse(grepl("DNA", rownames(volcano)),'DNA',
                                              'other')))))

volcano$feature_color<-factor(volcano$feature_color,
                              levels = c("Mito","AGP","Lipid","DNA", "other"))

#######delete NA´s####
volcano<- volcano[complete.cases(volcano), ]

#######set signifance level#####
volcano <- volcano %>%
  dplyr::mutate(feature_sig = 
                  ifelse(volcano$treat.pvalue  <= 0.05 &  volcano$q <= 0.001,'0.1%FDR',
                         ifelse(volcano$p  <= 0.05 &  volcano$q <= 0.01,'1%FDR',
                         
                         'not')))

######plot#####
library (ggplot2)

plot <- ggplot(volcano, aes(x= `t-test`, y=-log10(treat.pvalue ))) +
              geom_point(aes(color=feature_color,
                 size = feature_sig,
                 alpha = feature_sig ))+
            xlab("t-statistics") +
            ylab("-log10 p-value") +
            #ylim(0,9)+
            #xlim(-100,100)+
            scale_size_manual(name = "", 
                    values = c("0.1%FDR" = 1.5,"1%FDR" = 1, "not" = 0.5))+
            scale_alpha_manual(name = "", 
                     values = c("0.1%FDR" = 0.75,"1%FDR" = 0.4, "not" = 0.05))+
            scale_color_manual(name = "", 
                     values = c("Mito_feature" = "#f56464" ,
                                "AGP_feature" = "#f2cf41",
                                "bodipy_feature" = "#4dac26",
                                "DNA_feature" = "#65a6db",
                                "other_feature" = "#959ca3"))+
            theme_bw() +
            theme(strip.background = element_rect(colour = "black",
                                        fill = "#fdfff4"))

################
pdf("LP profile OA 21% oxygen hepatocytes.pdf")
plot
dev.off()
############################
write.table(volcano, file = "OA 21%oxygen hepatocytes.csv", sep = ",", col.names = NA,
            qmethod = "double")



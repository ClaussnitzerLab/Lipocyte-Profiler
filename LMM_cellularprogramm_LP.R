#####LMM of transcriptome and LipocyteProfile####
#####cellular programms####
setwd("/Users/slaber/Desktop/R")
#####load input files#####
df<-read.csv(file="AP_RNA_alldata_0510.csv")
df<-df[,-1]
####subset####
###depot###day###FFA-treatment#####
# sc_nonffa_day14<-subset(df, df$cellType == "sc" & df$FFA == "0" & df$Day == "14")
# sc_nonffa_day8<-subset(df, df$cellType == "sc" & df$FFA == "0" & df$Day == "8")
# sc_nonffa_day3<-subset(df, df$cellType == "sc" & df$FFA == "0" & df$Day == "3")
# sc_nonffa_day0<-subset(df, df$cellType == "sc" & df$FFA == "0" & df$Day == "0")
# 
# sc_ffa_day14<-subset(df, df$cellType == "sc" & df$FFA == "1" & df$Day == "14")
# sc_ffa_day8<-subset(df, df$cellType == "sc" & df$FFA == "1" & df$Day == "8")
# sc_ffa_day3<-subset(df, df$cellType == "sc" & df$FFA == "1" & df$Day == "3")
# 
# vc_nonffa_day14<-subset(df, df$cellType == "vc" & df$FFA == "0" & df$Day == "14")
# vc_nonffa_day8<-subset(df, df$cellType == "vc" & df$FFA == "0" & df$Day == "8")
# vc_nonffa_day3<-subset(df, df$cellType == "vc" & df$FFA == "0" & df$Day == "3")
# vc_nonffa_day0<-subset(df, df$cellType == "vc" & df$FFA == "0" & df$Day == "0")
# 
# vc_ffa_day14<-subset(df, df$cellType == "vc" & df$FFA == "1" & df$Day == "14")
# vc_ffa_day8<-subset(df, df$cellType == "vc" & df$FFA == "1" & df$Day == "8")
# vc_ffa_day3<-subset(df, df$cellType == "vc" & df$FFA == "1" & df$Day == "3")

#########define which subset you want to run#####
model_input<-sc_nonffa_day14

##############define AP-features you want to run#####
AP_features<-model_input[,52182:55183]
meta<-model_input[, 1:11]

#######structure of data ####
variables <-AP_features
#######make sure all variables are numeric####
variables[,] <- sapply(sapply(variables[,], as.numeric), as.numeric)

#######remove features by name you don't want#####
remove<-as.data.frame(t(variables))
var.selection.1<-remove[!grepl( "_Costes_" , rownames(remove) ) ,  ]
var.selection.2<-var.selection.1[!grepl( "_Manders_" , rownames(var.selection.1)) ,  ]
var.selection.3<-var.selection.2[!grepl( "_RWC_" , rownames(var.selection.2)) ,  ]
var.selection.4<-var.selection.3[!grepl( "SmallBODIPY" , rownames(var.selection.3) ) ,  ]
######remove 0 features#######
mean<-rowMeans(var.selection.4)
var<-cbind(mean, var.selection.4)
var.1<-subset(var, var$mean != "0")
######remove NA features#####
var.1a<- var.1[complete.cases(var.1), ]

#####set n numbers###
variables<-as.data.frame(t(var.1a[,2:27]))

####define with cellular programm you want to use###
####load HALLMARK genes into R#####
####download genset first from http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp #####
####transfer geneset from downloads to working directory###
library(sigPathway)

GSEA<-importGeneSets("geneset.apotosis.grp", verbose = TRUE)

hallmark_apoptosis<-GSEA[[1]][["probes"]]
hallmark_apoptosis<-hallmark_apoptosis[-1]

####load csv file with gene annotaions#####
gene_symbol<-read.csv(file="gene_symbol.csv")
gene_symbol<-gene_symbol[,-1]
colnames(gene_symbol)<-c("gene", "EN_number")

########map gene name with ensemble gene identifier####
genset_hallmark_apoptosis<-subset(gene_symbol, gene_symbol$gene %in% hallmark_apoptosis)

########selecte genes for LMM ####
apoptosis<-intersect(genset_hallmark_apoptosis$EN_number, colnames(model_input))

#########combine input file#####
data<-cbind(meta, model_input[, apoptosis], variables)
data[, 171:175]
#####define outcome and predicting variables#####
####outcome variables = AP features######
####adj out_start depending how many genes you use####
####adj out_end to last column of dataframe data####
# outcome
out_start=172
out_end= 3173
out_nvar=out_end-out_start+1

out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_pvalue=rep(NA, out_nvar)

######predicting variables = genes#####
####adj exp_end depending how may genes you analyse####
#exposure
exp_start=11
exp_end=171
exp_nvar=exp_end-exp_start+1

exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_se = rep(NA, out_nvar)
exp_pvalue=rep(NA, exp_nvar)

number=1
########model#######
library(lme4)
for (i in out_start:out_end){
  outcome = colnames(data)[i]
  for (j in exp_start:exp_end){
    exposure = colnames(data)[j]
    
    model <- lmer(get(outcome) ~ get(exposure) + age + BMI + sex + T2D + (1|batch),
                  na.action = na.exclude,
                  data=data)
    
    Vcov <- vcov(model, useScale = FALSE)
    beta <- fixef(model)
    se <- sqrt(diag(Vcov))
    zval <- beta / se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    
    out_beta[number] = as.numeric(beta[2])
    out_se[number] = as.numeric(se[2])
    out_pvalue[number] = as.numeric(pval[2])
    out_variable[number] = outcome
    number = number + 1
    
    exp_beta[number] = as.numeric(beta[2])
    exp_se[number] = as.numeric(se[2])
    exp_pvalue[number] = as.numeric(pval[2])
    exp_variable[number] = exposure
    number = number + 1
  }
}

outcome = data.frame(out_variable, out_beta, out_se, out_pvalue)
exposure = data.frame(exp_variable, exp_beta, exp_se, exp_pvalue)

library(tidyverse)
outcome = outcome %>% 
  rename(
    variable = out_variable,
    beta = out_beta,
    se = out_se,
    pvalue = out_pvalue
    
  )
exposure = exposure %>% 
  rename(
    variable = exp_variable,
    beta = exp_beta,
    se = exp_se,
    pvalue = exp_pvalue
    
  )

#########combine outcome and 
all = rbind(outcome, exposure)

all_outcome = na.omit(outcome)
all_exposure = na.omit(exposure)


model_output<-cbind(all_outcome, all_exposure)
model_output.1<-model_output[,1:5]

names(model_output.1)[names(model_output.1) == "variable.1"]<- "EN_number"

######calculate FDR#####
library(qvalue)
p<-model_output.1$pvalue
qvalues<-qvalue(p)
hist(qvalues)
q_value<-qvalues[["qvalues"]]

model_output_final<-cbind(model_output.1, q_value)
#########build matrix LP features  and genes ####
matrix<-model_output_final[, c("variable", "EN_number","beta" )]
matrix.1<-spread(matrix, key="EN_number", value= "beta", fill = NA, convert = FALSE, drop = TRUE, sep = NULL)
#######structure matrix#####
rownames(matrix.1)<-matrix.1$variable
matrix.1<-matrix.1[, -1]
matrix_final<-t(matrix.1)

#######add gene names to matrix####
EN_gene<-as.data.frame(rownames(matrix_final))
colnames(EN_gene)<-"EN_number"

data_output<-merge(EN_gene,gene_symbol, "EN_number", all.x =T, sort=F)

apoptosis_LMM<-cbind(data_output, matrix_final)

write.csv(apoptosis_LMM, file="apoptosis_LMM_subq -ffa.csv")

########


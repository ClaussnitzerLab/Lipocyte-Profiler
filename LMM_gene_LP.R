#####LMM of transcriptome and LipocyteProfile####
#####set directory####
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

######define AP-features you want to run#####
AP_features<-model_input[,52182:55183]

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

#######adj. n number####
variables<-as.data.frame(t(var.1a[,27]))

data<-cbind(model_input[,1:52181], variables)

#####define outcome and predicting variables#####
####outcome variables = AP features######
# outcome
out_start=52182
out_end= 54898
out_nvar=out_end-out_start+1

out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_pvalue=rep(NA, out_nvar)

######predicting variables = genes#####
#exposure
exp_start=11
exp_end=52181
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

#########combine outcome and and exposure#####
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
###############
write.csv(model_output_final, file="LMM subq D14 -ffa.csv")



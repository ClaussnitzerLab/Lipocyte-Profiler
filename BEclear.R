###Batch effect with BEclear package on Lipocyte Profiler data of on hWAT hBAT SGBS, across 3 batches with and without FFA treatment.
################################################
#####set working directory#####
setwd("/Users/sophiestrobel/Desktop/LP/data raw/LP data import")
####load raw data files#####
input<-read.csv(file = 'LP_hWAT_hBAT_SGBS_data.csv')
input<-input[,-1]
############subset data#####
nonffa<-subset(input, input$FFA == "0")

######split data in input and meta data#####
input<- as.matrix(t(nonffa[, 9:3010]))
########define variables####
batch_id<-as.character(nonffa$batch)
day<-as.character(nonffa$day)
sample_id<-as.character(nonffa$newcol)
celltype_id<-as.character(nonffa$patientID)
##########combine variables#####
samples<-as.data.frame(cbind(sample_id, batch_id, day, celltype_id))

colnames(input)<-sample_id

################
library(BEclear)
batchEffects<-calcBatchEffects(input, samples, adjusted=TRUE, method="fdr")
#####extract feature medians and p-values###
med <- batchEffects$med
pvals <- as.data.frame(batchEffects$pval)

####structure data for plotting####
med.1<-tidyr::gather(med)
colnames(med.1)<-c("key", "median")

pvals.1<-tidyr::gather(pvals)
colnames(pvals.1)<-c("key", "pvals")

data<-cbind(med.1, pvals.1[, 2])
colnames(data)<-c("key", "median","pvals")

#####plot#####
plot <- ggplot(data, aes(x=key ,y=median)) + 
  geom_boxplot()+
  geom_jitter(aes(fill=-log10(pvals)),position=position_dodge(0.8), size=1, shape=21)+
  theme_bw() 

##Summarize p-values and median differences for batch effect features#####
sum <- calcSummary(medians = med, pvalues = pvals)

write.table(sum, file = "batch effect features BEclear.csv", sep = ",", col.names = NA, qmethod = "double")
##Calculates the score table#####
score.table <- calcScore(data = input, samples = samples, summary = sum)
write.table(score.table, file = "score.table BEclear.csv", sep = ",", col.names = NA, qmethod = "double")

## Simple boxplot for the example data separated by samples####
makeBoxplot(data = input, samples = samples, score = score.table, bySamples = TRUE, main = "beta across samples")
##Simple boxplot for the example data separated by batch####
makeBoxplot(input, samples, score.table, bySamples = F, col = "standard", main = "batch detection data", xlab = "Batch", ylab = "feature median", scoreCol = TRUE)

###################

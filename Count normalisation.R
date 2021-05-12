#####normalization of raw counts (batch 3 of RNA-seq) using DESeq package####
#####set working directory#####
setwd("/Users/slaber/Desktop/github")
#####load featuresCounts######
RNA.new<-read.table(file = 'rnaseq_feature_counts_adipose.tsv', sep = '\t', header = TRUE)
##
rownames(RNA.new)<- RNA.new$Geneid
RNA.new<-RNA.new[, 2:82]

########exclude BAM file that didn't pass quality control#####
bad<-c("X280scD14_DMSO.bam","X256vcD3plusFFA.bam","X129vcD3.bam","X216scD14_DMSO.bam","X216scD14plusFFA.bam")
#
RNA_new.1<-as.data.frame(t(RNA.new))
RNA.new.good<-subset(RNA_new.1,! rownames(RNA_new.1) %in% bad)

#######move on with good quality reads#####
sample<-rownames(RNA.new.good)
RNA_new.2<-cbind(sample,RNA.new.good)

#####structure sample information data####
#####split xx 
z<-strsplit(gsub("(([A-Z])([0-9]*))([a-z]*)(([A-Z])([0-9]))", "\\1 \\2 \\3 \\4 \\5 \\6", RNA_new.2$sample), " ")
df <- data.frame(matrix(unlist(z), nrow=76, byrow=T), stringsAsFactors=FALSE)
df$X1<-gsub("X", "m",df$X1 )

df<-df[, -c(2:3)]
colnames(df)<-c("patientID", "cellType", "Day", "treat")

#######load sample meta data####
treat<-read.csv(file="treatement_layout.csv")
treat<-treat[,-1]

#######merge treatment layout and RNA-seq counts####
RNA.new.data<-cbind( treat[1:76,4],df[1:76,1:3],RNA_new.1[1:76,2:57821])
names(RNA.new.data)[names(RNA.new.data) == "treat[1:76, 4]"]<- "FFA"
rownames(RNA.new.data)<- rownames(RNA_new.1[1:76,])

#########normalize counts#######
########split into counts and meta data####
count<-RNA.new.data[, 5:57824]

meta<-RNA.new.data[, 1:4]
meta1<-unite(meta, newcol, c("patientID", "cellType", "FFA", "Day"))

##
count<-as.data.frame(t(count))

count[, ]<-sapply(sapply(count[,], as.character),as.numeric) 

countData<-as.matrix(count)
colnames(countData) <- NULL

condition<-meta1$newcol

library(DESeq2)
######create DEseq dataframe###
dds <- DESeqDataSetFromMatrix(round(countData), 
                              DataFrame(condition), ~ condition)
#
View(counts(dds))
#####create size factor#####
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
#####normalise####
normalized_counts <- counts(dds, normalized=TRUE)
colnames(normalized_counts)<-meta1$newcol
new<-as.data.frame(t(normalized_counts))

RNA.new.data.norm<-cbind(meta,new )

######structure column names####
en<-as.data.frame(rownames(normalized_counts))
colnames(en)<-"EN"
z<-strsplit(en$EN, ".", fixed = TRUE)
df.1 <- data.frame(matrix(unlist(z), nrow=57820, byrow=T), stringsAsFactors=FALSE)

col<-c(colnames(meta),df.1$X1)
colnames(RNA.new.data.norm)<-col

#########
write.table(RNA.new.data.norm, file="RNA-seq normalised batch3.csv", sep=",", col.names = NA, qmethod = "double")
#########
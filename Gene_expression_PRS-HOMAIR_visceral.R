###########
#####differences of gene expression PRS-HOMA####
#####set directory####
setwd("/Users/slaber/Desktop/")
#####load input files#####
data<-read.csv(file="RNA_features_022521_AP.csv")
data<-data[, -1]

######exlcude LP data and just keep RNA-seq data####
data_RNA<-data[, c(1:29, 3032:55206)]
######subset data to files you want to analyis#####
day14_ffa_vc<-subset(data_RNA, data_RNA$FFA == "1" & data_RNA$cellType == "vc" & data_RNA$Day == "14")
# day8_ffa_vc<-subset(data_RNA, data_RNA$FFA == "1" & data_RNA$cellType == "vc" & data_RNA$Day == "8")
# day3_ffa_vc<-subset(data_RNA, data_RNA$FFA == "1" & data_RNA$cellType == "vc" & data_RNA$Day == "3")
# day0_ffa_vc<-subset(data_RNA, data_RNA$FFA == "0" & data_RNA$cellType == "vc" & data_RNA$Day == "0")
###
input<-day14_ffa_vc
#####define set of genes you want to analyse#####
#####download geneset from GESA
#####eg. hallmarks adipoctes, lipid metabolism, glycolysis#####
library(sigPathway)
gene_symbol<-read.csv(file="gene_symbol.csv")
gene_symbol<-gene_symbol[, -1]
######
GSEA<-importGeneSets("geneset.grp", verbose = TRUE)
GSEA_lipid<-importGeneSets("geneset.lipid.grp", verbose = TRUE)
GSEA_gyco<-importGeneSets("genesetglyco.grp", verbose = TRUE)

lipid_panel<-GSEA_lipid[[1]][["probes"]]
panel_lipid<-lipid_panel[-1]

panel_glyco<-GSEA_gyco[[1]][["probes"]]
panel_glyco<-panel_glyco[-1]

panel_adipo<-GSEA[[1]][["probes"]]
panel_adipo<-panel_adipo[-1]

selection_glyco_gene<-subset(gene_symbol, gene_symbol$RNA_big.symbol %in% panel_glyco)
selection_adipo_gene<-subset(gene_symbol, gene_symbol$RNA_big.symbol %in% panel_adipo)
selection_lipid_gene<-subset(gene_symbol, gene_symbol$RNA_big.symbol %in% panel_lipid)

#######comnbine gene sets and remove dublicated genes###
gene_set<-rbind(selection_glyco_gene, selection_adipo_gene, selection_lipid_gene)
gene_set_uni<-gene_set[!duplicated(gene_set$RNA_big.symbol), ]

#######define genes to analyis#####
for(i in colvector){
  lm_model[[i]] <- summary.lm(lm(input[,i] ~ SCORESUM 
                                        + age
                                        + BMI
                                        + sex
                                        + batch, data = input))
  lm_anov[[i]]<-anova(lm(input[,i] ~ SCORESUM 
                                  + age
                                  + BMI
                                  + sex
                                  + batch, data = input))
}
########
############extract p-values
pvalue.l<-list()
for(i in colvector){
  pvalue.l[[i]]<-lm_model[[i]][["coefficients"]]
}
#########
##########structure pvalue dataframe
coef<-as.data.frame(sapply(pvalue.l,function(x) x[2]))
colnames(coef)<-"estimate"
pvalue_df<-sapply(pvalue.l,function(x) x[,4])
pvalue_df_1<-as.data.frame(t(pvalue_df))
colnames(pvalue_df_1)<-c("intersect", "score.pvalue", "age.pvalue", "BMI.pvalue","sex.pvalue", "batch.pvalue")

##################calculate q-values (FDR)
library(qvalue)
p<-pvalue_df_1$score.pvalue
qvalues<-qvalue(p)
hist(qvalues)
q.l<-qvalues[["qvalues"]]
######################
pvalue.a<-list()
for(i in colvector){
  pvalue.a[[i]]<-lm_anov[[i]][["Pr(>F)"]]
}
##########structure pvalue dataframe
pvalue_df.a<-do.call(cbind, pvalue.a)
pvalue_df.a_1<-as.data.frame(t(pvalue_df.a))
pvalue_df.a_1<-pvalue_df.a_1[, -6]
colnames(pvalue_df.a_1)<-c("score.anova.pvalue", "age.anova.pvalue", "BMI.anova.pvalue" , "sex.anova.pvalue" ,  "batch.anova.pvalue")
#“sex.anova.pvalue”,
########
p<-pvalue_df.a_1$score.anova.pvalue
qvalues<-qvalue(p)
hist(qvalues)
q.a<-qvalues[["qvalues"]]

########
stats<-list()
for(i in colvector){
  stats[[i]]<-anova_stats(lm_anov[[i]])
}
anova_stats(lm_anov[[3]])
effect<-sapply(stats,function(x) x[7])
effect_df<-do.call(cbind, effect)
effect_df_1<-as.data.frame(t(effect_df))
effect_df_1<-effect_df_1[, -6]
colnames(effect_df_1)<-c("score.anova.eta_sq", "age.anova.eta_sq", "BMI.anova.eta_sq", "sex.anova.eta_sq", "batch.anova.eta_sq")
##################
EN_number <-colvector

#####merge data files#####
pvalue_data<-cbind(pvalue_df_1, q.l, pvalue_df.a_1, q.a,coef, EN_number )

#####plot#####
volcano<-pvalue_data

volcano<-merge(volcano, gene_symbol, by="EN_number")

volcano<- volcano[complete.cases(volcano), ]

volcano <- volcano %>%
  dplyr::mutate(label =
                  ifelse(volcano$score.anova.pvalue <= 0.05 & volcano$q.a <= 0.05,"5%FDR",
                         ifelse(volcano$score.anova.pvalue <= 0.05 & volcano$q.a <= 0.1,"10%FDR",
                                "not")))
#######plot####
library (ggplot2)

plot<- ggplot(volcano, aes(x= estimate, y=-log10(score.anova.pvalue))) +
              geom_point(aes(size = as.character(label),
                 alpha = as.character(label), 
                 fill = as.character(label)), 
                 pch=21)+
                 xlab("estimate") +
                 ylab("-log10 p-value") +
                 scale_size_manual(name = "",
                    values = c("5%FDR" = 1.5, "10%FDR" = 1, "not" = 0.25))+
                scale_alpha_manual(name = "",
                     values = c("5%FDR" = 1, "10%FDR" = 1, "not" = 0.25))+
                scale_fill_manual(name = "",
                    values = c("5%FDR" = "#EF1665", "10%FDR" = "#166DBA", "not"="#57636e"))+
                theme_bw() +
                xlim(-25000, 25000)+
                geom_text(aes(label = ifelse(q.a <= 0.1 & estimate < -2500 | q.a <= 0.1 & estimate > 2500, as.character(RNA_big.symbol), ""),
                hjust=0, vjust=0, color = "black")) 


#######save plot#####
pdf("RNA seq PRS-HOMAIR D14.pdf")
plot
dev.off()

write.table(volcano, file="RNA seq PRS-HOMAIR D14.csv", sep=",", col.names = NA, qmethod = "double")


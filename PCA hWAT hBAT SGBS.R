setwd("/Users/sophiestrobel/Desktop/LP/data raw/LP data import")
####load raw data files#####
input<-read.csv(file = 'LP_hWAT_hBAT_SGBS_data.csv')
input<-as.data.frame(input[,-1])
############subset data#####
nonffa<-subset(input, input$FFA == "0")

######remove####
set<-c("SGBS_14_0_3_BR00112243", "hWAT_3_0_4_BR00112227", "SGBS_0_0_5_BR00112231")
nonffa.1<-subset(nonffa, ! nonffa$newcol %in% set)
######split data in input and meta data#####
input<- t(nonffa.1[, 9:3010])
#######PCA#####
########remove features with 0´s or NA´s####
mean<-rowMeans(input)
var<-as.data.frame(cbind(mean, input))

var.1<-subset(var, var$mean != "0")
var.1a<- var.1[complete.cases(var.1), ]

variables<-as.data.frame(t(var.1a[,2:74]))

pca<-prcomp(variables, scale = FALSE)
#######visualise####
label<-as.data.frame(nonffa.1[, c("newcol","batch" ,"patientID", "day")])
######extract PC´s###
plot<-pca[["x"]]
######merge PC data und label####
plot.1<-cbind(label, plot)
######plot#####
####by cell type####
celltype <- ggplot(plot.1, aes(x=PC1 ,y=PC2)) + 
            geom_point(aes(fill = patientID,
                           alpha = as.character(day)),
                           size = 2.5,
                           pch= 21)+
            scale_fill_manual(name= "none",
                    values= c("hWAT" ="#25ad1d", "hBAT" = "#ffc933", "SGBS" = "#3363ff"))+
            scale_alpha_manual(name= "none",
                     values= c("0" = 0.2, "3" = 0.5, "8" = 0.7, "14" = 1))+
            labs(x='PC1', y="PC2",title='PCA analysis') +
            ylim(-20, 20)+
            xlim(-25, 35)+
            theme_bw() 

####by batch####
batch <- ggplot(plot.1, aes(x=PC1 ,y=PC2)) + 
            geom_point(aes(fill = as.character(batch),
                        alpha = as.character(day)),
                        size = 2.5,
                        pch= 21)+
            scale_alpha_manual(name= "none",
                     values= c("0" = 0.2, "3" = 0.5, "8" = 0.7, "14" = 1))+
            labs(x='PC1', y="PC2",title='PCA analysis') +
            ylim(-20, 20)+
            xlim(-25, 35)+
            theme_bw() 


#####plot and save next to each other#####
library(ggpubr)
pdf("PCA hWAT hBAT SGBS.pdf")
ggarrange(celltype, batch, nrow= 1, ncol = 2)

dev.off()
##########


#####Visualising LP profiles using ComplexHeatmp########
####set working directory#####
setwd("/Users/sophiestrobel/Desktop/LP/network profiles/gene feature connections data/")
input<-read.csv(file = "profile_TIMM22.csv")
input<-input[,-1]

######optional####
####remove redundant features and keep features with largest effect size#####
########structure feature string#######
library(stringr)
input$variable<-str_remove(input$variable, "[0-9][0-9][_][0-9][0-9]")
input$variable<-str_remove(input$variable, "[0-9][_][0-9][0-9]")
input$variable<-str_remove(input$variable, "[0-9][_][0-9]")
input$variable<-str_remove(input$variable, "[_][0-9][_]")
input$variable<-str_remove(input$variable, "[0-9]*[0-9]")
#########Granularity features#######
granularity<-input[grepl("_Granularity_", input$raw_features),]
granularity$variable<-granularity$raw_features

input.a<-input[!grepl("_Granularity_", input$raw_features),]

input.b<-rbind(granularity, input.a)
#######order feature based on feature string and effect size######
input_ordered<-input.b[order(input.b[,"variable"], -abs(input.b[,"beta"])),]

#######redundant features are sorted in groups, where feature with highest effect size is on top
#######keep redundant feature with highest effetct size
nonredundant_features<-input_ordered[!duplicated(input_ordered$variable), ]

##########################################
######categories LP features in location in the cell, features channel and features measrurment######
library(dplyr)
pie_1 <- nonredundant_features %>% dplyr::mutate(feature_location = 
                                                               ifelse(grepl("Cells_", nonredundant_features$variable),'Cells', 
                                                                      ifelse(grepl("Cytoplasm_", nonredundant_feature$variable),'Cytoplasm', 
                                                                             ifelse(grepl("Nuclei", nonredundant_features$variable),'Nuclei',
                                                                                    'other_feature'))))

pie_2 <- pie_1 %>% dplyr::mutate(feature_org = 
                                   ifelse(grepl("BODIPY", pie_1$variable) & ! grepl("AGP", pie_1$variable) & ! grepl("Mito", pie_1$variable) & ! grepl("DNA", pie_1$variable), 'BODIPY', 
                                          ifelse(grepl("BODIPY", pie_1$variable) & grepl("AGP", pie_1$variable), 'cor',      
                                                 ifelse(grepl("BODIPY", pie_1$variable) & grepl("Mito", pie_1$variable), 'cor',
                                                        ifelse(grepl("BODIPY", pie_1$variable) & grepl("DNA", pie_1$variable),'cor',
                                                               ifelse(grepl("Mito", pie_1$variable) & ! grepl("AGP", pie_1$variable) & ! grepl("BODIPY", pie_1$variable) & ! grepl("DNA", pie_1$variable), 'Mito',       
                                                                      ifelse(grepl("Mito", pie_1$variable)  & grepl("AGP", pie_1$variable) ,'cor', 
                                                                             ifelse(grepl("Mito", pie_1$variable)  & grepl("DNA", pie_1$variable) ,'cor', 
                                                                                    ifelse(grepl("AGP", pie_1$variable) & ! grepl("BODIPY", pie_1$variable) & ! grepl("Mito", pie_1$variable) & ! grepl("DNA", pie_1$variable), 'AGP',                    
                                                                                           ifelse(grepl("AGP", pie_1$variable)  & grepl("DNA", pie_1$variable) ,'cor',
                                                                                                  ifelse(grepl("DNA", pie_1$variable) & ! grepl("AGP", pie_1$variable) & ! grepl("Mito", pie_1$variable) & ! grepl("BODIPY", pie_1$variable), 'DNA', 
                                                                                                         'other_feature')))))))))))


pie_3 <- pie_2 %>% dplyr::mutate(feature_measure = 
                                   ifelse(grepl("Texture", pie_2$variable),'Texture', 
                                          ifelse(grepl("Intensity", pie_2$variable),'Intensity', 
                                                 ifelse(grepl("Granularity", pie_2$variable),'Granularity',
                                                        'other_feature'))))




##########generate Heatmap########
library(ComplexHeatmap)
library(circlize)

#####effect sizes###
matrix<-as.data.frame(pie_3[,"beta"])
rownames(matrix)<-pie_3$raw_features
colnames(matrix)<-"TIMM22 profile"
matrix[,] <- sapply(sapply(matrix[,], as.character),as.numeric)

#####coloring####
features<-as.character(pie_3$raw_features)
feature_org<-pie_3$feature_org
feature_measure<-pie_3$feature_measure

meta<-as.data.frame(cbind(feature_org,feature_measure, features))

#######define color palette#####
f1 = colorRamp2(seq(min(matrix), max(matrix), length = 3),  c("blue", "#EEEEEE","red"))
#f2 = colorRamp2(seq(min(matrix), max(matrix), length = 3), c("#0DB0BB","white",  "#E9511E"))

#####plot####
profile = Heatmap(as.matrix(matrix),
            col = f2,
            name = "beta",
            show_row_names = T,
            row_names_gp = gpar(col = "black", fontsize = 2),
            right_annotation = rowAnnotation(channel =  meta$feature_org,
                                             measurement = meta$feature_measure,
                                             col= list(channel = c("Mito" = "#f56464" ,
                                                                   "AGP" = "#f2cf41",
                                                                   "BODIPY" = "#4dac26",
                                                                   "DNA" = "#3b8fe3",
                                                                   "other_feature" = "#959ca3",
                                                                   "cor" = "#d7dce0"),
                                                       measurement =c("Intensity" = "#b5bac7" ,
                                                                     "Granularity" = "#626a80",
                                                                     "Texture" = "#343f5c",
                                                                     "other_feature" = "#e4e7f0"))))

pdf("TIMM22 profile.pdf")
profile
dev.off()

#####extract data for LP profile of specific genes####
#####set directory####
setwd("/Users/slaber/Desktop/")
#####load input files#####
sig<-read.csv(file="Cells_significants_LMM_5FDR.csv")

significant_features<-sig[,-c(1, 2)]

###load csv file with gene annotation##
gene_symbol<-read.csv(file="gene_symbol.csv")
gene_symbol<-gene_symbol[,-1]
colnames(gene_symbol)<-c("gene", "EN_number")

######merge gene annotation and LMM data###
all_data<-merge(significant_features, gene_symbol, "EN_number", all.x =T, sort=F)

######remove redundant features###
raw_features<-as.character(all_data$variable)
significant_features<-cbind(all_data, raw_features)

########structure feature string####
library(stringr)
significant_features$variable<-str_remove(significant_features$variable, "[0-9][0-9][_][0-9][0-9]")
significant_features$variable<-str_remove(significant_features$variable, "[0-9][_][0-9][0-9]")
significant_features$variable<-str_remove(significant_features$variable, "[0-9][_][0-9]")
significant_features$variable<-str_remove(significant_features$variable, "[_][0-9][_]")
significant_features$variable<-str_remove(significant_features$variable, "[0-9]*[0-9]")

granularity<-significant_features[grepl( "_Granularity_" , significant_features$raw_features ) ,  ]
granularity$variable<-granularity$raw_features

significant_features<-significant_features[!grepl( "_Granularity_" , significant_features$raw_features ) ,  ]

significant_features<-rbind(significant_features, granularity)
#####order feature based on feature string and effect size#######
significant_features_ordered<-significant_features[order(significant_features[,"variable"], -abs(significant_features[,"beta"])),]

significant_features_ordered<-unite(significant_features_ordered, newcol, c("EN_number", "variable"), remove=FALSE)
####redundant features are sorted in groups, where feature with highest effect size is on the top######
###keep redundant feature gene association with highest effect size######
library(tidyr)
nonredundant_significant_features<-significant_features_ordered[!duplicated(significant_features_ordered$newcol), ]

#####annotaed features based on compartment, channel, measurment#####
pie_1 <- nonredundant_significant_features %>% dplyr::mutate(feature_location = 
                                                               ifelse(grepl("Cells_", nonredundant_significant_features$variable),'Cells', 
                                                                      ifelse(grepl("Cytoplasm_", nonredundant_significant_features$variable),'Cytoplasm', 
                                                                             ifelse(grepl("Nuclei", nonredundant_significant_features$variable),'Nuclei',
                                                                                    'other_feature'))))

pie_2 <- pie_1 %>% dplyr::mutate(feature_org = 
                                   ifelse(grepl("BODIPY", pie_1$variable) & ! grepl("AGP", pie_1$variable) & ! grepl("Mito", pie_1$variable) & ! grepl("DNA", pie_1$variable), 'BODIPY', 
                                          ifelse(grepl("BODIPY", pie_1$variable) & grepl("AGP", pie_1$variable), 'BODIPY_AGP',      
                                                 ifelse(grepl("BODIPY", pie_1$variable) & grepl("Mito", pie_1$variable), 'BODIPY_Mito',
                                                        ifelse(grepl("BODIPY", pie_1$variable) & grepl("DNA", pie_1$variable),'BODIPY_DNA',
                                                               ifelse(grepl("Mito", pie_1$variable) & ! grepl("AGP", pie_1$variable) & ! grepl("BODIPY", pie_1$variable) & ! grepl("DNA", pie_1$variable), 'Mito',       
                                                                      ifelse(grepl("Mito", pie_1$variable)  & grepl("AGP", pie_1$variable) ,'Mito_AGP', 
                                                                             ifelse(grepl("Mito", pie_1$variable)  & grepl("DNA", pie_1$variable) ,'Mito_DNA', 
                                                                                    ifelse(grepl("AGP", pie_1$variable) & ! grepl("BODIPY", pie_1$variable) & ! grepl("Mito", pie_1$variable) & ! grepl("DNA", pie_1$variable), 'AGP',                    
                                                                                           ifelse(grepl("AGP", pie_1$variable)  & grepl("DNA", pie_1$variable) ,'AGP_DNA',
                                                                                                  ifelse(grepl("DNA", pie_1$variable) & ! grepl("AGP", pie_1$variable) & ! grepl("Mito", pie_1$variable) & ! grepl("BODIPY", pie_1$variable), 'DNA', 
                                                                                                         'other_feature')))))))))))


pie_3 <- pie_2 %>% dplyr::mutate(feature_measure = 
                                   ifelse(grepl("Texture", pie_2$variable),'Texture', 
                                          ifelse(grepl("Intensity", pie_2$variable),'Intensity', 
                                                 ifelse(grepl("Granularity", pie_2$variable),'Granularity',
                                                        'other_feature'))))
####select gene for LP profile####
#### TIMM22 = ENSG00000177370
gene_profile<-subset(pie_3, pie_3$EN_number == "ENSG00000177370")

write.table(gene_profile, file="TIMM22profile.csv", sep=",", col.names = NA, qmethod = "double")

######

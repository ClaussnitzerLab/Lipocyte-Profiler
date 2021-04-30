#######piechart of significant features of morphological profiles####
##set working directory###
#setwd("/Users/sophiestrobel/Desktop/LP/data raw/LP data import")
####input data you want to filter####
####eg. df = nonredundant_significant_features##

#####assign categories to LP features by compartment, channel, measurement###
####compartment###
pie_1 <- nonredundant_significant_features %>% dplyr::mutate(feature_location = 
                                                               ifelse(grepl("Cells_", nonredundant_significant_features$features),'Cells', 
                                                                      ifelse(grepl("Cytoplasm_", nonredundant_significant_features$features),'Cytoplasm', 
                                                                             ifelse(grepl("Nuclei", nonredundant_significant_features$features),'Nuclei',
                                                                                    'other_feature'))))
####channel###
pie_2 <- pie_1 %>% dplyr::mutate(feature_org = 
                                   ifelse(grepl("BODIPY", pie_1$features) & ! grepl("AGP", pie_1$features) & ! grepl("Mito", pie_1$features) & ! grepl("DNA", pie_1$features), 'BODIPY', 
                                          ifelse(grepl("BODIPY", pie_1$features) & grepl("AGP", pie_1$features), 'BODIPY_AGP',      
                                                 ifelse(grepl("BODIPY", pie_1$features) & grepl("Mito", pie_1$features), 'BODIPY_Mito',
                                                        ifelse(grepl("BODIPY", pie_1$features) & grepl("DNA", pie_1$features),'BODIPY_DNA',
                                                               ifelse(grepl("Mito", pie_1$features) & ! grepl("AGP", pie_1$features) & ! grepl("BODIPY", pie_1$features) & ! grepl("DNA", pie_1$features), 'Mito',       
                                                                      ifelse(grepl("Mito", pie_1$features)  & grepl("AGP", pie_1$features) ,'Mito_AGP', 
                                                                             ifelse(grepl("Mito", pie_1$features)  & grepl("DNA", pie_1$features) ,'Mito_DNA', 
                                                                                    ifelse(grepl("AGP", pie_1$features) & ! grepl("BODIPY", pie_1$features) & ! grepl("Mito", pie_1$features) & ! grepl("DNA", pie_1$features), 'AGP',                    
                                                                                           ifelse(grepl("AGP", pie_1$features)  & grepl("DNA", pie_1$features) ,'AGP_DNA',
                                                                                                  ifelse(grepl("DNA", pie_1$features) & ! grepl("AGP", pie_1$features) & ! grepl("Mito", pie_1$features) & ! grepl("BODIPY", pie_1$features), 'DNA', 
                                                                                                         'other_feature')))))))))))

####measurement###
pie_3 <- pie_2 %>% dplyr::mutate(feature_measure = 
                                   ifelse(grepl("Texture", pie_2$features),'Texture', 
                                          ifelse(grepl("Intensity", pie_2$features),'Intensity', 
                                                 ifelse(grepl("Granularity", pie_2$features),'Granularity',
                                                        'other_feature'))))

###plot pie charts###
library(ggplot2)
library(moonBook)
library(webr)

###plot compartments and channels
PieDonut(pie_3,aes(pies= feature_location , donuts =  feature_org))

###plot channels and measurement
PieDonut(pie_3,aes(pies= feature_org , donuts =  feature_measure))

####################



####
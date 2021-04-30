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

pie_4 <- pie_3 %>% dplyr::mutate(significant = 
                                   ifelse(pie_3$q.a <= 0.05 & pie_3$estimate < 0,'negative', 
                                      ifelse(pie_3$q.a <= 0.05 & pie_3$estimate > 0,'positive', 
                                                        'other_feature')))

###plot heatmap of effect size###
library(ComplexHeatmap)
library(circlize)
#####transform effect size into matrix###
matrix<-as.matrix(pie_4$score.eta_sq)
colnames(mat)<-"eta_sq"
#####define color palette####
f2 = colorRamp2(seq(min(matrix), max(matrix), length = 6), c( "#fafafa","#f2d38a","#f2bf80", "#f29180","#f29180","#f01d1d"))
###plot###
plot = Heatmap(matrix,
             col = f2,
             name = "effect size",
             show_row_names = F, 
             #row_names_gp = gpar(col = "black", fontsize = 0.5),
             column_names_gp = gpar(col = "black", fontsize = 2),
             right_annotation = rowAnnotation(channel =  pie_4$feature_org,
                                              measurment = pie_4$feature_measure,
                                              significance = pie_4$significant,
                                              col= list(channel = c("Mito" = "#f56464" ,
                                                                    "AGP" = "#f2cf41",
                                                                    "BODIPY" = "#4dac26",
                                                                    "DNA" = "#3b8fe3",
                                                                    "other_feature" = "#959ca3",
                                                                    "AGP_DNA" = "#d7dce0",
                                                                    "BODIPY_AGP" = "#d7dce0",
                                                                    "Mito_AGP" = "#d7dce0",
                                                                    "BODIPY_Mito" = "#d7dce0",
                                                                    "BODIPY_DNA" = "#d7dce0",
                                                                     "Mito_DNA" = "#d7dce0"),
                                              measurment =c("Intensity" = "#b5bac7" ,
                                                            "Granularity" = "#626a80",
                                                            "Texture" = "#343f5c",
                                                            "other_feature" = "#e4e7f0"), 
                                              significance = c("negative"= "#9BC7E4", 
                                                               "positive"= "#B02636", 
                                                               "other_feature" = "white"))))

####save####
pdf("heatmap_profile_comparison.pdf")
plot
dev.off()
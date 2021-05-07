#####visualise newtork using igraph####
#####set directory####
setwd("/Users/slaber/Desktop/")
#####load input files#####
sig0.1<-read.csv(file="Cells_significants0.001.csv")

significant_features<-sig0.1[,-c(1, 2)]

###load csv file with gene annotation##
gene_symbol<-read.csv(file="gene_symbol.csv")
gene_symbol<-gene_symbol[,-1]
colnames(gene_symbol)<-c("gene", "EN_number")

######merge gene annotation and LMM data###
all_data<-merge(significant_features, gene_symbol, "EN_number", all.x =T, sort=F)

######remove redundant features###
raw_features<-as.character(all_data$variable)
significant_features<-cbind(all_data, raw_features)

#significant_features<-subset(significant_features,significant_features$variable %in% raw )
#significant_features<-cbind(raw_features,sig)

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

######generate network####
library(igraph)
library(dplyr)

#####filter data based on genes with min. xx connections, or features with min.xx connections###
#newd <-  pie_3 %>% group_by(EN_number) %>% filter(n() >= 5) 
#newd <-  pie_3 %>% group_by(raw_features) %>% filter(n() >= 5) 


######label nodes of the network in origin = LP features and destination = gene####
origin<-as.character(pie_3$raw_features)
destination<-as.character(pie_3$EN_number)
links <- data.frame(origin, destination)

######define genes and features name###
unique_origin<-unique(as.character(origin))
unique_dest<-unique(as.character(destination))

name<-c(as.character(unique_origin), as.character(unique_dest))

#####adj. numbers based on unique_origin and unique_dest####
carac<-c(rep("AP", 5), rep("RNA",36 ))

#####define categories of nodes####
nodes<-data.frame(name, carac)

nodes.1<- nodes %>% dplyr::mutate(feature_org = 
                                   ifelse(grepl("BODIPY", nodes$name) & ! grepl("AGP", nodes$name) & ! grepl("Mito", nodes$name) & ! grepl("DNA", nodes$name), 'BODIPY', 
                                          ifelse(grepl("BODIPY", nodes$name) & grepl("AGP", nodes$name), 'BODIPY_AGP',      
                                                 ifelse(grepl("BODIPY", nodes$name) & grepl("Mito", nodes$name), 'BODIPY_Mito',
                                                        ifelse(grepl("BODIPY", nodes$name) & grepl("DNA", nodes$name),'BODIPY_DNA',
                                                               ifelse(grepl("Mito", nodes$name) & ! grepl("AGP", nodes$name) & ! grepl("BODIPY", nodes$name) & ! grepl("DNA", nodes$name), 'Mito',       
                                                                      ifelse(grepl("Mito", nodes$name)  & grepl("AGP", nodes$name) ,'Mito_AGP', 
                                                                             ifelse(grepl("Mito", nodes$name)  & grepl("DNA", nodes$name) ,'Mito_DNA', 
                                                                                    ifelse(grepl("AGP", nodes$name) & ! grepl("BODIPY", nodes$name) & ! grepl("Mito", nodes$name) & ! grepl("DNA", nodes$name), 'AGP',                    
                                                                                           ifelse(grepl("AGP", nodes$name)  & grepl("DNA", nodes$name) ,'AGP_DNA',
                                                                                                  ifelse(grepl("DNA", nodes$name) & ! grepl("AGP", nodes$name) & ! grepl("Mito", nodes$name) & ! grepl("BODIPY", nodes$name), 'DNA',
                                                                                                         "other_feature")))))))))))
                                                                                                         
####Transform input data in a adjacency matrix#####
network <- graph_from_data_frame(d=links, vertices=nodes.1, directed=F) 

####set new margins to limit whitespace in plot####
par(mar=rep(.1, 4))

#####adj vertex and edge size based on number of connections####
deg <- degree(network, mode="all")

####define color and shape##
###circle = LP features###
###square = genes####
coul<-c("#F2CF41", "#959CA3", "#4DAC26", "#959CA3","#959CA3", "#959CA3","#36CEF7", "#F56464", "#959CA3", "#959CA3","#959CA3")
shape<-c("circle","square")

my_color <- coul[as.numeric(as.factor(V(network)$feature_org))]
my_shape <- shape[as.numeric(as.factor(V(network)$carac))]

####make the plot####
##adj deg*X based on connections###
##adj layout###
plot(network, vertex.color=adjustcolor(my_color, alpha.f = .75), 
              vertex.shape=my_shape, vertex.size=deg*1, 
              layout=layout.sphere, 
              vertex.label.cex=deg*0.5, 
              edge.width=1, 
              edge.color="#f0ebeb" )


pdf("newtork.pdf")
plot
dev.off()

#########

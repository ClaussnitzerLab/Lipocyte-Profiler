#####remove redundant features #####
####keep redundant feature with largest effect size#####
####input data you want to filter####
####eg. df = significant features##
raw_features<-as.character(significant_features$features)

significant_features<-cbind(significant_features, raw_features)
########structure feature string#######
library(stringr)
significant_features$features<-str_remove(significant_features$features, "[0-9][0-9][_][0-9][0-9]")
significant_features$features<-str_remove(significant_features$features, "[0-9][_][0-9][0-9]")
significant_features$features<-str_remove(significant_features$features, "[0-9][_][0-9]")
significant_features$features<-str_remove(significant_features$features, "[_][0-9][_]")
significant_features$features<-str_remove(significant_features$features, "[0-9]*[0-9]")
#########Granularity features#######
granularity<-significant_features[grepl("_Granularity_", significant_features$raw_features),]
granularity$features<-granularity$raw_features

significant_features.a<-significant_features[!grepl("_Granularity_", significant_features$raw_features),]

significant_features.b<-rbind(granularity, significant_features.a)
#######order feature based on feature string and effect size######
significant_features_ordered<-significant_features.b[order(significant_features.b[,"features"], -abs(significant_features.b[,"estimate"])),]

#######redundant features are sorted in groups, where feature with highest effect size is on top
#######keep redundant feature with highest effetct size
nonredundant_significant_features<-significant_features_ordered[!duplicated(significant_features_ordered$features), ]

#########

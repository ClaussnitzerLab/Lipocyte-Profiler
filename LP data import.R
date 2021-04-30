###Lipocyte Profiler on AMSCs of two adipose depots, subcutaneous and visceral, across 8 batches with and without FFA treatment.
#scripte for loading and merging of data inputs#
################################################
#####set working directory#####
setwd("/Users/sophiestrobel/Desktop/LP/data raw/LP data import")
####load raw data files#####
patient<-read.csv(file = 'patients.csv')
names(patient)[names(patient) == "ID"] <- "patientID"

layouts<-read.csv(file = 'layouts.csv')

plates<-read.csv(file = 'plates.csv')

#####file including data of sample handling
nerd_data<-read.csv(file = 'nerdData.csv')

meta_geno<-read.csv(file = "Metadata - Geno.csv")

PRS<-read.csv(file = "PRS_all_0403.csv")

merged_meta<-merge(x = meta_geno, y = PRS, by = "patientID", sort = FALSE)

meta<-read.csv(file = "MOBB_Phenotype.csv")
meta<-meta[,1:10]

####import plate normalised LP data##########
new.normalised<-read.csv(file = 'features_replicates_median_plate_norm.csv')

####merge meta data and feature data#######
#####################################################################
###use same column names for dataframes patient and nerd_data

names(nerd_data)[names(nerd_data) == "depot"] <- "cellType"

####use just m239 of batch 4, remove m239 of batch 6###
nerd_data.1<-subset(nerd_data, nerd_data$patientID != "m239_6")
nerd_data.1$patientID <- gsub("m239_4", "m239", nerd_data.1$patientID)

###########merge meta data######
merged_meta_all<-merge(x = merged_meta, y = meta, by = "patientID", sort = FALSE)

merged_meta_patient<-merge(x = patient, y = merged_meta_all, by = "patientID", all.x = T, sort = FALSE)
merged_meta_patient<-merged_meta_patient[!duplicated(merged_meta_patient$patientID) , ]

######structure layout, create new column named "column" for well position localisation###
layout<-layouts[!(is.na(layouts$patientID) | layouts$patientID==""), ]

library(stringr)
well<-str_split_fixed(layout$well, "[[:digit:]]", 2)
layout<-cbind(layout, well[,-1])
colnames(layout)<-c("plate", "well", "cellType", "patientID", "FFA", "column")

######
layout$cellType <- gsub("0", "sc", layout$cellType)
layout$cellType <- gsub("1", "vc", layout$cellType)

#######merge layout and plate -> merged layout
merged_layout <- merge(x = layout, y = plates, by = "plate", sort = FALSE)
merged_layout<-merged_layout[,-9]

#####median feature data aggregate by replicate
#####overwrite ER into BODIPY
colnames(new.normalised) <- gsub("_ER", "_BODIPY", colnames(new.normalised))

##########use just m239 of batch 4, remove m239 of batch 6
median_data<-subset(new.normalised, new.normalised$patientID != "m239_6")
median_data$patientID <- gsub("m239_4", "m239", median_data$patientID)

#######merge merged_layout and merged_plate by creating new column
library(tidyr)
merged_layout.1<-unite(merged_layout, newcol, c("patientID", "cellType"), remove=FALSE)

###############split dataframe into to batches sc anbd vc
sc_patient<-subset(merged_meta_patient, merged_meta_patient$sc == "TRUE")
cellType<-rep("sc", 33)
sc_patient<-cbind(sc_patient, cellType)

vc_patient<-subset(merged_meta_patient, merged_meta_patient$vc == "TRUE")
cellType<-rep("vc", 39)
vc_patient<-cbind(vc_patient, cellType)

patient<-rbind(sc_patient[, -c(2:3)], vc_patient[, -c(2:3)])

merged_patient<-unite(patient, newcol, c("patientID", "cellType"), remove=FALSE)

##########merge nerd data######
nerd_data.2<-unite(nerd_data.1, newcol, c("patientID", "cellType"), remove=FALSE)

merged_metadata_patient <- merge(x = merged_patient, y = nerd_data.2, by = "newcol",  all.x= T, sort = FALSE)
merged_metadata_patient.1<-merged_metadata_patient[, -c(27:28, 33:49)]

##############

merged_metadata.1<-subset(merged_layout, merged_layout$patientID != "m239_6")
merged_metadata.1$patientID <- gsub("m239_4", "m239", merged_metadata.1$patientID)
merged_metadata.1a<-unite(merged_metadata.1, newcol, c("patientID", "cellType"), remove=FALSE)

merged_meta_layout_patient<- merge(x=merged_metadata_patient.1 , y=merged_metadata.1a  , by = "newcol", all.x = T, all.y=T, sort = FALSE)

names(merged_meta_layout_patient)[names(merged_meta_layout_patient) == "patientID.x"] <- "patientID"
names(merged_meta_layout_patient)[names(merged_meta_layout_patient) == "cellType.x"] <- "cellType"

###################merge meta data and feature data by creating new column
library(tidyverse)

merged_meta_layout_patient<-unite(merged_meta_layout_patient, newcol, c("patientID", "day", "cellType", "FFA"), remove=FALSE)

merged_meta_layout_patient_a<-merged_meta_layout_patient[!duplicated(merged_meta_layout_patient$newcol), ]
#######################
median_data.1<-unite(median_data, newcol, c("patientID", "day", "cellType", "FFA"), remove=FALSE)

########
combined_data <- merge(x=merged_meta_layout_patient_a , y=median_data.1  , by = "newcol", all.x = F, all.y=T, sort = FALSE)

combined_data<-combined_data[,-c(31,37:40)]
###################structure and rename columns
names(combined_data)[names(combined_data) == "patientID.x"] <- "patientID"
names(combined_data)[names(combined_data) == "cellType.x"] <- "cellType"
names(combined_data)[names(combined_data) == "day.x"] <- "day"
names(combined_data)[names(combined_data) == "FFA.x"] <- "FFA"


#########delete rosi and isoproternol treatment plates D14#####
treatment_plates<-c("BR00112242",
                    "BR00112241",
                    "BR00112240",
                    "BR00112238",
                    "BR00112236",
                    "BR00112233",
                    "BR00112232",
                    "BR00112230"
)

combined_data<-subset(combined_data, ! combined_data$Metadata_Plate %in% treatment_plates)

write.table(combined_data, file = "LP_data.csv", sep = ",", col.names = NA,
            qmethod = "double")


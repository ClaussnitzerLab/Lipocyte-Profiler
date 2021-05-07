###Lipocyte Profiler on hWAT hBAT SGBS, across 3 batches with and without FFA treatment.
#scripte for loading and merging of data inputs#
################################################
#####set working directory#####
setwd("/Users/sophiestrobel/Desktop/LP/data raw/LP data import")
####load raw data files#####
layouts<-read.csv(file = 'layouts.csv')
plates<-read.csv(file = 'plates.csv')
#####LP data hWAT hBAT SGBS####
new.normalised<-read.csv(file = 'features_replicates_median_plate_norm_nonPACs.csv')
#####merge data files####
layout<-layouts[!(is.na(layouts$patientID) | layouts$patientID==""), ]

library(stringr)
well<-str_split_fixed(layout$well, "[[:digit:]]", 2)
layout<-cbind(layout, well[,-1])
colnames(layout)<-c("plate", "well", "cellType", "patientID", "FFA", "column")

#####structure data###
layout$cellType <- gsub("0", "sc", layout$cellType)
layout$cellType <- gsub("1", "vc", layout$cellType)

#######merge layout and plate -> merged layout
merged_layout <- merge(x = layout, y = plates, by = "plate", sort = FALSE)
merged_layout<-merged_layout[,-9]

#####overwrite ER into BODIPY####
colnames(new.normalised) <- gsub("_ER", "_BODIPY", colnames(new.normalised))

#######merge merged_layout and merged_plate by creating new column
library(tidyr)
new.normalised.1<-unite(new.normalised, newcol, c("patientID", "day", "FFA", "Metadata_Column", "Metadata_Plate"), remove=FALSE)

merged_layout.1<-unite(merged_layout, newcol, c("patientID", "day", "FFA", "column", "plate"), remove=FALSE)
merged_layout.1a<-merged_layout.1[!duplicated(merged_layout.1$newcol) , ]

merged_metadata <- merge(x = merged_layout.1a, y = new.normalised.1, by = "newcol",  all.y= T, sort = FALSE)

combined_data<-merged_metadata[, -c(2,4, 10:12, 14)]


#########structure and rename columns#######
names(combined_data)[names(combined_data) == "patientID.x"] <- "patientID"
names(combined_data)[names(combined_data) == "day.x"] <- "day"
names(combined_data)[names(combined_data) == "FFA.x"] <- "FFA"

combined_data$patientID <- gsub("hWAT-1", "hWAT", combined_data$patientID)
combined_data$patientID <- gsub("hWAT-2", "hWAT", combined_data$patientID)
combined_data$patientID <- gsub("hBAT-1", "hBAT", combined_data$patientID)
combined_data$patientID <- gsub("hBAT-2", "hBAT", combined_data$patientID)

#########delete rosi and isoproternol treatment plates D14#####
treatment_plates<-c("BR00112242",
                    "BR00112241",
                    "BR00112240",
                    "BR00112238",
                    "BR00112236",
                    "BR00112233",
                    "BR00112232",
                    "BR00112230")

combined_data<-subset(combined_data, ! combined_data$Metadata_Plate %in% treatment_plates)

###########################
write.table(combined_data, file = "LP_hWAT_hBAT_SGBS_data.csv", sep = ",", col.names = NA,
            qmethod = "double")



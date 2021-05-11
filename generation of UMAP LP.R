######UMAP of LP data#####
####set working directory#####
setwd("/Users/sophiestrobel/Desktop/LP/data raw/LP data import")

###load LP data set#####
combined_data<-read.csv(file = "LP_data.csv")

###subset dataset#####
#ffa<-subset(combined_data, combined_data$FFA == "1")
nonffa<-subset(combined_data, combined_data$FFA == "0")

########UMAP######
meta<-nonffa[, 1:36]
input<-nonffa[, 37:3038]

######tranform LP features into numeric#####
input[,] <- sapply(sapply(input[,], as.character),as.numeric) 

######remove blocklisted features and features with 0´s or NA´s across data###
variables<-input

variables[,] <- sapply(sapply(variables[,], as.numeric), as.numeric) 
back<-as.data.frame(t(variables))

var.selection.1<-back[!grepl( "_Costes_" , rownames(back) ) ,  ]
var.selection.2<-var.selection.1[!grepl( "_Manders_" , rownames(var.selection.1)) ,  ]
var.selection.3<-var.selection.2[!grepl( "_RWC_" , rownames(var.selection.2)) ,  ]
var.selection.5<-var.selection.3[!grepl( "SmallBODIPY" , rownames(var.selection.3) ) ,  ]

mean<-rowMeans(var.selection.5)
var<-cbind(mean, var.selection.5)

var.1<-subset(var, var$mean != "0")
var.1a<- var.1[complete.cases(var.1), ]

variables<-as.data.frame(t(var.1a[,2:244]))

#####merge meta and LP data###
#df<-cbind(meta, variables)

###########unsupervised clustering, distance measurment######
library(umap)
#####define parameters, change neighbors and min distance, first try with default settings#####
custom.config = umap.defaults
custom.config$n_neighbors = 21
custom.config$metric = "euclidean"
custom.config$min_dist= 0.05
custom.config$n_epoch = 200
custom.config$init= "spectral"

#####create umap matrix ####
umap <- umap(variables,custom.config)

####structure data for umap######
meta$day <- gsub("0", "1", meta$day)
meta$day <- gsub("3", "2", meta$day)
meta$day <- gsub("8", "3", meta$day)
meta$day <- gsub("14", "4",meta$day)


plot.data<-as.data.frame(unlist(umap[[1]]))
colnames(plot.data)<-c("x","y")
data_plot<-cbind(plot.data,meta)

data_plot$cellType<-factor(data_plot$cellType, level=c("sc", "vc"))

data_plot$cellType <- gsub("sc", "1", data_plot$cellType)
data_plot$cellType <- gsub("vc", "2", data_plot$cellType)

#####assign categories for day and depot###
library(dplyr)
data_plot_color <- data_plot %>% 
  dplyr::mutate(
    day_color = ifelse(
      data_plot$day == "4" & data_plot$cellType == "1", "day_14_sc", 
      ifelse(data_plot$day == "3" & data_plot$cellType == "1", "day_8_sc",
             ifelse(data_plot$day == "2" & data_plot$cellType == "1", "day_3_sc",
                    ifelse(data_plot$day == "1" & data_plot$cellType == "1", "day_0_sc", 
                           ifelse(data_plot$day == "4" & data_plot$cellType == "2", "day_14_vc", 
                                  ifelse(data_plot$day == "3" & data_plot$cellType == "2", "day_8_vc",
                                         ifelse(data_plot$day == "2" & data_plot$cellType == "2", "day_3_vc",
                                                ifelse(data_plot$day == "1" & data_plot$cellType == "2", "day_0_vc","others")))))))))


data_plot_color$day_color<-as.factor(data_plot_color$day_color)

#####plot####
library(ggplot2)
plot<- ggplot(data_plot_color, aes(x, y)) +
          geom_point(aes(color = paste(day_color),
                 alpha = as.numeric(paste(day)) ,
                 size =as.numeric(paste(day))),
                 pch= 20)+
            scale_color_manual(name = "Cell depot",
                    values = c("day_14_sc" = "#EF1665",
                               "day_8_sc" = "#ed6b9a",
                               "day_3_sc" = "#ed8aae",
                               "day_0_sc" = "#f0b4ca",
                               "day_14_vc" = "#166dba",
                               "day_8_vc" = "#3d7bb3",
                               "day_3_vc" = "#5c8bb5",
                               "day_0_vc" = "#849eb5")) +
            scale_alpha_continuous(name = "", 
                         range = c(0.5,0.8), 
                         guide = 'none')+
            scale_size_continuous(name = "", 
                        range = c(1.5,4), 
                        guide = 'none') +
            theme_bw() +
            xlab("UMAP (x)") +
            ylab("UMAP (y)")

pdf("UMAP LP profile.pdf")
plot
dev.off()
######################

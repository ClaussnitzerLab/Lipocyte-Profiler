####variance explained analysis####
###variance component analysis on adipocyte morphological and cellular traits 
#across 65 donor-derived AMSCs (D0-D14) to assess the contribution of intrinsic genetic variation compared to 
#the contribution of other possible confounding factors 
#such as batch, adipose depot, T2D status, age, sex, BMI, 
#cell density, month/year of sampling, and passage numbers
########set working directory#####
setwd("/Users/sophiestrobel/Desktop/LP/data raw/")

######import data####
#combined_data<-read.csv(file = 'LP data import.csv')
#combined_data<-combined_data[,-1]

###subset data#####
nonffa<-subset(combined_data, combined_data$FFA == "0") 

####input data you want to analyis'######
LP_subset<-nonffa

#####subset in meta and LP feature data######
meta<-LP_subset[, 1:36]
input<-LP_subset[, 37:3038]

##########transfrom LP feature data to numeric######
variables<-input
variables[,] <- sapply(sapply(variables[,], as.character),as.numeric) 

######### remove backlisted features########
back<-as.data.frame(t(variables))

var.selection.1<-back[!grepl( "_Costes_" , rownames(back)) ,  ]
var.selection.2<-var.selection.1[!grepl( "_Manders_" , rownames(var.selection.1)) ,  ]
var.selection.3<-var.selection.2[!grepl( "_RWC_" , rownames(var.selection.2)) ,  ]
var.selection.4<-var.selection.3[!grepl( "SmallBODIPY" , rownames(var.selection.3) ) ,  ]

###########remove features with "0" values or NA´s across data set####
mean<-rowMeans(var.selection.4)
var<-cbind(mean, var.selection.4)
var.1<-subset(var, var$mean != "0")
var.1a<- var.1[complete.cases(var.1), ]

###########set n number#####
variables<-as.data.frame(t(var.1a[,2:244]))

#########merge meta and LP data again#####
df<-cbind(meta, variables)

#######variance component analysis######
####define charater of variables you want to loop through######
df$cellType<-as.factor(df$cellType)
df$patientID<-as.factor(df$patientID)
df$batch<-as.factor(df$batch)
df$sex<-as.factor(df$sex)
df$T2D<-as.factor(df$T2D)

#####choose density variable "Cells_Neighbors_PercentTouching_Adjacent" as one of the confounding variables#######
y<-as.data.frame(t(variables))
input.1<-y[!grepl( "Cells_Neighbors_PercentTouching_Adjacent" , rownames(y)) ,  ]

#########define LP features you want to include in lm#####
colvector<-rownames(input.1)
########lm function#####
out<-list()

for (i in colvector) {
  
  Fm1<-lm(df[,i]~
            cellType +
            batch +
            T2D +
            as.numeric(BMI) +
            sex +
            as.numeric(passaging) +
            as.numeric(month) +
            as.numeric(year) +
            as.numeric(age) +
            as.numeric(day) +
            as.numeric(Cells_Neighbors_PercentTouching_Adjacent )+ 
            patientID ,
            data = df)
  
  out[[i]] <- anova(Fm1)
  
}

###extract contribution of predicting variables to a dataframe########
contrib<-sapply(out,function(x) x[2])
contrib.1<-do.call(cbind, contrib)
####define row and column names#######
rownames(contrib.1)<-c(
    "adipose depot",
    "batch",
    "T2D",
    "BMI",
    "sex",
    "passaging",
    "month",
    "year",
    "age",
    "day",
    "Neighbors",
    "patientID" , 
    "residuals")

colnames(contrib.1)<-rownames(input.1)

#####calculate precentage of contribution#######
percentage<-as.data.frame(matrix(nrow=1, ncol=1, byrow=T, dimnames=NULL))

rowvector<-rownames(contrib.1)

for (i in colvector) {
  for (j in rowvector) {
    
    percentage[j,i]<-(contrib.1[j,i]/sum(contrib.1[,i]))*100
    
  }}

percentage_explained<-percentage[-1,-1]

####plot data in boxplot format#####
####structure output data for plotting#######
data<-percentage_explained

data.1<-tidyr::gather(data)

group<-as.factor(rownames(data))

####n = amount of features######
group<-rep(group, 2683)

df.1<-cbind(group, data.1)

########reorder for plot#####
df.1$group <- factor(df.1$group, levels = c( 
  "batch",
  "day",
  "adipose depot",
  "patientID" ,
  "BMI",
  "T2D",
  "sex",
  "age",
  "passaging",
  "month",
  "year",
  "Neighbors",
  "residuals"))


#########define color code########
######highlight donor#####
library(tidyr)
df_plot <- df.1 %>%
  dplyr::mutate(highlight = 
                  ifelse(grepl("patientID", df.1$group),'mark', 
                         'other_feature'))

########plot######
library(hrbrthemes)
library(ggplot2)

#########plot######
plot<-ggplot(df_plot,  aes(x=group, y=value, fill=highlight, alpha=highlight) )+ 
      geom_boxplot() +
      scale_fill_manual(values=c(  "#EF1665", "grey")) +
      scale_alpha_manual(values=c(0.7,0.3)) +
      theme_ipsum() +
      theme(legend.position = "none") +
      xlab("")+
      ylab("% explained variance")+
      ylim(0, 100)
      #coord_flip()

#######save data and plot#####
write.table(data, file = "table_variance_explained.txt", sep = ",", col.names = NA,
            qmethod = "double")

pdf("variance_explained.pdf")
plot
dev.off()

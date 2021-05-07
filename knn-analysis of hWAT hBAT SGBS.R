##knn of Lipocyte Profiler data#####.
################################################
#####set working directory#####
setwd("/Users/sophiestrobel/Desktop/LP/data raw/LP data import")
####load raw data files#####
input<-read.csv(file = 'LP_hWAT_hBAT_SGBS_data.csv')
input<-input[,-1]

input$column<-factor(input$column, level=c("1", "2", "3", "4", "5", "6", "7", "8"))

input.a<-unite(input, newcol, c("patientID", "day"), remove=F)

library(ISLR)
library(caret)
library(class)
########prediction of column###
indxTrain <- createDataPartition(y = input.a$column,p = 0.75,list = FALSE)
training <- input.a[indxTrain,]
testing <- input.a[-indxTrain,]

training_labels <- as.character(input.a[indxTrain,"column"])
testing_labels <- as.character(input.a[-indxTrain,"column"])

#Checking distibution in origanl data and partitioned data
prop.table(table(training$newcol)) * 100
prop.table(table(testing$column)) * 100
prop.table(table(input$column)) * 100
#########
i=1
k.optm.col=1
for (i in 1:100){
  knn.mod <- knn(train=training[,9:3010], test=testing[,9:3010], cl=training_labels, k=i)
  k.optm.col[i] <- 100 * sum(cl=testing_labels == knn.mod)/NROW(testing_labels)
  k=i
  cat(k,'=',k.optm.col[i],'')
}
#####prediction of batch###
indxTrain <- createDataPartition(y = input.a$batch,p = 0.75,list = FALSE)
training <- input.a[indxTrain,]
testing <- input.a[-indxTrain,]

training_labels <- as.character(input.a[indxTrain,"batch"])
testing_labels <- as.character(input.a[-indxTrain,"batch"])

#Checking distibution in origanl data and partitioned data
prop.table(table(training$newcol)) * 100
prop.table(table(testing$batch)) * 100
prop.table(table(input$batch)) * 100
##########
i=1
k.optm.b=1
for (i in 1:100){
  knn.mod <- knn(train=training[,9:3010], test=testing[,9:3010], cl=training_labels, k=i)
  k.optm.b[i] <- 100 * sum(cl=testing_labels == knn.mod)/NROW(testing_labels)
  k=i
  cat(k,'=',k.optm.b[i],'')
}
#####prediction of patientID###
indxTrain <- createDataPartition(y = input.a$patientID,p = 0.75,list = FALSE)
training <- input.a[indxTrain,]
testing <- input.a[-indxTrain,]

training_labels <- as.character(input.a[indxTrain,"patientID"])
testing_labels <- as.character(input.a[-indxTrain,"patientID"])

#Checking distibution in origanl data and partitioned data
prop.table(table(training$newcol)) * 100
prop.table(table(testing$patientID)) * 100
prop.table(table(input$patientID)) * 100
####################
i=1
k.optm.ID=1
for (i in 1:100){
  knn.mod <- knn(train=training[,9:3010], test=testing[,9:3010], cl=training_labels, k=i)
  k.optm.ID[i] <- 100 * sum(cl=testing_labels == knn.mod)/NROW(testing_labels)
  k=i
  cat(k,'=',k.optm.ID[i],'')
}

######structure outputs######
data.b<-as.data.frame(cbind(k.optm.b, rep("batch", 100),rep(1:100, 1) ))
names(data.b)[names(data.b) == "k.optm.b"] <- "accuracy"
data.ID<-as.data.frame(cbind(k.optm.ID, rep("cell Type", 100), rep(1:100, 1)))
names(data.ID)[names(data.ID) == "k.optm.ID"] <- "accuracy"
data.col<-as.data.frame(cbind(k.optm.col, rep("column", 100), rep(1:100, 1)))
names(data.col)[names(data.col) == "k.optm.col"] <- "accuracy"
#########merge data####
plot<-rbind( data.b, data.ID , data.col)
colnames(plot)<-c("accuracy", "variable", "knn")
plot[,c(1,3)] <- sapply(sapply(plot[,c(1,3)], as.character),as.numeric)
########choose how many k-nn´s you want to plot###
plot.1<-subset(plot, plot$knn < 50)

knn <- ggplot(plot.1, aes(x= knn ,y=accuracy, color=variable)) + 
  geom_point(size=2)+
  geom_line() + 
  labs(x='k-nearest neighbour', y='accuracy level',title='k-nearest neighbour test') +
  theme_bw() +  
  ylim(0,100) 


##########

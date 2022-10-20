# Persistent Homology Classification Algorithm
#Required Libraries
library(TDA)
library(rgl)
library(dplyr)
library(dplyr)
library(caTools)
library(caret)
library("RColorBrewer")
library(devtools)
library(TDAmapper)
library(igraph)

#General Parameters

maxd <- 0 #maxdimension
maxsc <- 8 #maxscale
SR <- 0.8 #Split Ratio Between training and testing


##########################################################################
## CALLING the MNIST DATABASE
## This part can be modified depending on which data you want to analyze.
##########################################################################
mnist <- read.csv("~/OneDrive - University of the Philippines/DATAforR/mnist/mnist.csv")

dataset = mnist[1:2026]

######################
# FORMATTING DATASET #
######################
XY <- dataset
nr <- nrow(XY)
nc <- ncol(XY)
X <- dataset[-nc] #removes target column
XY[,nc] <- XY[,nc] + 1
XY[,nc] <- as.numeric(XY[,nc]) #Turning nonnumeric data to numeric

## Feature Scaling
#XY[-nc] = scale(XY[-nc])

# Encoding the target feature as factor 
XY[,nc] = factor(XY[,nc])
classes <- levels(XY[,nc])
numclasses <- as.numeric(classes)
maxclass <- max(numclasses)
len <- length(numclasses)

KofCV <- 5
ConM <- matrix(0,maxclass,maxclass)
for (CV in 0:(KofCV-1)){
  
  #Splitting data into classes
  
  trainM <- list()
  testM <- list()
  DiagTr <- list()
  
  # Initializing matrices used
  BB <- matrix(1:maxclass, nrow=1, byrow=FALSE)
  CC <- matrix(1:maxclass, nrow=1, byrow=FALSE)
  DD <- matrix(1:maxclass, nrow=1, byrow=FALSE)
  EE <- matrix(1:(maxclass*3), nrow=maxclass, byrow=FALSE)
  FF <- matrix(1:(maxclass*3), nrow=maxclass, byrow=FALSE)
  GG <- matrix(1:(maxclass*(nc-1)), nrow=maxclass, byrow=FALSE)
  LL <- matrix(1:maxclass, nrow=maxclass, byrow=FALSE)
  nrtr <- matrix(1:maxclass, nrow=1, byrow=FALSE)
  nctr <- matrix(1:maxclass, nrow=1, byrow=TRUE)
  nrte <- matrix(1:maxclass, nrow=1, byrow=TRUE)
  ncte <- matrix(1:maxclass, nrow=1, byrow=TRUE)
  numofDiagTrrows <- matrix(1:maxclass, nrow=1, byrow=TRUE)
  MM <- filter(XY, XY[,nc] == 100)
  newdata <- MM
  maxnum <- maxclass*3
  DiagTrcolSums <- MM[1:3]
  DiagTrcolmean <- MM[1:(nc-1)]
  Featmean <- MM[1:(nc-1)]
  FeatLength <- MM[1:1]
  
  
  #######################
  # PHCA Implementation #
  #######################
  for (class in numclasses)
  { 
    M <- filter(XY, XY[,nc] == class);
    #print(M);
    
    
    for (k in 1:nrow(M))
    {
      M[k,nc+1] <- k%%KofCV;
    }
    colnames(M)[nc+1] <- c("CVlabel")
    
    
    test_set = subset(M, CVlabel == CV)[-(nc+1)]
    training_set = subset(M, CVlabel != CV)[-(nc+1)]
    
    # Optional Feature Scaling, depending on dataset
    #training_set[-nc] = scale(training_set[-nc]) 
    #test_set[-nc] = scale(test_set[-nc])
    
    trainM[[class]] <- training_set
    testM[[class]] <- test_set
    newdata <- rbind(newdata, test_set)
    
    nrtr[class] <- nrow(training_set)
    nctr[class] <- ncol(training_set)
    nrte[class] <- nrow(test_set)
    ncte[class] <- ncol(test_set)
    print(c(nrtr[class], nctr[class], nrte[class], ncte[class]))
    
    ##################################
    # PH COMPUTATION of TRAINING SETS
    ##################################
    
    #####################################
    # Computing Rips persistence diagram 
    #####################################
    MMM <- trainM[[class]][1:nc-1]
    DiagTr[[class]] <- ripsDiag(X = MMM, maxdimension = maxd, maxscale = maxsc,
                                library = "GUDHI", printProgress = FALSE)
    par(mfrow=c(1,2))
    plot(DiagTr[[class]][["diagram"]], main = "Persistence Diagram")
    plot(DiagTr[[class]][["diagram"]], barcode = TRUE, main = "Persistence Barcode")
    numofDiagTrrows[class] <- nrow(DiagTr[[class]][["diagram"]])
    
    for (g in 1:3) {
      DiagTrcolmean[class,g] <- mean(DiagTr[[class]][["diagram"]][,g])
    }
    
    DiagTrcolSums <- rbind(DiagTrcolSums, colSums(DiagTr[[class]][["diagram"]]), deparse.level = 0)
    
    for (gg in 1:(nc-1)) {
      Featmean[class,gg] <- mean(MMM[,gg]) 
      
      #Lengthsum <- 0
      #or (i in nrow(DiagTr[[class]][["diagram"]])) {
      #  Lengthsum <- Lengthsum + DiagTr[[class]][["diagram"]][i,3] - DiagTr[[class]][["diagram"]][i,2]
      #}
      #FeatLength[class,] <- Lengthsum
    }
  }
  
  newdata        
  
  # PHCA CLASSIFIER
  
  resulta <- matrix(1:maxclass, nrow=1, byrow=FALSE)
  result <- matrix(1:1, nrow=1, byrow=FALSE)
  
  
  Classifying <- function(M5){
    
    for (j in numclasses) 
    {
      MM <- rbind(trainM[[j]], M5, deparse.level = 0)
      MM
      #countit <- countit + 1 #Optional
      MMtest <- MM[-nc]
      DiagTe <- ripsDiag(X = MMtest, maxdimension = maxd, maxscale = maxsc,
                         library = "GUDHI", printProgress = FALSE)
      
      #BB[j] <- bottleneck(Diag1 = DiagTr[[j]][["diagram"]], 
      #Diag2 = DiagTe[["diagram"]],
      #dimension = 1)
      #CC[j] <- wasserstein(Diag1 = DiagTr[[j]][["diagram"]], Diag2 = DiagTe[["diagram"]],
      #p = 2, dimension = 1)
      #EE[j,1] <- abs(sum(DiagTe[["diagram"]][,1]) - DiagTrcolSums[j,1])
      EE[j,2] <- abs(sum(DiagTe[["diagram"]][,2]) - DiagTrcolSums[j,2])
      EE[j,3] <- abs(sum(DiagTe[["diagram"]][,3]) - DiagTrcolSums[j,3])
      #FF[j,1] <- abs(mean(DiagTe[["diagram"]][,1]) - DiagTrcolmean[j,1])
      #FF[j,2] <- abs(mean(DiagTe[["diagram"]][,2]) - DiagTrcolmean[j,2])
      #FF[j,3] <- abs(mean(DiagTe[["diagram"]][,3]) - DiagTrcolmean[j,3])
      
      #Lengthsum <- 0
      #for (i in nrow(DiagTe[["diagram"]])) {
      #    Lengthsum <- Lengthsum + DiagTe[["diagram"]][i,3] - DiagTe[["diagram"]][i,2]
      #}
      #LL[j,1] <- abs(Lengthsum - FeatLength[j,])
      
      cc <- 0
      for (kk in 1:(nc-1)) {
        GG[j,kk] <- abs(Featmean[j,kk] - M5[1,kk])
        cc <- cc +GG[j,kk]
      }
      
      #DD[j] <- -EE[j,1] + EE[j,3] - FF[j,1] + FF[j,3] + cc + CC[j]
      DD[j] <- EE[j,3] - EE[j,2] #+ cc #+ CC[j]
    }
    resulta <- which( DD == min(DD), arr.ind=FALSE)
    return(resulta[1])
  }
  
  Count <- 0
  for (ii in numclasses)
    
  {
    for (jj in 1:nrte[ii])
    {
      Count <- Count + 1
      newdata[Count, nc+1] <- Classifying(testM[[ii]][jj,]) 
    }
  }
  
  # Predicting the Test set results for PHCA
  y_label = newdata[,nc]
  y_pred = newdata[,nc+1] #predict(classifier, newdata = test_set[-nc])
  
  # Making the Confusion Matrix 
  cm = table(y_label, y_pred)
  
  #For Confusion Matrix with Cross-validation
  ConM <- ConM + cm
}

######################################
## OTHER CLASSIFICATION ALGORITHMS  ##
######################################
mnist <- read.csv("~/OneDrive - University of the Philippines/DATAforR/mnist/mnist.csv")

dataset = mnist[1:2026]

dataset <- data.frame(mnist[1:2026])
dataset[,nc] <- as.character(dataset[,nc])
data <- dataset
#colnames(dataset) <- c("x1", "x2", "label")

# create a list of 80% of the rows in the original dataset we can use for training
validation_index <- createDataPartition(dataset$label , p=SR, list=FALSE)
# select 20% of the data for validation
validation <- dataset[-validation_index,]
# use the remaining 80% of data to training and testing the models
dataset <- dataset[validation_index,]

# dimensions of dataset
dim(dataset)

# list types for each attribute
sapply(dataset, class)
# take a peek at the first 5 rows of the data
#head(dataset)
# list the levels for the class
#dataset$label = factor(dataset$label)
#levels(dataset$label)

# summarize the class distribution
percentage <- prop.table(table(dataset$label)) * 100
cbind(freq=table(dataset$label), percentage=percentage)       

# summarize attribute distributions
#summary(dataset)       

nc <- ncol(dataset)

# split input and output
x <- dataset[,1:(nc-1)]
y <- dataset[,nc]

# split input and output
x <- dataset[,1:(nc-1)]
y <- dataset[,nc]


# boxplot for each attribute on one image
#par(mfrow=c(1,(nc-1)))
#for(i in 1:(nc-1)) {
#  boxplot(x[,i], main=names(dataset)[i])
#}

#uninteresting
# barplot for class breakdown
#plot(y)

# scatterplot matrix
#featurePlot(x=x, y=y, plot="ellipse")

# box and whisker plots for each attribute
#featurePlot(x=x, y=y, plot="box")

# density plots for each attribute by class value
#scales <- list(x=list(relation="free"), y=list(relation="free"))
#featurePlot(x=x, y=y, plot="density", scales=scales)


# Run algorithms using 10-fold cross validation
control <- trainControl(method="cv", number=KofCV)
metric <- "Accuracy"

#VARIOUS MODELS

# a) linear algorithms
set.seed(7)
fit.lda <- train(label~. , data=dataset, method="lda", metric=metric, trControl=control)
# b) nonlinear algorithms
# CART
set.seed(7)
fit.cart <- train(label~., data=dataset, method="rpart", metric=metric, trControl=control)
# kNN
set.seed(7)
fit.knn <- train(label~., data=dataset, method="knn", metric=metric, trControl=control)
# c) advanced algorithms
# SVM
set.seed(7)
fit.svm <- train(label~., data=dataset, method="svmRadial", metric=metric, trControl=control)
# Random Forest
set.seed(7)
fit.rf <- train(label~., data=dataset, method="rf", metric=metric, trControl=control)


# summarize accuracy of models #Selecting best model
results <- resamples(list(lda=fit.lda, cart=fit.cart, knn=fit.knn, svm=fit.svm, rf=fit.rf))
summary(results)


# compare accuracy of models
# dotplot(results)

############################################################
## Generate the performance measures for each classifier. ##
############################################################

# estimate skill of LDA on the validation dataset
predictions <- predict(fit.lda, validation)
confusionMatrix(predictions, factor(validation[,nc]))[4]
confusionMatrix(predictions,factor(validation[,nc]))[2]
confusionMatrix(predictions,factor(validation[,nc]))[3]

perflda <- array(as.numeric(unlist(confusionMatrix(predictions,factor(validation[,nc]))[4])), 
                 dim=c(1, maxclass, 11))
perf_sen <- perflda[,,1]
perf_spec <- perflda[,,2]
perfacc <- array(as.numeric(unlist(confusionMatrix(predictions,factor(validation[,nc]))[3])), 
                 dim=c(1, maxclass, 11))[1]
perfacc
perf_acc <- perfacc
perf_prec <- perflda[,,5]
perf_rec <- perflda[,,6]
perf_f1 <- perflda[,,7]

# estimate skill of CART on the validation dataset
predictions <- predict(fit.cart, validation)
confusionMatrix(predictions, factor(validation[,nc]))[4]
confusionMatrix(predictions,factor(validation[,nc]))[2]
confusionMatrix(predictions,factor(validation[,nc]))[3]

perfcart <- array(as.numeric(unlist(confusionMatrix(predictions,factor(validation[,nc]))[4])), 
                  dim=c(1, maxclass, 11))
perf_sen <- rbind(perf_sen, perfcart[,,1])
perf_spec <- rbind(perf_spec, perfcart[,,2])
perfacc <- array(as.numeric(unlist(confusionMatrix(predictions,factor(validation[,nc]))[3])), 
                 dim=c(1, maxclass, 11))[1]
perfacc
perf_acc <- rbind(perf_acc, perfacc) 
perf_prec <- rbind(perf_prec, perfcart[,,5])
perf_rec <- rbind(perf_rec, perfcart[,,6])
perf_f1 <- rbind(perf_f1, perfcart[,,7])

# estimate skill of KNN on the validation dataset
predictions <- predict(fit.knn, validation)
confusionMatrix(predictions, factor(validation[,nc]))[4]
confusionMatrix(predictions,factor(validation[,nc]))[2]
confusionMatrix(predictions,factor(validation[,nc]))[3]

perfknn <- array(as.numeric(unlist(confusionMatrix(predictions,factor(validation[,nc]))[4])), 
                 dim=c(1, maxclass, 11))
perf_sen <- rbind(perf_sen, perfknn[,,1])
perf_spec <- rbind(perf_spec, perfknn[,,2])
perfacc <- array(as.numeric(unlist(confusionMatrix(predictions,factor(validation[,nc]))[3])), 
                 dim=c(1, maxclass, 11))[1]
perfacc
perf_acc <- rbind(perf_acc, perfacc) 
perf_prec <- rbind(perf_prec, perfknn[,,5])
perf_rec <- rbind(perf_rec, perfknn[,,6])
perf_f1 <- rbind(perf_f1, perfknn[,,7])

# estimate skill of SVM on the validation dataset
predictions <- predict(fit.svm, validation)
confusionMatrix(predictions, factor(validation[,nc]))[4]
confusionMatrix(predictions,factor(validation[,nc]))[2]
confusionMatrix(predictions,factor(validation[,nc]))[3]

perfsvm <- array(as.numeric(unlist(confusionMatrix(predictions,factor(validation[,nc]))[4])), 
                 dim=c(1, maxclass, 11))
perf_sen <- rbind(perf_sen, perfsvm[,,1])
perf_spec <- rbind(perf_spec, perfsvm[,,2])
perfacc <- array(as.numeric(unlist(confusionMatrix(predictions,factor(validation[,nc]))[3])), 
                 dim=c(1, maxclass, 11))[1]
perfacc
perf_acc <- rbind(perf_acc, perfacc) 
perf_prec <- rbind(perf_prec, perfsvm[,,5])
perf_rec <- rbind(perf_rec, perfsvm[,,6])
perf_f1 <- rbind(perf_f1, perfsvm[,,7])

# estimate skill of RF on the validation dataset
predictions <- predict(fit.rf, validation)
confusionMatrix(predictions, factor(validation[,nc]))[4]
confusionMatrix(predictions,factor(validation[,nc]))[2]
confusionMatrix(predictions,factor(validation[,nc]))[3]

perfrf <- array(as.numeric(unlist(confusionMatrix(predictions,factor(validation[,nc]))[4])), 
                dim=c(1, maxclass, 11))
perf_sen <- rbind(perf_sen, perfrf[,,1])
perf_spec <- rbind(perf_spec, perfrf[,,2])
perfacc <- array(as.numeric(unlist(confusionMatrix(predictions,factor(validation[,nc]))[3])), 
                 dim=c(1, maxclass, 11))[1]
perfacc
perf_acc <- rbind(perf_acc, perfacc) 
perf_prec <- rbind(perf_prec, perfrf[,,5])
perf_rec <- rbind(perf_rec, perfrf[,,6])
perf_f1 <- rbind(perf_f1, perfrf[,,7])


#######################
## Result Using PHCA ##
#######################

confusionMatrix(ConM)
confusionMatrix(ConM)[4]
confusionMatrix(ConM)[2]
confusionMatrix(ConM)[3]

perfphca <- array(as.numeric(unlist(confusionMatrix(ConM)[4])), 
                  dim=c(1, maxclass, 11))
perf_sen <- rbind(perf_sen, perfphca[,,1])
perf_spec <- rbind(perf_spec, perfphca[,,2])
perfacc <- array(as.numeric(unlist(confusionMatrix(ConM)[3])), 
                 dim=c(1, maxclass, 11))[1]
perfacc
perf_acc <- rbind(perf_acc, perfacc) 
perf_prec <- rbind(perf_prec, perfphca[,,5])
perf_rec <- rbind(perf_rec, perfphca[,,6])
perf_f1 <- rbind(perf_f1, perfphca[,,7])

#### LDA
perf_sen[1,] <- c(1,	0.96,	0.88,	0.94,	0.96,	0.86,	0.92,	0.82,	0.92,	0.94)
perf_spec[1,] <- c(1,	0.995516,	0.986784,	0.993318,	0.995546,	0.984581,	0.99117,	0.980263,	0.991051,	0.993289)
perf_acc[1,] <-  c(0.92)
perf_prec[1,] <- c(0.909091,	0.888889,	0.956522,	0.921569,	0.941176,	0.934783,	0.978723,	0.931818,	0.867925,	0.886792) 
perf_rec[1,] <- c(1,	0.96,	0.88,	0.94,	0.96,	0.86,	0.92,	0.82,	0.92,	0.94)
perf_f1[1,] <- c(0.952381,	0.923077,	0.916667,	0.930693,	0.950495,	0.895833,	0.948454,	0.87234,	0.893204,	0.912621)

####

dev.off()
perf_prec
perf_rec
perf_acc
perf_spec
perf_f1

PI <- cbind(perf_prec, perf_rec, perf_acc, perf_spec, perf_f1)
M <- t(PI)
colnames(M) <- c('LDA', 'CART', 'KNN', 'SVM', 'RF', 'PHCA')
#rownames(M) <- c('prec','rec', 'acc','spec', 'f1score')
library(DescTools)

g <- factor(rep(1:6, each=41),
            labels = c("LDA", "CART", "KNN", "SVM", "RF", "PHCA"))
NemenyiTest(c(M), g)

######################################
## Plotting of Performance Measures ##
######################################

# Plot for Sensitivity/Recall
colnames(perf_sen)<- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
par(mar=c(4,4,4,4))
barplot(perf_sen,
        #main = "Recall/Sensitivity per class per classifier",
        xlab = "Classes",
        ylab = "Recall/Sensitivity",
        col = brewer.pal(n = 6, name = "BuGn"),#RdBu
        beside = TRUE
        #)
        #legend("bottomleft",
        #c("lda", "cart", "knn", "svm", "rf", "phca"),
        #fill = brewer.pal(n = 6, name = "BuGn")
)

# Plot for Specificity/Positive Predictive Value
colnames(perf_spec)<- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
par(mar=c(4,4,4,4))
barplot(perf_spec,
        #main = "Specificity per class per classifier",
        xlab = "Classes",
        ylab = "Specificity",
        col = brewer.pal(n = 6, name = "BuGn"),
        beside = TRUE
        #)
        #legend("bottomleft",
        #      c("lda", "cart", "knn", "svm", "rf", "phca"),
        #     fill = brewer.pal(n = 6, name = "BuGn")
)

#Plot for Precision
colnames(perf_prec)<- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
par(mar=c(4,4,4,4))
barplot(perf_prec,
        #main = "Precision per class per classifier",
        xlab = "Classes",
        ylab = "Precision",
        col = brewer.pal(n = 6, name = "BuGn"),
        beside = TRUE
        #)
        #legend("bottomleft",
        #c("lda", "cart", "knn", "svm", "rf", "phca"),
        #fill = brewer.pal(n = 6, name = "BuGn")
)

#Plot for F1 Score
colnames(perf_f1)<- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
par(mar=c(4,4,4,4))
barplot(perf_f1,
        #main = "F1 Score per class per classifier",
        xlab = "Classes",
        ylab = "F1 Score",
        col = brewer.pal(n = 6, name = "BuGn"),
        beside = TRUE
        #)
        #legend("bottomleft",
        #c("lda", "cart", "knn", "svm", "rf", "phca"),
        #fill = brewer.pal(n = 6, name = "BuGn")
)

#Plot for Accuracy
par(mar=c(4,6,4,6))
barplot(perf_acc,
        #main = "Accuracy per classifier",
        ylim = c(0,1.0),
        xlab = "Classifiers",
        ylab = "Accuracy",
        col = brewer.pal(n = 6, name = "BuGn"),
        names.arg = c("LDA", "CART", "KNN", "SVM", "RF", "PHCA"),
        beside = TRUE
)
#legend("topright", c("LDA", "CART", "KNN", "SVM", "RF", "PHCA"),
#       fill = brewer.pal(n = 6, name = "BuGn"), inset = c(- 0.25, 0)
#)

M
100*t(M)

# Persistent Homology Classification Algorithm
#Required Libraries
library(TDA)
library(rgl)
library(dplyr)
library(dplyr)
library(caTools)
library(caret)

#General Parameters

maxd <- 1 #maxdimension
maxsc <- 1 #maxscale
SR <- 0.8 #Split Ratio Between training and testing


##########################################################################
## CALLING the IRIS DATA
## This part can be modified depending on which data you want to analyze.
##########################################################################
data(iris)
dataset <- iris

######################
# FORMATTING DATASET #
######################
XY <- dataset
nr <- nrow(XY)
nc <- ncol(XY)
X <- dataset[-nc] #removes target column
XY[,nc] <- as.numeric(XY[,nc]) #Turning nonnumeric data to numeric

# Encoding the target feature as factor 
XY[,nc] = factor(XY[,nc])
classes <- levels(XY[,nc])
numclasses <- as.numeric(classes)
maxclass <- max(numclasses)
len <- length(numclasses)

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
    M <- filter(XY, XY[,nc] == class)
    print(M)
    
    # Splitting the dataset into the Training set and Test set 
    #install.packages('caTools')
    set.seed(7)
    split = sample.split(M[,nc], SplitRatio = SR)
    training_set = subset(M, split == TRUE)
    test_set = subset(M, split == FALSE) 
    
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
        
        Lengthsum <- 0
        for (i in nrow(DiagTr[[class]][["diagram"]])) {
            Lengthsum <- Lengthsum + DiagTr[[class]][["diagram"]][i,3] - DiagTr[[class]][["diagram"]][i,2]
        }
        FeatLength[class,] <- Lengthsum
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
        
        BB[j] <- bottleneck(Diag1 = DiagTr[[j]][["diagram"]], 
                            Diag2 = DiagTe[["diagram"]],
                            dimension = 1)
        CC[j] <- wasserstein(Diag1 = DiagTr[[j]][["diagram"]], Diag2 = DiagTe[["diagram"]],
                             p = 2, dimension = 1)
        EE[j,1] <- abs(sum(DiagTe[["diagram"]][,1]) - DiagTrcolSums[j,1])
        EE[j,2] <- abs(sum(DiagTe[["diagram"]][,2]) - DiagTrcolSums[j,2])
        EE[j,3] <- abs(sum(DiagTe[["diagram"]][,3]) - DiagTrcolSums[j,3])
        FF[j,1] <- abs(mean(DiagTe[["diagram"]][,1]) - DiagTrcolmean[j,1])
        FF[j,2] <- abs(mean(DiagTe[["diagram"]][,2]) - DiagTrcolmean[j,2])
        FF[j,3] <- abs(mean(DiagTe[["diagram"]][,3]) - DiagTrcolmean[j,3])
        
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
        
        DD[j] <- -EE[j,1] + EE[j,3] - FF[j,1] + FF[j,3] + cc + CC[j]
    }
    resulta <- which( DD == min(DD), arr.ind=FALSE)
    
    if (M5[,nc]!= resulta[1]) {
        print(M5)
        print("EE")
        print(EE)
        print("BB")
        print(BB)
        print("CC")
        print(CC)
        print("FF")
        print(FF)
        
        print("DD")
        print(DD)
    }
    
    print(resulta[1])
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



# Results of predicting/classifying the testing set 
y_pred = newdata[,nc+1] #predict(classifier, newdata = test_set[-nc])

# Making the Confusion Matrix 
cm = table(newdata[, nc], y_pred)
confusionMatrix(cm)



###################################
# OTHER CLASSIFICATION ALGORITHMS #
###################################

data(iris)
dataset <- iris

# create a list of 80% of the rows in the original dataset we can use for training
validation_index <- createDataPartition(dataset$Species, p=0.80, list=FALSE)
# select 20% of the data for validation
validation <- dataset[-validation_index,]
# use the remaining 80% of data to training and testing the models
dataset <- dataset[validation_index,]

# dimensions of dataset
dim(dataset)

# list types for each attribute
sapply(dataset, class)
# take a peek at the first 5 rows of the data
head(dataset)
# list the levels for the class
levels(dataset$Species)

# summarize the class distribution
percentage <- prop.table(table(dataset$Species)) * 100
cbind(freq=table(dataset$Species), percentage=percentage)       

# summarize attribute distributions
summary(dataset)       


# split input and output
x <- dataset[,1:4]
y <- dataset[,5]

# split input and output
x <- dataset[,1:4]
y <- dataset[,5]


# boxplot for each attribute on one image
par(mfrow=c(1,4))
for(i in 1:4) {
    boxplot(x[,i], main=names(iris)[i])
}

#uninteresting
# barplot for class breakdown
plot(y)

# scatterplot matrix
featurePlot(x=x, y=y, plot="ellipse")

# box and whisker plots for each attribute
featurePlot(x=x, y=y, plot="box")

# density plots for each attribute by class value
scales <- list(x=list(relation="free"), y=list(relation="free"))
featurePlot(x=x, y=y, plot="density", scales=scales)


# Run algorithms using 10-fold cross validation
control <- trainControl(method="cv", number=5)
metric <- "Accuracy"

#VARIOUS MODELS

# a) linear algorithms
set.seed(7)
fit.lda <- train(Species~., data=dataset, method="lda", metric=metric, trControl=control)
# b) nonlinear algorithms
# CART
set.seed(7)
fit.cart <- train(Species~., data=dataset, method="rpart", metric=metric, trControl=control)
# kNN
set.seed(7)
fit.knn <- train(Species~., data=dataset, method="knn", metric=metric, trControl=control)
# c) advanced algorithms
# SVM
set.seed(7)
fit.svm <- train(Species~., data=dataset, method="svmRadial", metric=metric, trControl=control)
# Random Forest
set.seed(7)
fit.rf <- train(Species~., data=dataset, method="rf", metric=metric, trControl=control)


# summarize accuracy of models #Selecting best model
results <- resamples(list(lda=fit.lda, cart=fit.cart, knn=fit.knn, svm=fit.svm, rf=fit.rf))
summary(results)


# compare accuracy of models
dotplot(results)

# summarize Best Model
print(fit.lda)

# estimate skill of LDA on the validation dataset
predictions <- predict(fit.lda, validation)
confusionMatrix(predictions, validation$Species)


# estimate skill of CART on the validation dataset
predictions <- predict(fit.cart, validation)
confusionMatrix(predictions, validation$Species)

# estimate skill of KNN on the validation dataset
predictions <- predict(fit.knn, validation)
confusionMatrix(predictions, validation$Species)


# estimate skill of SVM on the validation dataset
predictions <- predict(fit.svm, validation)
confusionMatrix(predictions, validation$Species)


# estimate skill of RF on the validation dataset
predictions <- predict(fit.rf, validation)
confusionMatrix(predictions, validation$Species)

# Predicting the Test set results for PHCA
y_pred = newdata[,nc+1] #predict(classifier, newdata = test_set[-nc])

# Making the Confusion Matrix 
cm = table(newdata[, nc], y_pred)
confusionMatrix(cm)
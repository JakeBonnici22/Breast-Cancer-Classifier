install.packages("mlbench")
install.packages("mice")
install.packages("bestglm")
install.packages("glmnet")
install.packages("nclSLR", repos="http://R-Forge.R-project.org")
install.packages("MASS")
install.packages("leaps")
library(mlbench)
library(mice)
library(bestglm)
library(glmnet)
library(nclSLR)
library(MASS)
library(leaps)

# Understanding the data
data(BreastCancer)
dim(BreastCancer)
head(BreastCancer)
summary(BreastCancer)
table(is.na(BreastCancer))
#C  leaning the data
class(BreastCancer)
str(BreastCancer)


# Changing from factors to quantitative(numeric)
cols.num <- c("Cl.thickness","Cell.size", "Cell.shape", "Marg.adhesion", "Epith.c.size", "Bare.nuclei", "Bl.cromatin", "Normal.nucleoli", "Mitoses")
BreastCancer[cols.num] <- sapply(BreastCancer[cols.num],as.numeric)
sapply(BreastCancer, class)
head(BreastCancer)
BreastCancer
# Replacing all rows with missing values
?mice
dataset_impute <- mice(BreastCancer[,2:10], print=FALSE)
head(dataset_impute)
BreastCancer <- cbind(BreastCancer[,11, drop = FALSE], mice::complete(dataset_impute, 1))
head(BreastCancer)
str(BreastCancer)

str(BreastCancer[,-1])
?mice::complete
pairs(BreastCancer)
?pairs
pairs(BreastCancer[,2:10], col=BreastCancer[,1], cex.labels=1.7)
cor(BreastCancer[,-1])
heatmap(data.matrix(BreastCancer), Rowv=NA, Colv=NA, labRow = NA)
colMeans(BreastCancer[,-1])
var(BreastCancer[,-1])
num_breastcancer = BreastCancer[,-1]
apply(num_breastcancer, 2, var)

  
# plot(BreastCancer[,4], BreastCancer[,3])
# plot(BreastCancer[,2], BreastCancer[,5])
# plot(BreastCancer[,8], BreastCancer[,10])

s = var(BreastCancer[,-1])
s_sq = diag(s)
total_variation = sum(s_sq)




#############################
#############################
  ## LOGISTIC REGRESSION ##
#############################
#############################
  
#Predictor variables
head(BreastCancer)
x1orig = BreastCancer[,2:10]
?scale
x1 = scale(x1orig)
# Response Variable
y = as.integer(BreastCancer$Class)-1
##Combine to create a new data frame.
Breast_canc = data.frame(x1, y)
head(Breast_canc)
Breast_canc



## Store n and p
n = nrow(Breast_canc); p = ncol(Breast_canc) - 1
logreg_fit = glm(y~.,data=Breast_canc, family="binomial")
summary(logreg_fit)
bss_fit_AIC = bestglm(Breast_canc, family=binomial, IC="AIC")
bss_fit_BIC  = bestglm(Breast_canc, family=binomial, IC="BIC")

bss_fit_AIC$Subsets
bss_fit_BIC$Subsets
(best_AIC = bss_fit_AIC$ModelReport$Bestk)
(best_BIC = bss_fit_BIC$ModelReport$Bestk)


## Create multi-panel plotting device
par(mfrow=c(1,2))
## Produce plots, highlighting optimal value of k
plot(0:p, bss_fit_AIC$Subsets$AIC, main="Best subset selection fit using AIC.",
     xlab="Number of predictors", ylab="AIC", type="b", lwd = 2)
points(best_AIC, bss_fit_AIC$Subsets$AIC[best_AIC+1], col="red", pch=16, cex=2)

plot(0:p, bss_fit_BIC$Subsets$BIC, main="Best subset selection fit using BIC.",
     xlab="Number of predictors", ylab="BIC", type="b", lwd = 2)
points(best_BIC, bss_fit_BIC$Subsets$BIC[best_BIC+1], col="red", pch=16, cex=2)

pstar = 7
## AIC 
## Check which predictors are in the 6-predictor model
bss_fit_AIC$Subsets[pstar+1,]
## Construct a reduced data set containing only the selected predictor
(indices = as.logical(bss_fit_AIC$Subsets[pstar+1, 2:(p+1)]))

Breast_canc_red = data.frame(x1[,indices], y)
## Obtain regression coefficients for this model
logreg1_fit = glm(y ~ ., data=Breast_canc_red, family="binomial")
summary(logreg1_fit)

#BIC
pstar = 5
## Check which predictors are in the 5-predictor model
bss_fit_BIC$Subsets[pstar+1,]
## Construct a reduced data set containing only the selected predictor
(indices = as.logical(bss_fit_BIC$Subsets[pstar+1, 2:(p+1)]))

Breast_canc_red = data.frame(x1[,indices], y)
## Obtain regression coefficients for this model
logreg1_fit = glm(y ~ ., data=Breast_canc_red, family="binomial")
logreg1_fit
summary(logreg1_fit)





#############################
#############################
      ## LASSO/RIDGE ##
#############################
#############################

# par(mfrow=c(1,1))
## Choose grid of values for the tuning parameter
grid = 10^seq(-4, -1, length.out=100)
grid
## Fit a model with LASSO/RIDGE penalty for each value of the tuning parameter
#LASSO
lasso_fit = glmnet(x1, y, family="binomial", alpha=1, standardize=FALSE, lambda=grid)
plot(lasso_fit, xvar="lambda", col=rainbow(p), label=TRUE)
title("Lasso", line = 2.5)

lasso_cv_fit = cv.glmnet(x1, y, family="binomial", alpha=1, standardize=FALSE, lambda=grid,
                         type.measure="class")
plot(lasso_cv_fit)
title("Cross Validation for the Lasso", line = 2.5)




#RIDGE
ridge_fit = glmnet(x1, y, family="binomial", alpha=0, standardize=FALSE, lambda=grid)
plot(ridge_fit, xvar="lambda", col=rainbow(p), label=TRUE)
title("Ridge", line = 2.5)
ridge_cv_fit = cv.glmnet(x1, y, family="binomial", alpha=0, standardize=FALSE, lambda=grid,
                         type.measure="class")
ridge_cv_fit
plot(ridge_cv_fit)
title("Cross Validation for the Ridge", line = 2.5)




## Identify the optimal value for the tuning parameter
#LASSO
(lambda_lasso_min = lasso_cv_fit$lambda.min)
which_lambda_lasso = which(lasso_cv_fit$lambda == lambda_lasso_min)
#RIDGE
(lambda_ridge_min = ridge_cv_fit$lambda.min)
which_lambda_ridge = which(ridge_cv_fit$lambda == lambda_ridge_min)
## Find the parameter estimates associated with optimal value of the tuning parameter
coef(lasso_fit, s=lambda_lasso_min)
coef(ridge_fit, s=lambda_ridge_min)


## LASSO TEST ERROR
phat = predict(lasso_fit, x1, s=lambda_lasso_min, type="response")
## Compute fitted (i.e. predicted) values:
yhat = ifelse(phat > 0.5, 1, 0)
## Calculate confusion matrix:
(confusion = table(Observed=y, Predicted=yhat))
1 - mean(y==yhat)


## RIDGE TEST ERROR
phat = predict(ridge_fit, x1, s=lambda_ridge_min, type="response")
## Compute fitted (i.e. predicted) values:
yhat = ifelse(phat > 0.5, 1, 0)
## Calculate confusion matrix:
(confusion = table(Observed=y, Predicted=yhat))
1 - mean(y==yhat)

##############################
        ##TEST ERROR##
##############################
## Fit logistic regression model to training data:
logreg1_train = glm(y ~ ., data=Breast_canc_red[1:500,], family="binomial")
logreg1_train

## Compute fitted values for the validation data:
phat_test = predict(logreg1_train, Breast_canc_red[501:699,], type="response")#500
yhat_test = ifelse(phat_test > 0.5, 1, 0)
yhat_test#500
## Compute test error
1 - mean(y[501:699] == yhat_test)




################################
##LINEAR DISCRIMINANT ANALYSIS##
################################

## Perform LDA on the training data - note that we need to convert the vector of predictors
## into a matrix because the linDA function expects its variables argument to be a matrix
## or data frame

##############
##Cell.shape##
##############
linDA(variables=matrix(Breast_canc$Cell.shape[1:500], ncol=1),
      group=Breast_canc$y[1:500])


(lda_fit = linDA(variables = x1, group=y))
lda_fit$error_rate
lda_fit$functions

x2 = seq(-6, 6, 0.001)
Q1 = -0.947 - 1.551*x2
Q2 = -2.493 + 2.902*x2
plot(x2, Q1, type="l", col="red", xlab="x2", ylab="Discriminant function")
lines(x2, Q2, col="blue")
abline(v=0.3472, col="grey", lty="dashed")
title("LDA for Cell Shape")


## Fit the LDA classifier using the training data
lda_train = lda(y~Cell.shape, data=Breast_canc[1:500,])
## Compute fitted values for the validation data
lda_test = predict(lda_train, Breast_canc[501:699,])
yhat_test = lda_test$class
yhat_test
## Calculate (test) confusion matrix
(confusion = table(Observed=Breast_canc$y[501:699], Predicted=yhat_test))

1 - mean(Breast_canc$y[Breast_canc$y[501:699]] == yhat_test)




################
##Cl.thickness##
################
linDA(variables=matrix(Breast_canc$Cl.thickness[1:500], ncol=1),
      group=Breast_canc$y[1:500])

(lda_fit = linDA(variables = x1, group=y))
lda_fit$error_rate
lda_fit$functions


x2 = seq(-6, 6, 0.001)
Q1 = -0.802 - 1.079*x2
Q2 = -1.990 + 2.2023*x2
plot(x2, Q1, type="l", col="red", xlab="x2", ylab="Discriminant function")
lines(x2, Q2, col="blue")
abline(v=0.362, col="grey", lty="dashed")
title("LDA for Cell Thickness")


## Fit the LDA classifier using the training data
lda_train = lda(y~Cl.thickness, data=Breast_canc[1:500,])
## Compute fitted values for the validation data
lda_test = predict(lda_train, Breast_canc[501:699,])
yhat_test = lda_test$class
yhat_test
## Calculate (test) confusion matrix
(confusion = table(Observed=Breast_canc$y[501:699], Predicted=yhat_test))

1 - mean(Breast_canc$y[Breast_canc$y[501:699]] == yhat_test)




####################################
## Quadratic discriminant analysis##
####################################

## Perform QDA on the training data - note that we need to convert the vector of predictors
## into a matrix because the linDA function expects its variables argument to be a matrix
## or data frame
##############
##Cell.shape##
##############

quaDA(variables=matrix(Breast_canc$Cell.shape[1:500], ncol=1),
      group=Breast_canc$y[1:500], functions=TRUE)

x2 = seq(-10, 10, 0.001)
Q1 = -0.7212 - 4.2403 * x2 - 3.6867 * x2^2
Q2 = -1.5661 + 1.4676 * x2 - 0.6818 * x2^2
plot(x2, Q1, type="l", col="red", xlab="x2", ylab="Discriminant function")
lines(x2, Q2, col="blue")
title("QDA for Cell Shape")


## Fit the QDA classifier using the training data
qda_train = qda(y~Cell.shape, data=Breast_canc[1:500,])
## Compute fitted values for the validation data
qda_test = predict(qda_train, Breast_canc[501:699,])
yhat_test = qda_test$class
## Calculate (test) confusion matrix
(confusion = table(Observed=Breast_canc$y[501:699], Predicted=yhat_test))

## Calculate the test error
1 - mean(Breast_canc$y[Breast_canc$y[501:699]] == yhat_test)



################
##Cl.thickness##
################

quaDA(variables=matrix(Breast_canc$Cl.thickness[1:500], ncol=1),
      group=Breast_canc$y[1:500], functions=TRUE)

x2 = seq(-10, 10, 0.001)
Q1 = -0.4250 - 1.5169 * x2 - 1.3578 * x2^2
Q2 = -1.5190 + 1.4003 * x2 - 0.6688 * x2^2
plot(x2, Q1, type="l", col="red", xlab="x2", ylab="Discriminant function")
lines(x2, Q2, col="blue")
title("QDA for Cell Thickness")



## Fit the QDA classifier using the training data
qda_train = qda(y~Cl.thickness, data=Breast_canc[1:500,])
## Compute fitted values for the validation data
qda_test = predict(qda_train, Breast_canc[501:699,])
yhat_test = qda_test$class
## Calculate (test) confusion matrix
(confusion = table(Observed=Breast_canc$y[501:699], Predicted=yhat_test))

## Calculate the test error
1 - mean(Breast_canc$y[Breast_canc$y[501:699]] == yhat_test)


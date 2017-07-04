rm(list=ls())

#' Centers data dividing by the euclidean norm as in Fu or by sd
#' @param data A matrix to center and normalize columwise
#' @return The normalized matrix
prepareDataFu <- function(data, func = sd){
  # centering as sugested in http://gastonsanchez.com/visually-enforced/how-to/2014/01/15/Center-data-in-R/
  centered <- data - rep(colMeans(data), rep.int(nrow(data),ncol(data)))
  # divide by given function
  prepared <- centered/apply(centered, 2, func)
  return(prepared)
}

#' Normalizes data columwise (R-version)
#' @param data A matrix to standarize
#' @return The standarized matrix
prepareData <- function(data){
  return(scale(data))
}

#' Computes the euclidean norm
#' @param x a number or vector
#' @return the euclidean norm
#' @examples
#' euclidnorm(3)
#' euclidnorm(3:4)
#' euclidnorm(c(3,4))
euclidnorm <- function(x){
  return(sqrt(sum(x^2)))
}


norm1 <- function(x){
  return(sum(abs(x)))
}

normInf <- function(x){
  return(max(abs(x)))
}

#' Shooting function at 0
#' @param j column index
#' @param beta weights
#' @param X nxp matrix of data
#' @param y n-vector of responses
shoot <- function(j, beta, X, y) {
  tmpSum <- 0
  for(i in 1:ncol(X)){
    if(i!=j){
      tmpSum <- tmpSum + sum(2*t(X[,j])*X[,i]*beta[i])
    }
  }
  #LHS <- sum(2*t(X[,j])*X[,j]*beta[j])+tmpSum-2*sum(t(X[,j]*y))
  LHS <- tmpSum-2*sum(t(X[,j]*y))
  return(LHS)
}

#' Shooting algorithm as in Fu
#' @param X design matrix
#' @param y vector of responses
#' @param lambda tuning parameter
#' @param eps tolerance
#' @return The vector of parameters, steps performed and indicator of convergence
shooting_fu <- function(X, y, lambda, eps = 1e-6, max_steps = 1000){

  p <- ncol(X)

  # Start with $\beta_OLS$
  myFormula <- as.formula(paste("y~",paste(colnames(X),collapse ="+")))
  OLS <- lm(myFormula, data=as.data.frame(X))
  beta_hat <- tail(OLS$coefficients,-1)

  converged <- FALSE
  step <- 0

  while (!converged & (step < max_steps)){

    beta_old <- beta_hat

    # Calculate
    for(j in 1:p){
      # shoot at 0
      S_0 <- shoot(j, beta_hat, X, y)
      if(S_0 > lambda){
        beta_hat[j] <- (lambda - S_0)/(2*euclidnorm(X[,j])^2)
      } else if (S_0 < -lambda) {
        beta_hat[j] <- (-lambda - S_0)/(2*euclidnorm(X[,j])^2)
      } else {
        beta_hat[j] <- 0
      }
    }

    step <- step + 1
    converged <- euclidnorm(beta_hat-beta_old) < eps
  }
  output <- list(beta = beta_hat, step = step, converged = converged)
  return(output)
}


shooting_fu2 <- function(X, y, lambda, eps = 1e-6, max_steps = 1000){

  p <- ncol(X)
  two_XX = 2*t(X)%*%X # Equivalent: 2*crossprod(X)
  two_Xy = 2*crossprod(X,y)

  # Start with $\beta_OLS$
  myFormula <- as.formula(paste("y~",paste(colnames(X),collapse ="+")))
  OLS <- lm(myFormula, data=as.data.frame(X))
  beta_hat <- tail(OLS$coefficients,-1)

  converged <- FALSE
  step <- 0

  while (!converged & (step < max_steps)){

    beta_old <- beta_hat

    # Calculate
    for (j in 1:p){
      # shoot at 0
      S_0 <- sum(two_XX[j,] %*% beta_hat) - two_Xy[j] - beta_hat[j] * two_XX[j,j]
      if (S_0 > lambda){
        beta_hat[j] <- (lambda - S_0)/two_XX[j,j]
      } else if (S_0 < -lambda) {
        beta_hat[j] <- (-lambda - S_0)/two_XX[j,j]
      } else {
        beta_hat[j] <- 0
      }
    }

    step <- step + 1
    converged <- norm1(beta_hat-beta_old) < eps
  }
  output <- list(beta = beta_hat, step = step, converged = converged)
  return(output)
}

#' Shooting SCD algorithm (Shalev-Shwartz and Tewari)
#' @param X design matrix
#' @param y vector of responses
#' @param lambda tuning parameter
#' @param eps tolerance
#' @return The vector of parameters, steps performed and indicator of convergence
shooting_scd <- function(X, y, lambda, eps = 1e-6, max_steps = 1000){

  p <- ncol(X)
  m <- nrow(X)

  converged <- FALSE
  step <- 0
  beta_hat <- rep(0,p)

  while (!converged & (step < max_steps)){

    beta_old <- beta_hat

    # Sample j uniformly at random from {1...p}
    j <- sample(1:p,1)
    # calculate j-derivative
    tmp <- 0
    for (i in 1:m){
      tmp <- tmp + (sum(beta_hat*X[i,])-y[i])*X[i,j]
      #tmp <- tmp + (crossprod(beta_hat,X[i,])-y[i])*X[i,j]
    }
    g_j <- tmp/m

    # Thresholding
    w_j_g_j <- beta_hat[j]-g_j

    if (w_j_g_j > lambda){
      beta_hat[j] <- (w_j_g_j-lambda)
    } else if (w_j_g_j < -lambda) {
      beta_hat[j] <- (w_j_g_j+lambda)
    } else {
      beta_hat[j] <- 0
    }
    step <- step + 1
    converged <- euclidnorm(beta_hat-beta_old) < eps
  }
  output <- list(beta = beta_hat, step = step, converged = converged)

  return(output)
}
###############################################################

# load data (Prostate data set)
#library(ElemStatLearn) # Prostate data set
library("lasso2", lib.loc="~/R/win-library/3.4")
library("lasso2", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
data(Prostate)

# skip last column
#prostate=prostate[,-ncol(prostate)]
prostate <- Prostate
# Normalized data
prostateNormalized <- as.data.frame(prepareData(prostate[,1:8]))
# Normalized data Fu
prostateNormalizedFuEuclid <- prepareDataFu(prostate[,1:8], euclidnorm)
prostateNormalizedFuSd <- prepareDataFu(prostate[,1:8], sd)
prostateNormalizedFuAbs <- prepareDataFu(prostate[,1:8], abs)

lpsa <- prostate[,9]

##### OLS calculation
# Raw data
beta_OLS_raw <- lm(lpsa~., data=prostate)
# Normalized data
beta_OLS_normalized <- lm(lpsa~., data=cbind(prostateNormalized,lpsa))
# Normalized Fu
beta_OLS_normalizedFuSd <- lm(lpsa~., data=cbind(prostateNormalizedFuSd,lpsa))
beta_OLS_normalizedFuEuclid <- lm(lpsa~., data=cbind(prostateNormalizedFuEuclid,lpsa))
beta_OLS_normalizedFuAbs <- lm(lpsa~., data=cbind(prostateNormalizedFuAbs,lpsa))

OLS_table <- round(matrix(c(coef(beta_OLS_raw), coef(beta_OLS_normalized),
                      coef(beta_OLS_normalizedFuSd), coef(beta_OLS_normalizedFuEuclid),
                      coef(beta_OLS_normalizedFuAbs)),
                    ncol = 5),digits=3)
row.names(OLS_table) <- names(coef(beta_OLS_raw))
colnames(OLS_table) <- c("raw","normalized", "normalizedFuSd", "normalizedFuEuclid",
                         "normalizedFuAbs")
OLS_table

# Table 2. Correlation Matrix. OK with FU
round(cor(prostateNormalized[,1:8]),digits=3)
# Condition number (should be 16.9?)
condNumber <- kappa(round(cor(prostateNormalized[,1:8]),digits=3)) #13.17081

# Raw data
y <- prostate$lpsa
X <- as.matrix(prostateNormalized[,1:8])
# Centered data
#y <- prostateCentered$lpsa
#X <- as.matrix(prostateCentered[,1:8])

### LASSO as seen in class
library(glmnet)
lasso_glmnet <- glmnet(X, y, alpha = 1, lambda = 0.05) #alpha = 1 for LASSO
lasso_Master <- as.matrix(round(coef(lasso_glmnet),4))
t(lasso_Master)

eps <- 1.e-6
lambda <- 7.2

ptm <- proc.time()
out <- shooting_fu(X,y,lambda)
proc.time() - ptm

ptm <- proc.time()
out2 <- shooting_fu2(X,y,lambda)
proc.time() - ptm

ptm <- proc.time()
lambdaSCD <- 0.0005
out3 <- shooting_scd(X,y,lambdaSCD)
proc.time() - ptm

beta_fu  <- out$beta
beta_fu2 <- out2$beta
beta_scd <- out3$beta

lassoTable <- matrix(c(tail(lasso_Master,-1), beta_fu, beta_fu2, beta_scd),
                     ncol=4)
row.names(lassoTable) <- names(beta_fu)
colnames(lassoTable) <- c("glmnet","shooting","shooting2","scd")
print(lassoTable)

rm(list=ls())
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

signo <- function(w, t){
  out <- NULL
  if(w > t){
    out <- w-t
  } else if (w < -t){
    out <- w+t
  } else {
    out  <- 0
  }
  return(out)
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
    w_j_g_j <- (beta_hat[j]-g_j)
    beta_hat[j] <- signo(w_j_g_j, lambda)

    #}
    step <- step + 1
    converged <- euclidnorm(beta_hat-beta_old) < eps
  }
  output <- list(beta = beta_hat, step = step, converged = converged)
  return(output)
}

#' Shooting SCD algorithm (Bradley, Kyrola, Bickson and Guestrin)
#' @param X design matrix
#' @param y vector of responses
#' @param lambda tuning parameter
#' @param eps tolerance
#' @return The vector of parameters, steps performed and indicator of convergence
shooting_scd2 <- function(X, y, lambda, eps = 1e-6, max_steps = 1000){
  x = cbind(X, -X)
  p <- ncol(X)
  m <- nrow(X)

  converged <- FALSE
  step <- 0
  beta_hat <- rep(0,p)

  while (!converged & (step < max_steps)){

    beta_old <- beta_hat

    # Calculate
    j <- sample(1:p,1)

    tmp <- 0
    for (i in 1:m){
      tmp <- tmp + (sum(beta_hat*X[i,])-y[i])*X[i,j]
    }

    g_j <- tmp/m

    delta_j <- max(-beta_hat[j], -g_j)

    beta_hat[j] <- beta_hat[j] + delta_j

    step <- step + 1
    converged <- euclidnorm(beta_hat-beta_old) < eps
  }
  output <- list(beta = beta_hat, step = step, converged = converged)
  return(output)
}

#' Shooting SCD algorithm as in Fu with stochastic selection
#' @param X design matrix
#' @param y vector of responses
#' @param lambda tuning parameter
#' @param eps tolerance
#' @return The vector of parameters, steps performed and indicator of convergence
shooting_fu_scd <- function(X, y, lambda, eps = 1e-6, max_steps = 1000){

  p <- ncol(X)
  # Products done just once for reuse
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
    j <- sample(1:p,1)
    # shoot at 0
    S_0 <- sum(two_XX[j,] %*% beta_hat) - two_Xy[j] - beta_hat[j] * two_XX[j,j]
    if (S_0 > lambda){
      beta_hat[j] <- (lambda - S_0)/two_XX[j,j]
    } else if (S_0 < -lambda) {
      beta_hat[j] <- (-lambda - S_0)/two_XX[j,j]
    } else {
      beta_hat[j] <- 0
    }


    step <- step + 1
    converged <- euclidnorm(beta_hat-beta_old) < eps
  }
  # Intercept <- mean(y-crossprod(t(X),beta_hat))
  output <- list(beta = beta_hat, step = step, converged = converged)
  return(output)
}

prepareData <- function(data){
  return(scale(data))
}

#myF <- function(X,y,beta, lambda){
#  return(1/2*(euclidnorm(X%*%beta-y))^2 - lambda*sum(abs(beta)))
#}

# Raw data
#library("lasso2", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
library(lasso2)
data(Prostate)
prostate <- Prostate
# Normalized data
prostateNormalized <- as.data.frame(prepareData(prostate[,1:8]))
lpsa <- prostate[,9]

y <- prostate$lpsa
X <- as.matrix(prostateNormalized[,1:8])

outSCD <- shooting_scd(X,y,0.005)
outSCD$beta

outSCD2 <- shooting_scd2(X,y,0.005)
outSCD2$beta

outSCD3 <- shooting_fu_scd(X,y,7.2)
outSCD3$beta

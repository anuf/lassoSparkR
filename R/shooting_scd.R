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

#' Shooting SCD algorithm
#' @param X design matrix
#' @param y vector of responses
#' @param lambda tuning parameter
#' @param eps tolerance
#' @return The vector of parameters, steps performed and indicator of convergence
shooting_scd <- function(X, y, lambda, eps = 1e-6, max_steps = 1000){

  p <- ncol(X)
  m <- nrow(X)
  # Products done just once for reuse
  two_XX = 2*t(X)%*%X # Equivalent: 2*crossprod(X)
  two_Xy = 2*crossprod(X,y)

  converged <- FALSE
  step <- 0
  beta_hat <- rep(0,p)

  while (!converged & (step < max_steps)){

    beta_old <- beta_hat

    # Calculate
    j <- sample(1:p,1)
    #for(j in 1:p){
    print(j)
    # shoot
    #g_j <- sum(two_XX[j,] %*% beta_hat) - two_Xy[j] - beta_hat[j] * two_XX[j,j]
    tmp <- 0
    for (i in 1:m){
      tmp <- tmp + (crossprod(beta_hat,X[i,])-y[i])*X[i,j]
    }
    print(tmp)
    g_j <- tmp/m

    print(g_j)
    w_j_g_j <- (beta_hat[j]-g_j)
    beta_hat[j] <- signo(w_j_g_j, lambda)

    #}
    step <- step + 1
    #converged <- euclidnorm(beta_hat-beta_old) < eps
    print(beta_old)
    print(beta_hat)
    print(myF(X,y,beta_hat,lambda))
    print(myF(X,y,beta_old,lambda))
    converged <- sum(abs(myF(X,y,beta_hat,lambda)-myF(X,y,beta_old,lambda))) < eps
  }
  output <- list(beta = beta_hat, step = step, converged = converged)
  return(output)
}

prepareData <- function(data){
  return(scale(data))
}

myF <- function(X,y,beta, lambda){
  return(1/2*(euclidnorm(X%*%beta-y))^2 - lambda*sum(abs(beta)))
}
# Raw data
library("lasso2", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
data(Prostate)
prostate <- Prostate
# Normalized data
prostateNormalized <- as.data.frame(prepareData(prostate[,1:8]))
lpsa <- prostate[,9]

y <- prostate$lpsa
X <- as.matrix(prostateNormalized[,1:8])

out3 <- shooting_scd(X,y,7.2)
out3

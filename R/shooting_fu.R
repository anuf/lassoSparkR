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

#' Shooting algorithm as in Fu
#' @param X design matrix
#' @param y vector of responses
#' @param lambda tuning parameter
#' @param eps tolerance
#' @return The vector of parameters, steps performed and indicator of convergence
shooting_fu2 <- function(X, y, lambda, eps = 1e-6, max_steps = 1000){

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
    converged <- euclidnorm(beta_hat-beta_old) < eps
  }
  output <- list(beta = beta_hat, step = step, converged = converged)
  return(output)
}

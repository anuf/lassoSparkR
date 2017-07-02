rm(list=ls())

prepareData <- function(data){
  return(scale(data))
}

spark_path <- strsplit(system("brew info apache-spark",intern=T)[6],' ')[[1]][1] # Get your spark path
.libPaths(c(file.path(spark_path,"libexec", "R", "lib"), .libPaths())) # Navigate to SparkR folder
library(SparkR) # Load the library
sparkR.session(master = "local[*]", sparkConfig = list(spark.driver.memory = "2g"))
###############################################################

# load data (Prostate data set)
#library(ElemStatLearn) # Prostate data set
library("lasso2", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
data(Prostate)

# skip last column
#prostate=prostate[,-ncol(prostate)]
prostate <- Prostate
# Normalized data
prostateNormalized <- as.data.frame(prepareData(prostate[,1:8]))
lpsa <- prostate[,9]
localProstata <- as.data.frame(cbind(prostateNormalized,lpsa))
df <- as.DataFrame(localProstata)

# Raw data
y <- prostate$lpsa
X <- prepareData(prostate[,1:8])
two_XX = 2*t(X)%*%X # Equivalent: 2*crossprod(X)
two_Xy = 2*crossprod(X,y)
# Start with $\beta_OLS$
myFormula <- as.formula(paste("y~",paste(colnames(X),collapse ="+")))
OLS <- lm(myFormula, data=as.data.frame(X))
beta_hat <- tail(OLS$coefficients,-1)
shooting_fu2 <- function(X, y, lambda, eps = 1e-6, max_steps = 1000){
  
  p <- ncol(X)
  # Products done just once for reuse
  
  
  
  
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
  # Intercept <- mean(y-crossprod(t(X),beta_hat))
  output <- list(coefficients = beta_hat, step = step, converged = converged)
  return(output)
}
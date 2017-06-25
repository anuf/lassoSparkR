# LASSO sparklyr ml_lib
# fom http://spark.rstudio.com/mllib.html#example-workflow

rm(list=ls())
#' Normalizes data columwise (R-version)
#' @param data A matrix to standarize
#' @return The standarized matrix
prepareData <- function(data){
  return(scale(data))
}


library("lasso2", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
data(Prostate)

# skip last column
#prostate=prostate[,-ncol(prostate)]
prostate <- Prostate
rm(Prostate)
# Normalized data
prostateNormalized <- prepareData(prostate[,1:8])
lpsa <- prostate[,9]
localProstata <- as.data.frame(cbind(prostateNormalized,lpsa))

#install.packages("sparklyr")
library(sparklyr)
# Define spark connection
sc <- spark_connect(master = "local")

library(dbplyr)
prostate_tbl <- copy_to(sc, localProstata)
sparkLR <- ml_linear_regression(lpsa~., data=prostate_tbl, alpha = 1, lambda = 0.05)
coef(sparkLR)

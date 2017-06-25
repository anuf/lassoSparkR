# LASSO sparklyr ml_lib
# fom http://spark.rstudio.com/mllib.html#example-workflow

install.packages("sparklyr")
library(sparklyr)
# Define spark context
sc <- spark_connect(master = "local")

library(dbplyr)
prostate_tbl <- copy_to(sc, Prostate)
ml_linear_regression(lpsa~., data=prostate_tbl, alpha = 1, lambda = 0.05)

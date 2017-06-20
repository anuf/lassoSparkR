# lassoSparkR

A R package that implements several algorithms of the LASSO method for estimation in linear models to be run in Apache Spark engine.

LASSO stands for 'Least Absolute Shrinkage ans Selection Operator'<sup>1</sup>.

Three implementations will be covered:
* Shooting algorithm<sup>2</sup>
* Secuential Stochastic Coordinate Descent (SCD)<sup>3</sup>
* Parallel Stochastic Coordinate Descent (Shotgun)<sup>3</sup>

[1] Tibshirani, R., "Regression Shrinkage and Selection via the Lasso" [1996]

[2] Fu, W. J., "Penalized Regressions: The Bridge Versus the Lasso" [1998]

[3] Bradley, J. C., Kyrola, A., Bickson D., Guestrin, C., "Parallel Coordinate Descent for L1-Regularized Loss Minimization"  [2011]

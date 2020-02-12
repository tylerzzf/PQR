# PQR
code for Shrinkage Quantile Regression for Panel Data with Multiple Structural Breaks

Here we provide a copy of code for Shrinkage Quantile Regression for Panel Data with Multiple Structural Breaks.
The most important function file is `mult-quantile.R`, which contains three functions, `rq.fit.lasso2`，`PQRcpt.func`，`IC.func`. 

`PQRcpt.func` is used to estimate regression coefficients, individual effects and structural breaks. 
The input variable `X` of PQRcpt.func is designed as $$ (x_{1,1}^\top,x_{1,2}^\top,\cdots,x_{1,T}^\top, \cdots, x_{N,T}^\top)^\top $$.
The input variable `y` is the corresponding variable, which is a $$ NT \times 1 $$ vector.  
The input variable `taus` denotes the quantile level(s).
The output of `PQRcpt.func` contains `coef`, `error`, `cpt` and `fix`, which denotes the estimation of regression coefficients, 
value of check loss function, position of structural break(s) and estimation of fixed effects, respectively.

`IC.func` is used to choose tuning parameter from the input.

`DGP.R`, `hausdorff.R` and `nobreak.sim.func` are designed to generate data, compute hausdorff distance and percentage of correctly detect 
structural break.

`b0.R` is main file to show our simulation results of Scenario 2 in the paper.

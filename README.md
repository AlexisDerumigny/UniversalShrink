# How to install

The development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("AlexisDerumigny/UniversalShrink")
```


# 1. Functions for estimation of the covariance matrix


- `cov_analytical_NL_shrinkage()` and `cov_quadratic_inverse_shrinkage()`: perform
  estimation of the covariance matrix using non-linear shrinkage. Both estimators
  are optimal for the Frobenius norm (asymptotically).



# 2. Functions for estimation of the precision matrix


## 2.1. Moore-Penrose-type estimators

- `Moore_Penrose()`: Moore-Penrose estimator of the precision matrix, obtained
  by computing the Moore-Penrose inverse of the sample covariance matrix.

- `Moore_Penrose_shrinkage()` and `Moore_Penrose_shrinkage_toIP`: perform a
  first-order shrinkage of the Moore-Penrose estimator of the precision matrix,
  respectively towards a general target and towards the identity matrix.
  
- `Moore_Penrose_higher_order_shrinkage()`: estimate the precision matrix via
  a polynomial in the Moore-Penrose estimator of the precision matrix.


## 2.2. Ridge-type estimators

- `ridge_target_identity()`, `ridge_target_identity_semioptimal()` and 
  `ridge_target_identity_optimal` perform first-order shrinkage of the Ridge
  estimator towards the identity matrix.
  
- `ridge_target_general()`, `ridge_target_general_semioptimal()`,
  `ridge_target_general_optimal()`  perform first-order shrinkage of the Ridge
  estimator towards a general target.
  
TODO: 
- implement an optimized in $t$ version of `ridge_no_shrinkage`

TODO: clean documentation of the following functions:
- `ridge_no_shrinkage`
- `ridge_higher_order_shrinkage`
- `ridge_higher_order_shrinkage_optimal`
- `ridge_shrinkage_rescaled_optimal`


## 2.3. Moore-Penrose-Ridge "hybrid" estimators

- `MPR_no_shrinkage()`

- `MPR_target()`


# 3. Functions for estimation of optimal portfolio weights


- From a given precision matrix, one can get the optimal portfolio weights
  (*GMV* for Global Minimum Variance) using the function `GMV_PlugIn()`. This can
  also be used by giving as input an estimator of the precision matrix,
  which then will give as output an estimator of the optimal portfolio weights.

- `GMV_Moore_Penrose()`: using the Moore-Penrose inverse `Moore_Penrose()` of
  the precision matrix as an input for the plug-in estimation.

- `GMV_Moore_Penrose_target()`: performs a first-order shrinkage of the
  estimator given by `GMV_Moore_Penrose()` towards an arbitrary (fixed) portfolio
  such as the equally weighted portfolio.




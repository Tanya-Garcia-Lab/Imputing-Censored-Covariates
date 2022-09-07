# Correcting conditional mean imputation for censored covariates and improving usability
## Lotspeich, Grosser & Garcia (2021+)
### Statistical methods to impute censored covariates. 

# Installation 

The `R` package `imputeCensoRd` which implements statistical methods to impute censored covariates can be installed using the `devtools` package in `R` as follows. 

```{r}
# Run once: install.packages("devtools")
devtools::install_github(repo = "Tanya-Garcia-Lab/Imputing-Censored-Covariates", subdir = "imputeCensoRd")
library(imputeCensoRd)
```

# Example
## Simulate Data
![](images/Sim-Setup.png)

```{r}
set.seed(114)

# Set parameters 
N <- 1000
lambda <- -2
beta0 <- 1
beta1 <- 1
beta2 <- 0.25

# Simulate data
z <- rbinom(n = N, size = 1, prob = 0.25)
```

There is a built-in function, `imputeCensoRd::cox_simulation()` which generates the covariate X as described. It takes in the following parameters and is implemented for our simulation in the code chunk below: 

- `n`: sample size
- `logHR`: log hazard ratio coefficients for linear predictor of Cox model 
- `A`: matrix of auxiliary covariates for linear predictor of Cox model
- `dist`: desired distribution, with choices `"Exponential"`, `"Weibull"`, or `"Gompertz"`. Default is `"Exponential"`.
- `lambda`: (<img src="https://render.githubusercontent.com/render/math?math=\lambda">) For `"Exponential"`, `"Weibull"`, and `"Gompertz"`, `lambda` is the scale parameter. Must be positive. Default is `1`.
- `nu`: (<img src="https://render.githubusercontent.com/render/math?math=\nu">) For `"Weibull"`, `nu` is the shape parameter. Must be positive. Default is `NULL`.
- `alpha`: (<img src="https://render.githubusercontent.com/render/math?math=\alpha">) For `"Gompertz"`, `alpha` is the shape parameter. Default is `NULL`.

The baseline density functions for the three `dist` options are shown below.

- `Exponential`: <img src="https://render.githubusercontent.com/render/math?math=f_0(x) = \lambda \exp(-\lambda x)">
- `Weibull`: <img src="https://render.githubusercontent.com/render/math?math=f_0(x) = \lambda \nu \x^{\nu - 1} \exp(-\lambda x^{\nu})">
- `Gompertz`: <img src="https://render.githubusercontent.com/render/math?math=f_0(x) = \lambda \exp(\alpha x) \exp(\frac{\lambda}{\alpha}(1 - \exp(\alpha x)))">

```{r}
x <- imputeCensoRd::cox_simulation(n = N, logHR = lambda, A = matrix(z, ncol = 1), dist = "Exponential", lambda = 5)
e <- rnorm(n = N, mean = 0, sd = 1)
y <- beta0 + beta1 * x + beta2 * z + e
c <- rexp(n = N, rate = 4)
delta <- as.numeric(x <= c)
w <- pmin(x, c)
x[delta == 0] <- NA
sim_dat <- data.frame(y, x, w, z, delta)
```

## Single Imputation

The function `imputeCensoRd::condl_mean_impute()` imputes censored covariates with their conditional mean given censored value and additional covariates (where supplied). This is conditional mean single imputation. We can use it to impute censored `x` in the simulated data and then fit the model for `y ~ x + z` to the imputed dataset. This function takes in the following parameters: 

- `fit`: A `coxph` or `survfit` imputation model object (from the `survival` package).
- `obs`: String column name for the censored covariate.
- `event`: String column name for the censoring indicator of the covariate.
- `addl_covar`: (Optional) string or vector of strings for the additional fully-observed covariates. Default is `NULL`.
- `data`: Datafrane containing columns `obs`, `event`, and (if provided) `addl_covar`.
- `approx_beyond`: Choice of approximation used to extrapolate the survival function beyond the last observed event time. Default is `"expo"` for the exponential approximation from Brown, Hollander, and Kowar (1974). Other choices include `"zero"`, which immediately goes to zero (Efron, 1967), or `"carryforward"`, which carries forward the survival at last event time (Gill, 1980).
- `sample_lambda': (Optional) A logical value. If TRUE, then lambda will be randomly sampled from its estimated asymptotic distribution according to the Cox model estimate.

```{r}
# Fit the imputation model for x ~ z 
imp_mod <- survival::coxph(formula = survival::Surv(time = w, event = delta) ~ z, data = sim_dat)
# Impute censored x in sim_dat
sim_dat_imp <- imputeCensoRd::condl_mean_impute(fit = imp_mod, obs = "w", event = "delta", addl_covar = "z", data = sim_dat, approx_beyond = "expo")
```

The single imputation values are illustrated below, where the x-axis is the observed value `t` and the y-axis is the imputed value. Note: for uncensored subjects, there is no need for imputation so observed and imputed are the same. 

![](images/Imputed-Observed-SI.png)

**Fig. 1.** Illustration of conditional mean single imputation values for a censored covariate.

With the imputed dataset, `sim_dat_imp`, we can now fit the desired analysis model. Since outcome `y` is continuous, we fit a normal linear regression model with covariates `imp` (in place of `x`) and `z`. 

```{r}
lm(formula = y ~ imp + z, data = sim_dat_imp)

Call:
lm(formula = y ~ imp + z, data = sim_dat_imp)

Coefficients:
(Intercept)          imp            z  
     0.8448       1.7176       0.5007  
```

While they might offer bias corrections, single imputation approaches like this are known to underestimate the variability due to the imputed values (since they are treated with the same certainty as the actual observed `x` values). This means that statistical inference based on single imputation will be invalid, so we instead turn to the following multiple imputation approach to correct for this. 

## Multiple Imputation 

The single imputation procedure can be broken down into two steps: (1) conditional mean imputation using a semi- or non-parametric survival function and (2) analysis wherein parameter estimatse are obtained with standard complete data procedures. With multiple imputation, we use bootstrap resampling to draw new data in each iteration and then we apply these same two steps to each of the `M` datasets. We then pool these analyses into one set of parameters which are expected to be unbiased and with appropriate variability estimates. This process is briefly illustrated in the following diagram. 

![](images/MI_Diagram.png)
**Fig. 2.** Multiple imputation: Overview of the steps. *This figure was adapted from Figure 1.6. in Buuren, S. (2012). Flexible imputation of missing data. Boca Raton, FL: CRC Press.*

### Imputing and Analyzing

The function `imputeCensoRd::condl_mean_impute_bootstrap()` imputes censored covariates with their conditional mean given censored value and additional covariates (where supplied) using bootstrap resamples from the supplied dataset. This is conditional mean multiple imputation. We can use it to impute censored `x` in `M` bootstrap samples of the simulated data and then fit the model for `y ~ x + z` to the `M` completed datasets. This function takes in the following parameters: 

- `obs`: String column name for the censored covariate.
- `event`: String column name for the censoring indicator of the covariate.
- `addl_covar`: (Optional) string or vector of strings for the additional fully-observed covariates. Default is `NULL`.
- `data`: Datafrane containing columns `obs`, `event`, and (if provided) `addl_covar`.
- `approx_beyond`: Choice of approximation used to extrapolate the survival function beyond the last observed event time. Default is `"expo"` for the exponential approximation from Brown, Hollander, and Kowar (1974). Other choices include `"zero"`, which immediately goes to zero (Efron, 1967), or `"carryforward"`, which carries forward the survival at last event time (Gill, 1980).
- `M`: an integer number of bootstrap samples to be taken from `data`.

```{r}
# Multiple imputation
sim_dat_imp <- imputeCensoRd::condl_mean_impute_bootstrap(obs = "w", event = "delta", addl_covar = "z", 
							  data = sim_dat, approx_beyond = "expo", M = 5)
```

In this case, `sim_dat_imp` is actually a list of length `M` containing the imputed datasets from each imputation. Individual datasets can be accessed as follows: 

```{r}
head(sim_dat_imp[[1]])
```
```{r}
               w          y            x z delta m        hr      surv          imp
421 0.0000681393  0.3608872 0.0000681393 0     1 1 1.0000000 0.9986917 0.0000681393
999 0.0001519918  3.4185226           NA 1     0 1 0.1565266 0.9967313 0.7831441289
21  0.0002868707 -0.6203306 0.0002868707 0     1 1 1.0000000 0.9947710 0.0002868707
22  0.0002868707 -0.6203306 0.0002868707 0     1 1 1.0000000 0.9947710 0.0002868707
23  0.0002868707 -0.6203306 0.0002868707 0     1 1 1.0000000 0.9947710 0.0002868707
324 0.0005863644  1.6800965 0.0005863644 0     1 1 1.0000000 0.9934624 0.0005863644
```

### Pooling the Results

We can now fit a linear model to each imputed dataset in `sim_dat_imp`. The function `fit_lm_to_imputed_list()` fits the `base R` function `lm()` using a user-specified formula to each element in a list of imputed dataframes. The function then pools the results of parameter estimation for each linear model using Rubin's rules. This functions takes the following two arguments:

- `imputed_list`: A list with dataframe elements that have been completed via conditional mean imputation, returned by `condl_mean_impute_bootstrap()`.
- `formula`: A formula object used to fit a linear model to element of `imputed_list`.

```{r}
# Pooling Analysis Results
pooled_lm_res <- imputeCensoRd::fit_lm_to_imputed_list(imputed_list = sim_dat_imp, formula = y ~ imp + z)
```

The function returns a list with the following two vectors each of length `p` (the number of regression parameters in `formula`):

- `Coef`: The average coefficient estimates obtained from fitting `formula` to each dataframe in `imputed_list`
- `Pooled_Var`: The pooled variance estimates, calculated using Rubin's rules

In our simulation, this looks like

```{r}
$Coef
(Intercept)         imp           z 
  0.8348876   1.8622895   0.6744820 

$Pooled_Var
(Intercept)         imp           z 
0.008711573 0.096049517 0.102880558 
```

Finally, from this output we have the following interpretations about the association between the outcome, `y`, censored covariate, `x`, after controlling for fully-observed covariate `z`. If we were to compare two populations of subjects with the same value for `z` but whose `x` values differed by 1, the group with the higher `x` is expected to have a `y` value 1.86 higher, as well. Using `pooled_lm_res$Pooled_Var`, we can estimate the standard error (`se`) for this association, which could be used for significance testing or to build confidence intervals. 

```{r}
se <- sqrt(pooled_lm_res$Pooled_Var)[2]
```

# References

Bender, R., Augustin, T., and Blettner, M. (2005). Generating survival times to simulate Cox proportional hazards models. *Statistics in Medicine*, 24:1713–1723.

Brown, J. B. W., Hollander, M., and Korwar, R. M. (1974) Nonparametric Tests of Independence for Censored Data, with Applications to Heart Transplant Studies. *Reliability and Biometry: Statistical Analysis of Lifelength*, F. Proschan and R.
J. Serfling, eds. Philadelphia: SIAM, pp. 327-354.

Efron, B. (1967) The Two Sample Problem with Censored Data. *Proceedings of the Fifth Berkeley Symposium On Mathematical Statistics and Probability.* New
York: Prentice-Hall, 4:831-853.

Gill, R.D. (1980) Censoring and Stochastic Integrals. *Mathematical Centre Tracts*. Amsterdam: Mathematisch Centrum, 124.

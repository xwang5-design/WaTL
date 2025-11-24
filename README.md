# WaTL: Wasserstein Transfer Learning

<!-- badges: start -->
<!-- badges: end -->

An R package implementing the **Wasserstein Transfer Learning (WaTL)** algorithm for distribution regression, as described in Zhang et al. (2025).

## Overview

The WaTL package provides methods for combining information from target and source datasets to improve prediction accuracy when dealing with distributional responses. The algorithm consists of three main steps:

1. **Weighted Auxiliary Estimator**: Aggregates information from target and source datasets using global Fréchet weights
2. **Bias Correction**: Uses target data to correct bias via regularized optimization
3. **Projection**: Projects to Wasserstein space (enforces monotonicity)

## Installation

You can install the released version of WaTL from CRAN with:

```r
install.packages("WaTL")
```

Or install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("xwang5-design/WaTL")
```

## Usage

### Basic Example

```r
library(WaTL)

# Generate simulated data (Setting 1)
data <- simulate_data_setting1(
  n_t = 200,                           # Target sample size
  n_vec = c(100, 200, 300, 400, 500),  # Source sample sizes
  M = 100,                             # Grid size
  K = 5                                # Number of source datasets
)

# Compute weighted auxiliary estimator (Step 1)
f1_hat <- compute_f1_hat(
  source_data_list = data$source_data,
  X_value = 0.5,
  M = 100,
  X_t = data$X_t,
  Y_t = data$Y_t
)

# Select regularization parameter via cross-validation
cv_results <- cv_search_lambda(
  Y_t = data$Y_t,
  s_vec = s_vec,  # Fréchet weights
  X_t = data$X_t,
  source_d = data$source_data,
  M = 100,
  lambda_grid = seq(0, 3, by = 0.1)
)

# Compute bias-corrected estimator (Step 2)
f_hat <- compute_f_L2(
  Y_t = data$Y_t,
  s_vec = s_vec,
  f1_hat = f1_hat,
  lambda = cv_results$best_lambda,
  M = 100
)
```

### Global Fréchet Regression

```r
# Using the grem function for baseline comparison
n <- 100
x <- runif(n)
y <- lapply(1:n, function(i) rnorm(50, mean = x[i], sd = 0.1))

# Fit global Fréchet regression
fit <- grem(y, x, xOut = c(0.3, 0.5, 0.7))
```

## Original Simulation and Real Data Scripts

The package also includes the original simulation and real data analysis scripts in the `Simulation/` and `RealData/` directories for reproducibility.

## Citation

```
@article{zhang2025wasserstein,
  title={Wasserstein Transfer Learning},
  author={Zhang, Kaicheng and Zhang, Sinian and Zhou, Doudou and Zhou, Yidong},
  journal={arXiv preprint arXiv:2505.17404},
  year={2025}
}
```



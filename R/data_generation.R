#' Simulate Data - Setting 1
#'
#' Generates target and source data according to Setting 1 specifications from
#' the Wasserstein Transfer Learning paper.
#'
#' @details
#' Target Population (k=0):
#' \itemize{
#'   \item X^(0) ~ U(0,1)
#'   \item F^{-1}_{nu^(0)}(u) = w^(0)(1-u)u + (1-X^(0))u + X^(0) F^{-1}_{Z^(0)}(u)
#'   \item Z^(0) ~ N(0.5, 1)|_{(0,1)} (truncated normal)
#'   \item w^(0) ~ N(0, 1)|_{(-0.5, 0.5)} (truncated normal)
#' }
#'
#' Source Population (k=1,...,K):
#' \itemize{
#'   \item psi_k = 0.1k (similarity parameter)
#'   \item X^(k) ~ U(0,1)
#'   \item F^{-1}_{nu^(k)}(u) = w^(k)(1-u)u + (1-X^(k))u + X^(k) F^{-1}_{Z^(k)}(u)
#'   \item Z^(k) ~ N(0.5, 1-psi_k)|_{(0,1)}
#'   \item w^(k) ~ N(0, 1)|_{(-0.5, 0.5)}
#' }
#'
#' @param n_t Target sample size.
#' @param n_vec Vector of source sample sizes (length K).
#' @param M Grid size for quantile functions.
#' @param K Number of source datasets (default 5).
#' @return A list with components:
#'   \item{X_t}{Target predictors (n_t x 1)}
#'   \item{Y_t}{Target quantile functions (n_t x M)}
#'   \item{source_data}{List of K source datasets, each containing X_s, Y_s, and d_k}
#' @export
#' @examples
#' \dontrun{
#' data <- simulate_data_setting1(n_t = 200, n_vec = c(100, 200, 300, 400, 500), M = 100, K = 5)
#' }
simulate_data_setting1 <- function(n_t, n_vec, M, K=5)
{

  # Define quantile function for target distribution Z^(0)
  # Z^(0) ~ N(0.5, 1)|_{(0,1)}
  F_inv_Phi <- function(z) {
    truncnorm::qtruncnorm(z, a=0, b=1, mean=0.5, sd=1)
  }

  # Define quantile function for source distribution Z^(k)
  # Z^(k) ~ N(0.5, 1-psi_k)|_{(0,1)}
  F_inv_phi <- function(z, d_k){
    truncnorm::qtruncnorm(z, a=0, b=1, mean=0.5, sd=sqrt(1 - d_k))
  }

  #============================================================================
  # Generate TARGET data (k=0)
  #============================================================================
  X_t <- runif(n_t, 0, 1)
  w_t <- truncnorm::rtruncnorm(n_t, a=-0.5, b=0.5, mean=0, sd=1)

  z_grid <- seq(1 / M, 1 - 1/ M, length.out=M)

  Y_t <- t(sapply(seq_len(n_t), function(i) {
    w_t[i]*(1-z_grid)*z_grid + (1 - X_t[i])*z_grid + X_t[i]*F_inv_Phi(z_grid)
  }))

  #============================================================================
  # Generate SOURCE data (k=1,...,K)
  #============================================================================
  source_data <- list()
  for(k_idx in seq_len(K)){
    d_k <- k_idx * 0.1
    n_k <- n_vec[k_idx]

    X_s_k <- runif(n_k, 0, 1)
    w_s_k <- truncnorm::rtruncnorm(n_k, a=-0.5, b=0.5, mean=0, sd=1)

    Y_s_k <- t(sapply(seq_len(n_k), function(i) {
      w_s_k[i]*(1-z_grid)*z_grid + (1 - X_s_k[i])*z_grid + X_s_k[i]*F_inv_phi(z_grid, d_k)
    }))

    source_data[[k_idx]] <- list(
      X_s = X_s_k,
      Y_s = Y_s_k,
      d_k = d_k
    )
  }

  return(list(
    X_t = X_t,
    Y_t = Y_t,
    source_data = source_data
  ))
}


#' Simulate Data - Setting 2
#'
#' Generates target and source data according to Setting 2 specifications from
#' the Wasserstein Transfer Learning paper.
#'
#' @param n_t Target sample size.
#' @param n_vec Vector of source sample sizes (length K).
#' @param M Grid size for quantile functions.
#' @param K Number of source datasets (default 5).
#' @return A list with components:
#'   \item{X_t}{Target predictors (n_t x 1)}
#'   \item{Y_t}{Target quantile functions (n_t x M)}
#'   \item{source_data}{List of K source datasets, each containing X_s, Y_s, and d_k}
#' @export
#' @examples
#' \dontrun{
#' data <- simulate_data_setting2(n_t = 200, n_vec = c(100, 200, 300, 400, 500), M = 100, K = 5)
#' }
simulate_data_setting2 <- function(n_t, n_vec, M, K=5){

  X_t <- runif(n_t, -1, 1)
  mu_t <- rnorm(n_t, mean=3*X_t, sd=0.5)

  shape_t <- (3 + 3 * X_t)^2
  rate_t <- 3 + 3 * X_t
  sigma_t <- rgamma(n_t, shape=shape_t, rate=rate_t)


  F_inv_Phi <- function(z){
    truncnorm::qtruncnorm(z, a=0, b=1, mean=0.5, sd=1)
  }

  z_grid <- seq(1 / M, 1 - 1/ M, length.out=M)

  Y_t<- t(sapply(seq_len(n_t), function(i){
    mu_t[i] + sigma_t[i] * F_inv_Phi(z_grid)
  }))

  source_data <- list()
  for(k_idx in seq_len(K))
  {
    d_k <- k_idx * 0.05
    n_k <- n_vec[k_idx]

    X_s_k  <- runif(n_k, -1, 1)
    mu_s_k <- rnorm(n_k, mean=3*X_s_k, sd=0.5)
    shape_s_k <- (3 + d_k + (3 + d_k)*X_s_k)^2
    rate_s_k  <- 3 + d_k + (3 + d_k)*X_s_k
    sigma_s_k <- rgamma(n_k, shape=shape_s_k, rate=rate_s_k)

    Y_s_k <- t(sapply(seq_len(n_k), function(i){
      mu_s_k[i] + sigma_s_k[i] * F_inv_Phi(z_grid)
    }))

    source_data[[k_idx]] <- list(
      X_s = X_s_k,
      Y_s = Y_s_k,
      d_k = d_k
    )
  }

  return(list(
    X_t = X_t,
    Y_t = Y_t,
    source_data = source_data
  ))
}

#' Compute Weighted Auxiliary Estimator (Step 1 of WaTL)
#'
#' Implements the weighted auxiliary estimator that aggregates information
#' from both target (k=0) and source (k=1,...,K) datasets using global Frechet weights.
#'
#' @details
#' The weighted auxiliary estimator is computed as:
#' \deqn{\hat{f}(x) = \frac{1}{n_0 + n_A} \sum_{k=0}^K n_k \hat{f}^{(k)}(x)}
#'
#' where:
#' \itemize{
#'   \item \eqn{\hat{f}^{(k)}(x) = \frac{1}{n_k} \sum_{i=1}^{n_k} s_G^{(k)}(x_i) F^{-1}_{\nu_i^{(k)}}}
#'   \item \eqn{s_G^{(k)}(x) = 1 + (X^{(k)} - \theta_k)^T \Sigma_k^{-1} (x - \theta_k)} (global Frechet weight)
#'   \item \eqn{n_A = \sum_{k=1}^K n_k} (total source sample size)
#' }
#'
#' @param source_data_list List of K source datasets, each containing X_s and Y_s.
#' @param X_value Query point where to evaluate the estimator.
#' @param M Grid size for quantile functions.
#' @param X_t Target predictors (optional, but recommended).
#' @param Y_t Target quantile functions (optional, but recommended).
#' @return The weighted auxiliary estimator evaluated at X_value (vector of length M).
#' @export
#' @examples
#' \dontrun{
#' data <- simulate_data_setting1(n_t = 200, n_vec = c(100, 200, 300), M = 100, K = 3)
#' f1_hat <- compute_f1_hat(data$source_data, X_value = 0.5, M = 100,
#'                          X_t = data$X_t, Y_t = data$Y_t)
#' }
compute_f1_hat <- function(source_data_list,
                       X_value,
                       M,
                       X_t = NULL,
                       Y_t = NULL)
{

  K <- length(source_data_list)

  z_grid <- seq(1 / M, 1 - 1/ M, length.out=M)
  bigF <- rep(0, length(z_grid))
  n_total <- 0

  #============================================================================
  # Include TARGET data (k=0) in the weighted estimator
  #============================================================================
  if (!is.null(X_t) && !is.null(Y_t)) {
    n_t <- length(X_t)
    X_t_mat <- to_matrix(X_t)

    xbar_t <- colMeans(X_t_mat)
    Sigma_t <- cov(X_t_mat) * (nrow(X_t_mat)-1)/nrow(X_t_mat)
    invSigma_t <- solve(Sigma_t)

    s_vec_t <- sapply(1:nrow(X_t_mat), function(i){
      as.numeric(1 + (X_t_mat[i,] - xbar_t) %*% invSigma_t %*% (X_value - xbar_t))
    })

    for(i in seq_len(n_t)){
      bigF <- bigF + s_vec_t[i] * Y_t[i,]
    }
    n_total <- n_total + n_t
  }

  #============================================================================
  # Add SOURCE data (k=1,...,K) contributions
  #============================================================================
  for(k_idx in seq_len(K))
  {
    X_s   <- source_data_list[[k_idx]]$X_s
    Y_s <- source_data_list[[k_idx]]$Y_s
    n_k   <- length(X_s)

    X_s_mat  <- to_matrix(X_s)

    xbar <- colMeans(X_s_mat)
    Sigma <- cov(X_s_mat) * (nrow(X_s_mat)-1)/nrow(X_s_mat)
    invSigma <- solve(Sigma)

    s_vec <- sapply(1:nrow(X_s_mat), function(i){
      as.numeric(1 + (X_s_mat[i,] - xbar) %*% invSigma %*% (X_value - xbar))
    })

    for(i in seq_len(n_k)){
      bigF <- bigF + s_vec[i] * Y_s[i,]
    }
    n_total <- n_total + n_k
  }

  f1 <- bigF / n_total

  return(f1)
}


#' Compute Bias-Corrected Estimator (Step 2 of WaTL)
#'
#' Implements the bias-corrected estimator using gradient descent to minimize
#' a regularized objective function.
#'
#' @details
#' The bias-corrected estimator solves:
#' \deqn{\hat{f}_0(x) = \arg\min_{g \in L^2(0,1)} \frac{1}{n_0} \sum_{i=1}^{n_0} s_i ||F^{-1}_{\nu_i^{(0)}} - g||_2^2 + \lambda ||g - \hat{f}(x)||_2}
#'
#' where:
#' \itemize{
#'   \item First term: Target data fidelity (uses global Frechet weights s_i)
#'   \item Second term: Regularization toward auxiliary estimator from Step 1
#'   \item lambda: Regularization parameter (selected via cross-validation)
#' }
#'
#' @param Y_t Target quantile functions (n_0 x M matrix).
#' @param s_vec Global Frechet weights for target (length n_0).
#' @param f1_hat Weighted auxiliary estimator from Step 1 (length M).
#' @param lambda Regularization parameter.
#' @param M Grid size.
#' @param max_iter Maximum iterations for gradient descent (default 1000).
#' @param step_size Step size for gradient descent (default 0.5).
#' @param tol Convergence tolerance (default 1e-8).
#' @return The bias-corrected estimator (vector of length M).
#' @export
#' @examples
#' \dontrun{
#' data <- simulate_data_setting1(n_t = 200, n_vec = c(100, 200), M = 100, K = 2)
#' f1_hat <- compute_f1_hat(data$source_data, X_value = 0.5, M = 100,
#'                          X_t = data$X_t, Y_t = data$Y_t)
#' X_t_mat <- to_matrix(data$X_t)
#' xbar <- colMeans(X_t_mat)
#' Sigma <- cov(X_t_mat) * (nrow(X_t_mat)-1)/nrow(X_t_mat)
#' invSigma <- solve(Sigma)
#' s_vec <- sapply(1:nrow(X_t_mat), function(i){
#'   as.numeric(1 + (X_t_mat[i,] - xbar) %*% invSigma %*% (0.5 - xbar))
#' })
#' f_hat <- compute_f_L2(data$Y_t, s_vec, f1_hat, lambda = 0.1, M = 100)
#' }
compute_f_L2 <- function(Y_t,
                         s_vec,
                         f1_hat,
                         lambda,
                         M,
                         max_iter = 1000,
                         step_size = 0.5,
                         tol = 1e-8) {
  # Input validation
  n0 <- nrow(Y_t)
  if (length(s_vec) != n0) {
    stop("s_vec length must match row number of Y_t.")
  }
  if (length(f1_hat) != M) {
    stop("f1_hat length must be M.")
  }

  # Initialize at auxiliary estimator
  f <- as.numeric(f1_hat)

  #============================================================================
  # Gradient Descent Optimization
  #============================================================================
  for (iter in 1:max_iter)
  {
    diff_mat <- sweep(Y_t, 2, f, FUN = "-")
    diff_mat <- -diff_mat
    weighted_diff_mat <- sweep(diff_mat, 1, s_vec, FUN = "*")
    gradient_target <- colSums(weighted_diff_mat)

    gradient_reg <- lambda * (f - f1_hat)

    grad <- gradient_target + gradient_reg

    f_new <- f - step_size * grad

    diff_norm <- sqrt(sum((f_new - f)^2))
    if (diff_norm < tol) {
      f <- f_new
      break
    }

    f <- f_new
  }

  return(as.numeric(f))
}

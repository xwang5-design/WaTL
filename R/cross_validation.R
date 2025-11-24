#' Cross-Validation for Lambda Selection
#'
#' Performs k-fold cross-validation to select the optimal regularization parameter
#' lambda for the WaTL algorithm.
#'
#' @param Y_t Target quantile functions (n_0 x M matrix).
#' @param s_vec Global Frechet weights for target (length n_0).
#' @param X_t Target predictors.
#' @param source_d List of source datasets.
#' @param M Grid size for quantile functions.
#' @param lambda_grid Vector of candidate lambda values to test.
#' @param n_folds Number of folds for cross-validation (default 5).
#' @param max_iter Maximum iterations for gradient descent (default 200).
#' @param step_size Step size for gradient descent (default 0.5).
#' @param tol Convergence tolerance (default 1e-6).
#' @param n_cores Number of cores for parallel processing (default NULL, uses detectCores() - 1).
#' @param val_sample_size Validation sample size per fold (default 30).
#' @return A list with components:
#'   \item{results}{Data frame with lambda values and corresponding CV errors}
#'   \item{best_lambda}{The lambda value with minimum CV error}
#'   \item{total_time}{Total time elapsed for cross-validation}
#' @export
#' @examples
#' \dontrun{
#' data <- simulate_data_setting1(n_t = 200, n_vec = c(100, 200), M = 100, K = 2)
#' X_t_mat <- to_matrix(data$X_t)
#' xbar <- colMeans(X_t_mat)
#' Sigma <- cov(X_t_mat) * (nrow(X_t_mat)-1)/nrow(X_t_mat)
#' invSigma <- solve(Sigma)
#' s_vec <- sapply(1:nrow(X_t_mat), function(i){
#'   as.numeric(1 + (X_t_mat[i,] - xbar) %*% invSigma %*% (0.5 - xbar))
#' })
#' cv_out <- cv_search_lambda(Y_t = data$Y_t, s_vec = s_vec, X_t = data$X_t,
#'                            source_d = data$source_data, M = 100,
#'                            lambda_grid = seq(0, 1, by = 0.1))
#' }
cv_search_lambda <- function(Y_t,
                             s_vec,
                             X_t,
                             source_d,
                             M,
                             lambda_grid = NULL,
                             n_folds = 5,
                             max_iter = 200,
                             step_size = 0.5,
                             tol = 1e-6,
                             n_cores = NULL,
                             val_sample_size = 30
)
{
  start_time <- Sys.time()

  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }

  n <- nrow(Y_t)
  total_lambdas <- length(lambda_grid)

  cat(sprintf("\nStarting cross-validation at %s\n", format(start_time)))
  cat(sprintf("Total lambdas to test: %d\n", total_lambdas))
  cat(sprintf("Number of cores: %d\n", n_cores))
  cat(sprintf("Validation sample size: %d\n", val_sample_size))
  cat("----------------------------------------\n\n")

  shuffle_idx <- sample.int(n)
  fold_id <- cut(seq_len(n), breaks = n_folds, labels = FALSE)

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  pb <- utils::txtProgressBar(max=total_lambdas, style=3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress=progress)

  parallel::clusterExport(
    cl,
    c("compute_f1_hat",
      "compute_f_L2",
      "to_matrix",
      "compute_L2_distance",
      "Y_t",
      "s_vec",
      "X_t",
      "source_d",
      "M",
      "lambda_grid",
      "n_folds",
      "max_iter",
      "step_size",
      "tol",
      "val_sample_size",
      "shuffle_idx",
      "fold_id"
    ),
    envir = environment()
  )

  progress_counter <- 0

  results <- foreach::foreach(l_idx = seq_along(lambda_grid),
                     .combine = 'rbind',
                     .options.snow = opts) %dopar% {
                       lambda_val <- lambda_grid[l_idx]
                       fold_errors <- numeric(n_folds)

                       for (k in seq_len(n_folds)) {
                         val_idx_full <- shuffle_idx[fold_id == k]

                         if (length(val_idx_full) > val_sample_size) {
                           val_idx <- sample(val_idx_full, val_sample_size)
                         } else {
                           val_idx <- val_idx_full
                         }

                         train_idx <- shuffle_idx[fold_id != k]

                         Y_train <- Y_t[train_idx, , drop = FALSE]
                         s_vec_train <- s_vec[train_idx]
                         X_train <- X_t[train_idx, ]

                         Y_val <- Y_t[val_idx, , drop = FALSE]
                         s_vec_val <- s_vec[val_idx]
                         X_val <- X_t[val_idx, ]

                         val_errors_fold <- numeric(length(val_idx))

                         for (i in seq_along(val_idx)) {
                           X_value <- as.vector(as.numeric(X_val[i,]))

                           f1_hat_i <- compute_f1_hat(
                             source_data_list = source_d,
                             X_value = X_value,
                             M = M,
                             X_t = X_train,
                             Y_t = Y_train
                           )

                           f_train <- tryCatch({
                             compute_f_L2(
                               Y_t = Y_train,
                               s_vec = s_vec_train,
                               f1_hat = f1_hat_i,
                               lambda = lambda_val,
                               M = M,
                               max_iter = max_iter,
                               step_size = step_size,
                               tol = tol
                             )
                           }, error = function(e) NULL)

                           if (!is.null(f_train)) {
                             val_errors_fold[i] <- mean((f_train - Y_val[i,])^2)
                           } else {
                             val_errors_fold[i] <- NA
                           }
                         }

                         valid_errors <- !is.na(val_errors_fold)
                         if (sum(valid_errors) > 0) {
                           fold_errors[k] <- mean(val_errors_fold[valid_errors])
                         } else {
                           fold_errors[k] <- NA
                         }
                       }

                       c(lambda_val, mean(fold_errors, na.rm = TRUE))
                     }

  parallel::stopCluster(cl)
  close(pb)

  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "mins")

  cat("\n\n----------------------------------------\n")
  cat(sprintf("Cross-validation completed at %s\n", format(end_time)))
  cat(sprintf("Total time elapsed: %.2f minutes\n", as.numeric(total_time)))
  cat("----------------------------------------\n")

  results <- as.data.frame(results)
  colnames(results) <- c("lambda", "cv_error")

  best_lambda <- results$lambda[which.min(results$cv_error)]
  best_error <- min(results$cv_error)

  return(list(
    results = results,
    best_lambda = best_lambda,
    total_time = total_time
  ))
}

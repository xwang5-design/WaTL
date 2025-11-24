#' @keywords internal
#' @importFrom stats cov var runif rnorm rgamma weighted.mean
#' @importFrom foreach %dopar%
#' @importFrom parallel makeCluster stopCluster clusterExport detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils txtProgressBar setTxtProgressBar globalVariables
#' @importFrom Matrix Matrix
#' @importFrom osqp osqp osqpSettings solve_osqp
#' @importFrom truncnorm qtruncnorm rtruncnorm
"_PACKAGE"

# Fix for foreach global variable binding NOTE
utils::globalVariables(c("l_idx"))

#' WaTL: Wasserstein Transfer Learning
#'
#' The WaTL package implements the Wasserstein Transfer Learning algorithm for
#' distribution regression. It provides methods for combining information from
#' target and source datasets to improve prediction accuracy when dealing with
#' distributional responses.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{compute_f1_hat}}: Compute weighted auxiliary estimator (Step 1)
#'   \item \code{\link{compute_f_L2}}: Compute bias-corrected estimator (Step 2)
#'   \item \code{\link{cv_search_lambda}}: Cross-validation for lambda selection
#'   \item \code{\link{simulate_data_setting1}}: Generate data for Setting 1
#'   \item \code{\link{simulate_data_setting2}}: Generate data for Setting 2
#'   \item \code{\link{grem}}: Global regression for empirical measures
#' }
#'
#' @section References:
#' Zhang, K., Zhang, S., Zhou, D., & Zhou, Y. (2025).
#' Wasserstein Transfer Learning. arXiv preprint arXiv:2505.17404.
#'
#' @name WaTL-package
NULL

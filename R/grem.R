#' Global Regression for Empirical Measures (GREM)
#'
#' Performs global Frechet regression using the Wasserstein metric for
#' distributional responses.
#'
#' @param y A list of empirical distributions (each element is a vector of observations).
#' @param x A matrix, data frame, or vector of predictors.
#' @param xOut Optional matrix of out-of-sample predictor values for prediction.
#' @param optns Optional list of control parameters (e.g., upper, lower bounds).
#' @return An object of class "rem" containing:
#'   \item{qf}{Fitted quantile functions (n x M matrix)}
#'   \item{qfSupp}{Support points for quantile functions}
#'   \item{qp}{Predicted quantile functions at xOut (if provided)}
#'   \item{qpSupp}{Support points for predictions}
#'   \item{RSquare}{R-squared value}
#'   \item{adjRSquare}{Adjusted R-squared value}
#'   \item{residuals}{Vector of residuals}
#'   \item{y}{Input empirical distributions}
#'   \item{x}{Input predictors}
#'   \item{xOut}{Out-of-sample predictors}
#'   \item{optns}{Options used}
#' @export
#' @examples
#' \dontrun{
#' # Generate some data
#' n <- 100
#' x <- runif(n)
#' y <- lapply(1:n, function(i) rnorm(50, mean = x[i], sd = 0.1))
#' fit <- grem(y, x, xOut = c(0.3, 0.5, 0.7))
#' }
grem <- function(y = NULL,
                 x = NULL,
                 xOut = NULL,
                 optns = list()) {
  if (is.null(y) | is.null(x)) {
    stop("requires the input of both y and x")
  }
  if (!is.matrix(x)) {
    if (is.data.frame(x) | is.vector(x)) {
      x <- as.matrix(x)
    } else {
      stop("x must be a matrix or a data frame or a vector")
    }
  }
  n <- nrow(x)
  p <- ncol(x)
  if (!is.list(y)) {
    stop("y must be a list")
  }
  if (length(y) != n) {
    stop("the number of rows in x must be the same as the number of empirical measures in y")
  }
  if (!is.null(xOut)) {
    if (!is.matrix(xOut)) {
      if (is.data.frame(xOut)) {
        xOut <- as.matrix(xOut)
      } else if (is.vector(xOut)) {
        if (p == 1) {
          xOut <- as.matrix(xOut)
        } else {
          xOut <- t(xOut)
        }
      } else {
        stop("xOut must be a matrix or a data frame or a vector")
      }
    }
    if (ncol(xOut) != p) {
      stop("x and xOut must have the same number of columns")
    }
    nOut <- nrow(xOut)
  } else {
    nOut <- 0
  }

  N <- sapply(y, length)
  y <- lapply(1:n, function(i) {
    sort(y[[i]])
  })

  M <- min(plcm(N), n * max(N), 5000)
  yM <- t(sapply(1:n, function(i) {
    residual <- M %% N[i]
    if(residual) {
      sort(c(rep(y[[i]], each = M %/% N[i]), sample(y[[i]], residual)))
    } else {
      rep(y[[i]], each = M %/% N[i])
    }
  }))

  # initialization of OSQP solver
  A <- cbind(diag(M), rep(0, M)) + cbind(rep(0, M), -diag(M))
  if (!is.null(optns$upper) &
      !is.null(optns$lower)) {
    l <- c(optns$lower, rep(0, M - 1), -optns$upper)
  } else if (!is.null(optns$upper)) {
    A <- A[, -1]
    l <- c(rep(0, M - 1), -optns$upper)
  } else if (!is.null(optns$lower)) {
    A <- A[, -ncol(A)]
    l <- c(optns$lower, rep(0, M - 1))
  } else {
    A <- A[, -c(1, ncol(A))]
    l <- rep(0, M - 1)
  }
  P <- diag(M)
  A <- t(A)
  q <- rep(0, M)
  u <- rep(Inf, length(l))
  model <-
    osqp::osqp(
      P = P,
      q = q,
      A = A,
      l = l,
      u = u,
      osqp::osqpSettings(max_iter = 1e05, eps_abs = 1e-05, eps_rel = 1e-05, verbose = FALSE)
    )

  xMean <- colMeans(x)
  invVa <- solve(var(x) * (n - 1) / n)
  wc <-
    t(apply(x, 1, function(xi) {
      t(xi - xMean) %*% invVa
    }))
  if (nrow(wc) != n) {
    wc <- t(wc)
  }

  qf <- matrix(nrow = n, ncol = M)
  residuals <- rep.int(0, n)
  totVa <- sum((scale(yM, scale = FALSE))^2) / M
  for (i in 1:n) {
    w <- apply(wc, 1, function(wci) {
      1 + t(wci) %*% (x[i, ] - xMean)
    })
    qNew <- apply(yM, 2, weighted.mean, w)
    if (any(w < 0)) {
      model$Update(q = -qNew)
      qNew <- sort(model$Solve()$x)
    }
    if (!is.null(optns$upper)) {
      qNew <- pmin(qNew, optns$upper)
    }
    if (!is.null(optns$lower)) {
      qNew <- pmax(qNew, optns$lower)
    }
    qf[i, ] <- qNew
    residuals[i] <- sqrt(sum((yM[i, ] - qf[i, ])^2) / M)
  }
  qfSupp <- 1:M / M
  resVa <- sum(residuals^2)
  RSquare <- 1 - resVa / totVa
  adjRSquare <- RSquare - (1 - RSquare) * p / (n - p - 1)

  if (nOut > 0) {
    qp <- matrix(nrow = nOut, ncol = M)
    for (i in 1:nOut) {
      w <- apply(wc, 1, function(wci) {
        1 + t(wci) %*% (xOut[i, ] - xMean)
      })
      qNew <- apply(yM, 2, weighted.mean, w)
      if (any(w < 0)) {
        model$Update(q = -qNew)
        qNew <- sort(model$Solve()$x)
      }
      if (!is.null(optns$upper)) {
        qNew <- pmin(qNew, optns$upper)
      }
      if (!is.null(optns$lower)) {
        qNew <- pmax(qNew, optns$lower)
      }
      qp[i, ] <- qNew
    }
    qpSupp <- 1:M / M

    res <-
      list(
        qf = qf,
        qfSupp = qfSupp,
        qp = qp,
        qpSupp = qpSupp,
        RSquare = RSquare,
        adjRSquare = adjRSquare,
        residuals = residuals,
        y = y,
        x = x,
        xOut = xOut,
        optns = optns
      )
  } else {
    res <- list(
      qf = qf,
      qfSupp = qfSupp,
      RSquare = RSquare,
      adjRSquare = adjRSquare,
      residuals = residuals,
      y = y,
      x = x,
      optns = optns
    )
  }

  class(res) <- "rem"
  res
}

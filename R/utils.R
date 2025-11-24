#' Convert to Matrix
#'
#' Helper function to convert vectors or data frames to matrices.
#'
#' @param x A vector, data frame, or matrix.
#' @return A matrix representation of x.
#' @keywords internal
#' @export
to_matrix <- function(x) {
  if (is.vector(x)) {
    matrix(x, ncol=1)
  } else if (is.data.frame(x)) {
    as.matrix(x)
  } else {
    x
  }
}

#' Compute L2 Distance
#'
#' Computes the L2 distance between two vectors.
#'
#' @param vec1 First vector.
#' @param vec2 Second vector.
#' @return The L2 distance between vec1 and vec2.
#' @keywords internal
#' @export
compute_L2_distance <- function(vec1, vec2){
  sqrt(mean( (vec1 - vec2)^2 ))
}

#' Practical Least Common Multiple
#'
#' Computes the least common multiple of a vector of integers,
#' with a practical upper bound.
#'
#' @param x A numeric vector of integers.
#' @return The least common multiple, bounded by practical limits.
#' @keywords internal
#' @export
plcm <- function(x) {
  stopifnot(is.numeric(x))

  x <- x[x != 0]
  n <- length(x)
  if (n == 0) {
    l <- 0
  } else if (n == 1) {
    l <- x
  } else if (n == 2) {
    l <- lcm(x[1], x[2])
  } else {
    l <- lcm(x[1], x[2])
    for (i in 3:n) {
      l <- lcm(l, x[i])
    }
  }
  return(l)
}

#' Least Common Multiple
#'
#' Computes the least common multiple of two integers.
#'
#' @param n First integer.
#' @param m Second integer.
#' @return The least common multiple of n and m.
#' @keywords internal
#' @export
lcm <- function(n, m) {
  stopifnot(is.numeric(n), is.numeric(m))
  if (length(n) != 1 || floor(n) != ceiling(n) ||
      length(m) != 1 || floor(m) != ceiling(m)) {
    stop("Arguments 'n', 'm' must be integer scalars.")
  }
  if (n == 0 && m == 0) {
    return(0)
  }

  return(n / gcd(n, m) * m)
}

#' Greatest Common Divisor
#'
#' Computes the greatest common divisor of two integers.
#'
#' @param n First integer.
#' @param m Second integer.
#' @return The greatest common divisor of n and m.
#' @keywords internal
#' @export
gcd <- function(n, m) {
  stopifnot(is.numeric(n), is.numeric(m))
  if (length(n) != 1 || floor(n) != ceiling(n) ||
      length(m) != 1 || floor(m) != ceiling(m)) {
    stop("Arguments 'n', 'm' must be integer scalars.")
  }
  if (n == 0 && m == 0) {
    return(0)
  }

  n <- abs(n)
  m <- abs(m)
  if (m > n) {
    t <- n
    n <- m
    m <- t
  }
  while (m > 0) {
    t <- n
    n <- m
    m <- t %% m
  }
  return(n)
}

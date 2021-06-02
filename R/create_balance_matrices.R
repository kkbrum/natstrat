#' Create matrix of balance constraints for linear program
#'
#' This creates the matrix of constraints that seek covariate balance for use in \code{\link{balance_LP}()}
#' which creates the linear program used by \code{\link{optimize_controls}()} to choose which controls
#' to select.
#'
#' @inheritParams optimize_controls
#' @inheritParams balance_LP
#' @param return one of "all", "A", or "X", denoting whether all matrices should be
#' returned, or just A or just the X matrix blocks.
#'
#' @return A list containing up to three elements:
#' \describe{
#' \item{\code{A}}{The matrix of coefficients that will set the epsilons equal to the inbalances.}
#' \item{\code{x_blk}}{Just the coefficients equal to the inbalances for the main comparisons.}
#' \item{\code{x_blk2}}{Just the coefficients equal to the inbalances for the supplement.}
#' }
#'
#' @keywords internal
#' @import slam

# Set up balance constraints
create_balance_matrices <- function(X, z, N, nvars, kc2, q_s, q_star_s, return = "all") {
  X[is.na(X)] <- 0
  if (is.null(q_star_s)) {
    zero_blk <- simple_triplet_zero_matrix(nrow = nvars, ncol = N)
    zero_eps_blk <- simple_triplet_zero_matrix(nrow = nvars, ncol = 2 * nvars * kc2)
  } else {
    zero_blk <- simple_triplet_zero_matrix(nrow = nvars, ncol = 2 * N)
    zero_eps_blk <- simple_triplet_zero_matrix(nrow = nvars, ncol = 4 * nvars * kc2)
    zero_eps_blk2 <- simple_triplet_zero_matrix(nrow = nvars, ncol = 4 * nvars * kc2)
  }
  A <- NULL
  full_x_blk <- NULL
  full_x_blk2 <- NULL
  x_blk2 <- NULL
  pairs <- combn(unique(z), 2)
  for (pair_num in 1:kc2) {
    group1 <- pairs[1, pair_num]
    group2 <- pairs[2, pair_num]
    Q1 <- sum(q_s[levels(z) == group1, ])
    Q2 <- sum(q_s[levels(z) == group2, ])
    Q_star1 <- sum(q_star_s[levels(z) == group1, ])
    Q_star2 <- sum(q_star_s[levels(z) == group2, ])
    x_blk <- zero_blk
    eps_blk <- zero_eps_blk
    if (is.null(q_star_s) & Q1 > 0 & Q2 > 0) {
      x_blk[, which(z == group1)] <- t(X[z == group1, ] / Q1)
      x_blk[, which(z == group2)] <- -t(X[z == group2, ] / Q2)
    } else if (!is.null(q_star_s)) {
      x_blk2 <- zero_blk
      if (Q1 > 0 & Q2 > 0) {
        x_blk[, which(z == group1)] <- t(X[z == group1, ] / Q1)
        x_blk[, which(z == group2)] <- -t(X[z == group2, ] / Q2)
      }
      if (Q_star1 > 0 & Q_star2 > 0) {
        x_blk2[, (N + which(z == group1))] <- t(X[z == group1, ] / Q_star1)
        x_blk2[, (N + which(z == group2))] <- -t(X[z == group2, ] / Q_star2)
      }
    }
    eps_blk[, ((pair_num - 1) * nvars + 1):(pair_num * nvars)] <-
      simple_triplet_diag_matrix(rep(1, nvars))
    eps_blk[, (nvars * kc2 + (pair_num - 1) * nvars + 1):(nvars * kc2 + pair_num * nvars)] <-
      simple_triplet_diag_matrix(rep(-1, nvars))
    if (!is.null(q_star_s)) {
      eps_blk2 <- zero_eps_blk2
      eps_blk2[, (2 * nvars * kc2 + (pair_num - 1) * nvars + 1):(2 * nvars * kc2 + pair_num * nvars)] <-
        simple_triplet_diag_matrix(rep(1, nvars))
      eps_blk2[, (3 * nvars * kc2 + (pair_num - 1) * nvars + 1):(3 * nvars * kc2 + pair_num * nvars)] <-
        simple_triplet_diag_matrix(rep(-1, nvars))
    }
    new_A <- cbind(x_blk, eps_blk)
    A <- rbind(A, new_A)
    full_x_blk <- rbind(full_x_blk, x_blk)
    if (!is.null(q_star_s)) {
      new_A2 <- cbind(x_blk2, eps_blk2)
      A <- rbind(A, new_A2)
      full_x_blk2 <- rbind(full_x_blk2, x_blk2)
    }
  }
  if (return == "A") {
    return(list(A = A))
  }
  if (return == "X") {
    return(list(x_blk = full_x_blk, x_blk2 = full_x_blk2))
  } else {
    return(list(A = A, x_blk = full_x_blk, x_blk2 = full_x_blk2))
  }
}

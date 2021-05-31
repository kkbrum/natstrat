# Set up balance constraints
create_balance_matrices <- function(X, z, N, nvars, kc2, q_s, q_star_s) {
  X[is.na(X)] <- 0
  if (is.null(q_star_s)) {
    zero_blk <- matrix(0, nvars, N)
    zero_eps_blk <- matrix(0, nrow = nvars, ncol = 2 * nvars * kc2)
  } else {
    zero_blk <- matrix(0, nvars, 2 * N)
    zero_eps_blk <- matrix(0, nrow = nvars, ncol = 4 * nvars * kc2)
    zero_eps_blk2 <- matrix(0, nrow = nvars, ncol = 4 * nvars * kc2)
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
    x_blk <- zero_blk
    eps_blk <- zero_eps_blk
    Q_star1 <- sum(q_star_s[levels(z) == group1, ])
    Q_star2 <- sum(q_star_s[levels(z) == group2, ])
    if (is.null(q_star_s)) {
      x_blk[, z == group1] <- t(X[z == group1, ] / Q1)
      x_blk[, z == group2] <- -t(X[z == group2, ] / Q2)
    } else {
      x_blk2 <- zero_blk
      x_blk[, c(z == group1, rep(FALSE, N))] <- t(X[z == group1, ] / Q1)
      x_blk2[, c(rep(FALSE, N), z == group1)] <- t(X[z == group1, ] / Q_star1)
      x_blk[, c(z == group2, rep(FALSE, N))] <- -t(X[z == group2, ] / Q2)
      x_blk2[, c(rep(FALSE, N), z == group2)] <- -t(X[z == group2, ] / Q_star2)
    }
    eps_blk[, ((pair_num - 1) * nvars + 1):(pair_num * nvars)] <-
      diag(1, nvars, nvars)
    eps_blk[, (nvars * kc2 + (pair_num - 1) * nvars + 1):(nvars * kc2 + pair_num * nvars)] <-
      diag(-1, nvars, nvars)
    if (!is.null(q_star_s)) {
      eps_blk2 <- zero_eps_blk2
      eps_blk2[, (2 * nvars * kc2 + (pair_num - 1) * nvars + 1):(2 * nvars * kc2 + pair_num * nvars)] <-
        diag(1, nvars, nvars)
      eps_blk2[, (3 * nvars * kc2 + (pair_num - 1) * nvars + 1):(3 * nvars * kc2 + pair_num * nvars)] <-
        diag(-1, nvars, nvars)
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
  return(list(A = A, x_blk = full_x_blk, x_blk2 = full_x_blk2))
}

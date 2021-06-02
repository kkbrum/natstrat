#' Linear program that selects which controls to use in order to optimize balance
#'
#' This linear program is used by \code{\link{optimize_controls}()} to choose which controls
#' to select.
#'
#' @inheritParams optimize_controls
#' @param q_s a named vector or matrix indicating how many control units are to be selected from each stratum.
#'   If there is one control group and all treated units are desired, this can be a vector; otherwise,
#'   this should have one row per treatment group, where the order of the rows matches the order of
#'   the levels of \code{z}, including the treated level.
#' @param st_vals the unique stratum levels contained in \code{st}.
#' @param S the number of unique stratum levels contained in \code{st}.
#' @param N the total number of available controls in the data.
#'
#' @return A list containing two elements:
#' \describe{
#' \item{\code{lpdetails}}{The output of either \code{gurobi()} or \code{\link[Rglpk]{Rglpk_solve_LP}()},
#' except that if \code{gurobi()} is used, the elements \code{objval} and \code{x}
#' are renamed \code{optimum} and \code{solution}
#' to be consistent with the output of \code{\link[Rglpk]{Rglpk_solve_LP}()}.}
#' \item{\code{o}}{The original output of either \code{gurobi()} or \code{\link[Rglpk]{Rglpk_solve_LP}()}.}
#' }
#'
#' @keywords internal
#' @import ramify
#' @import slam

balance_LP <- function(z, X, importances, st, st_vals, S, q_s, N,
                       solver, integer, time_limit, q_star_s = NULL, weight_star = 1) {

  if (solver == "gurobi" && !requireNamespace("gurobi", quietly = TRUE)) {
    stop("Package \'gurobi\' needed if \"solver\" parameter set to \"gurobi\". Please
         install it or switch the \"solver\" parameter to \"Rglpk\".",
         call. = FALSE)
  }

  k <- length(unique(z))
  kc2 <- choose(k, 2)

  # Set up and solve the linear program
  model <- list()
  params <- list(TimeLimit = time_limit, OutputFlag = 0)

  nvars <- dim(X)[2]  # number of variables
  if (is.null(q_star_s)) {
    model$obj <- c(rep(0, N), rep(rep(importances, 2), kc2))
  } else {
    model$obj <- c(rep(0, 2 * N), rep(rep(importances, 2), kc2), rep(rep(importances * weight_star, 2), kc2))
  }

  model$A <- create_balance_matrices(X = X, z = z, N = N, nvars = nvars,
                          kc2 = kc2, q_s = q_s, q_star_s = q_star_s, return = "A")$A

  groups <- unique(z)
  # Now, append stratum size constraints for comparison 1
  st_mats <- simple_triplet_zero_matrix(nrow = k * S, ncol = N)
  for (group_num in 1:length(groups)) {
    group <- groups[group_num]
    st_mats[((group_num - 1) * S + 1):(group_num * S), which(z == group)] <- 1 * outer(st_vals, st[z == group], "==")
  }
  if (is.null(q_star_s)) {
    model$A <- rbind(model$A,
                     cbind(st_mats, simple_triplet_zero_matrix(nrow = k * S, ncol = 2 * kc2 * nvars)))
  } else {
    model$A <- rbind(model$A,
                     cbind(st_mats, simple_triplet_zero_matrix(nrow = k * S, ncol = N + 4 * kc2 * nvars)))
    # Also append stratum size constraints for comparison 2 if exists
    model$A <- rbind(model$A,
                     cbind(simple_triplet_zero_matrix(nrow = k * S, ncol = N),
                           st_mats, simple_triplet_zero_matrix(nrow = k * S, ncol = 4 * kc2 * nvars)))
  }

  # Now, if two comparisons, add constraint that two as for a unit add to <= 1
  # (so that one unit is not chosen for both comparisons)
  if (!is.null(q_star_s)) {
    mat<- cbind(simple_triplet_diag_matrix(rep(1, N)), simple_triplet_diag_matrix(rep(1, N)))
    model$A <- rbind(model$A, cbind(mat, simple_triplet_zero_matrix(nrow = N, ncol = 4 * kc2 * nvars)))
  }

  # Constraints for eps are equalities, number of controls per strata are equalities
  # Constraints for units only counting in one comparison are <=
  if (is.null(q_star_s)) {
    model$sense <- c(rep("==", kc2 * nvars), rep("==", k * S))
  } else {
    model$sense <- c(rep("==", 2 * kc2 * nvars), rep("==", 2 * k * S), rep("<=", N))
  }

  # right hand side of constraints
  if (is.null(q_star_s)) {
    model$rhs <- c(rep(0, kc2 * nvars), ramify::flatten(q_s))
  } else {
    model$rhs <- c(rep(0, 2 * kc2 * nvars),
                   ramify::flatten(q_s), ramify::flatten(q_star_s),
                   rep(1, N))
  }

  if (is.null(q_star_s)) {
    ndecv <- as.integer(N + (2 * kc2 * nvars))  # number of decision variables
    model$ub <- c(rep(1, N), rep(Inf, 2 * kc2 * nvars))
  } else {
    ndecv <- as.integer(2 * N + (4 * kc2 * nvars))  # number of decision variables
    model$ub <- c(rep(1, 2 * N), rep(Inf, 4 * kc2 * nvars))
  }
  model$lb <- rep(0, ndecv)
  bounds <- list(lower = list(ind = 1:ndecv, val = model$lb),
                 upper = list(ind = 1:ndecv, val = model$ub))
  if (integer) {
    if (is.null(q_star_s)) {
      model$vtype <- c(rep("B", N), rep("C", 2 * kc2 * nvars))
    } else {
      model$vtype <- c(rep("B", 2 * N), rep("C", 4 * kc2 * nvars))
    }
  } else {
    model$vtype <- rep("C", ndecv)
  }

  if (solver == "Rglpk") {
    if (params$TimeLimit < Inf) {
      params$TimeLimit <- params$TimeLimit * 1000
    } else {
      params$TimeLimit <- 0
    }
    o <- Rglpk::Rglpk_solve_LP(model$obj, model$A, model$sense, model$rhs, bounds = bounds,
                               types = model$vtype, control = list(
                                 canonicalize_status = FALSE, tm_limit = params$TimeLimit))
    if (o$status != 5) {
      warning("No optimal solution found for the linear program.")
      return(NULL)
    }
    lpdetails <- o
  }
  if (solver == "gurobi") {
    # Note that for gurobi, all inequalities are interpreted to be "or equal to"
    if (is.null(q_star_s)) {
      model$sense <- c(rep("=", kc2 * nvars), rep("=", k * S))
    } else {
      model$sense <- c(rep("=", 2 * kc2 * nvars), rep("=", 2 * k * S),
                       rep("<", N))
    }
    o <- gurobi::gurobi(model, params)
    if (o$status != "OPTIMAL") {
      warning("No solution found for the linear program.")
      return(NULL)
    }
    lpdetails <- o
    names(o)[c(6,7)] <- c("optimum", "solution")
  }

  return(list(lpdetails = lpdetails, o = o))
}

#' Linear program that selects which controls to use in order to optimize balance
#'
#' This linear program is used by \code{\link{optimize_controls}()} to choose which controls
#' to select.
#'
#' @inheritParams optimize_controls
#' @param q_s a named vector indicating how many control units are to be selected from each stratum.
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

balance_LP <- function(z, X, importances, st, st_vals, S, q_s, N,
                       solver, integer, time_limit) {
  if (solver == "gurobi" && !requireNamespace("gurobi", quietly = TRUE)) {
    stop("Package \'gurobi\' needed if \"solver\" parameter set to \"gurobi\". Please
         install it or switch the \"solver\" parameter to \"Rglpk\".",
         call. = FALSE)
  }

  # Set up and solve the linear program
  model <- list()
  params <- list(TimeLimit = time_limit, OutputFlag = 0)

  nvars <- dim(X)[2]  # number of variables
  X[is.na(X)] <- 0
  X0 <- X[z == 0, ]
  model$obj <- c(rep(0, N), rep(importances, 2))
  blk1 <- t(X0)
  ident <- diag(1, nvars, nvars)  # identity matrix
  model$A <- cbind(blk1 / sum(q_s), ident, -ident)  # constraints, individual vars

  # Now, append stratum size constraints
  model$A <- rbind(model$A, cbind(1 * outer(st_vals, st[z == 0], "=="), matrix(0, S,
                                                                               2 * nvars)))
  # Constraints for eps are equalities, number of controls per strata are equalities
  model$sense <- c(rep("==", nvars), rep("==", S))
  model$rhs <- c(rep(0, nvars), q_s)  # Right hand side of constraints

  ndecv <- as.integer(N + (2 * nvars))  # number of decision variables
  model$lb <- rep(0, ndecv)
  model$ub <- c(rep(1, N), rep(Inf, 2 * nvars))
  bounds <- list(lower = list(ind = 1:ndecv, val = model$lb),
                 upper = list(ind = 1:ndecv, val = model$ub))
  if (integer) {
    model$vtype <- c(rep("B", N), rep("C", 2 * nvars))
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
    model$sense <- c(rep("=", nvars), rep("=", S))
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

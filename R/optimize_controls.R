#' Select control units that optimize covariate balance
#'
#' Select control units within strata that optimize covariate balance.
#' Uses randomized rounding of a linear program or a mixed
#' integer linear program.
#'
#' @inheritParams stand
#' @param st a stratum vector with the \code{i}th entry equal to the
#'   stratum of unit \code{i}. This should have the same order of units and length
#'   as \code{z}.
#' @param X a matrix or data frame containing constraints in the columns. The number
#'   of rows should equal the length of \code{z}. Balance is achieved when a constraint
#'   sums to 0, such that numbers closer to 0 are better. When a constraint
#'   does not apply to a particular unit, the entry should be \code{NA}.
#'   This should typically be generated using \code{\link{generate_constraints}()}.
#' @param treated_star which treatment value should be considered the treated units
#' for the supplemental comparison. This
#' must be one of the values of \code{z}.
#' @param ratio a numeric or vector specifying the desired ratio of controls to `treated` in
#'   each stratum. If there is one control group and all treated units should be included,
#'   this can be a numeric. Otherwise, this should be
#'   a vector with one entry per treatment group, in the same order as the levels of
#'   \code{z}, including the treated level. If \code{NULL}, \code{q_s} should be specified.
#' @param ratio_star a numeric or vector specifying the desired ratio of supplemental units to `treated_star` in
#'   each stratum. This should be
#'   a vector with one entry per treatment group, in the same order as the levels of
#'   \code{z}, including the treated level. If \code{NULL}, \code{q_star_s} should be specified.
#' @param q_s a named vector or matrix indicating how many units are to be selected from each stratum.
#'   If there is one control group and all treated units are desired, this can be a vector; otherwise,
#'   this should have one row per treatment group, where the order of the rows matches the order of
#'   the levels of \code{z}, including the treated level.
#'   If \code{NULL}, \code{ratio} should be specified. If both are specified, \code{q_s} will take priority.
#'   Typically, if the desired ratio is not feasible for every stratum, \code{q_s} should be generated
#'   using \code{\link{generate_qs}()}.
#' @param q_star_s a named vector or matrix indicating how many supplemental units are to be selected from each stratum.
#'   This should have one row per treatment group, where the order of the rows matches the order of
#'   the levels of \code{z}, including the treated level.
#'   If \code{NULL}, \code{ratio_star} should be specified. If both are specified, \code{q_star_s} will take priority.
#'   Typically, if the desired ratio is not feasible for every stratum, \code{q_star_s} should be generated
#'   using \code{\link{generate_qs}()}.
#' @param weight_star a numeric stating how much to prioritize balance between the supplemental units as
#' compared to balance between the main units.
#' @param importances a vector with length equal to the number of constraints or columns
#'   in \code{X}. This can be generated using \code{\link{generate_constraints}()} and each nonnegative value
#'   denotes how much to prioritize each constraint, with the default being 1
#'   for all constraints.
#' @param integer a logical stating whether to use a mixed integer programming solver
#'   instead of randomized rounding. Default is \code{FALSE}.
#' @param solver a character stating which solver to use to run the linear program.
#'   Options are "Rglpk" (default) or "gurobi". You must have the 'gurobi' package
#'   installed to use the "gurobi" option. If available, this is the recommended solver.
#' @param seed the seed to use when doing the randomized rounding of the linear program.
#'   This will allow results to be reproduced if desired. The default is \code{NULL},
#'   which will choose a random seed to use and return.
#' @param runs the number of times to run randomized rounding of the linear solution.
#'   The objective values of all runs will be reported, but the detailed results
#'   will only be reported for the run with the lowest objective value. The default is 1.
#' @param time_limit numeric stating maximum amount of seconds for which the
#'   program is allowed to run before aborting. Default is \code{Inf} for no time limit.
#' @param correct_sizes boolean stating whether the desired sample sizes should
#'  be exactly correct (if \code{correct_sizes = TRUE}) or only need to be correct
#'  in expectation. For nested comparisons, sample sizes may only be
#'  correct in expectation.
#' @param low_memory boolean stating whether some outputs should not be included
#'  due to the scale of the problem being too large compared to memory space.
#'  If \code{TRUE}, \code{eps} and \code{eps_star} will not be reported. Inbalances
#'  can be computed post hoc using the \code{\link{check_balance}()} instead.
#' @param threads The maximum number of threads that should be used. This is only
#'  applicable if \code{solver = 'gurobi'}.
#'
#' @return List containing:
#' \describe{
#'   \item{\code{objective}}{objective value of the randomized rounding or mixed integer
#'     linear program solution.}
#'   \item{\code{objective_wo_importances}}{objective value of the randomized rounding or mixed integer
#'     linear program solution not weighted by the variable importances.}
#'   \item{\code{eps}}{the amount of imbalance obtained in each constraint from the linear program.
#'   The row names specify the covariate, the population of interest, and, if there are
#'   more than two comparison groups, which groups are being compared.}
#'   \item{\code{eps_star}}{same as \code{eps} but for the supplemental units instead of the units
#'   in the main comparison.}
#'   \item{\code{importances}}{the importance of each on the balance constraints.}
#'   \item{\code{selected}}{whether each unit was selected for the main comparison.}
#'   \item{\code{selected_star}}{whether each unit was selected for the supplement.}
#'   \item{\code{pr}}{the linear program weight assigned to each unit for the main comparison.}
#'   \item{\code{pr_star}}{the linear program weight assigned to each unit for the supplement.}
#'   \item{\code{rrdetails}}{A list containing:
#'   \describe{
#'   \item{\code{seed}}{the seed used before commencing the random sampling.}
#'   \item{\code{run_objectives}}{the objective values for each run of randomized rounding.}
#'   \item{\code{run_objectives_wo_importances}}{the objective values for each run of randomized rounding,
#'   not scaled by constraint importances.}
#'   }}
#'   \item{\code{lpdetails}}{the full return of the function \code{\link[Rglpk]{Rglpk_solve_LP}()}
#'     or \code{gurobi()} plus information about the epsilons and objective values
#'     for the linear program solution.}
#' }
#'
#' @importFrom stats complete.cases
#' @export
#'
#' @examples
#'
#' data('nh0506')
#'
#' # Create strata
#' age_cat <- cut(nh0506$age,
#'                breaks = c(19, 39, 50, 85),
#'                labels = c('< 40 years', '40 - 50 years', '> 50 years'))
#' strata <- age_cat : nh0506$sex
#'
#' # Balance age, race, education, poverty ratio, and bmi both across and within the levels of strata
#' constraints <- generate_constraints(
#'                  balance_formulas = list(age + race + education + povertyr + bmi ~ 1 + strata),
#'                  z = nh0506$z,
#'                  data = nh0506)
#'
#' # Choose one control for every treated unit in each stratum,
#' # balancing the covariates as described by the constraints
#' results <- optimize_controls(z = nh0506$z,
#'                              X = constraints$X,
#'                              st = strata,
#'                              importances = constraints$importances,
#'                              ratio = 1)
#'
#' # If you want to use a ratio that's not feasible,
#' # you can supply a vector of the desired number of controls per stratum, q_s,
#' # typically generated by creating a distance matrix between strata and using
#' # generate_qs():
#'
#' age_dist <- matrix(data = c(0, 1, 2, 1, 0, 1, 2, 1, 0),
#'                    nrow = 3,
#'                    byrow = TRUE,
#'                    dimnames = list(levels(age_cat), levels(age_cat)))
#'
#' sex_dist <- matrix(data = c(0, 1, 1, 0),
#'                    nrow = 2,
#'                    dimnames = list(levels(nh0506$sex), levels(nh0506$sex)))
#'
#' strata_dist <- create_dist_matrix(age_dist, sex_dist)
#'
#' qs <- generate_qs(z = nh0506$z,
#'                   st = strata,
#'                   ratio = 2.5,
#'                   max_ratio = 2.6,
#'                   max_extra_s = 0,
#'                   strata_dist = strata_dist)
#'
#' results <- optimize_controls(z = nh0506$z,
#'                              X = constraints$X,
#'                              st = strata,
#'                              importances = constraints$importances,
#'                              q_s = qs)

# TODO: Use ratio_star input if provided instead of q_star_s
# TODO: Process q_star_s if not matrix

optimize_controls <- function(z, X, st, importances = NULL, treated = 1,
                              ratio = NULL, q_s = NULL, treated_star = NULL,
                              ratio_star = NULL, q_star_s = NULL, weight_star = 1,
                              integer = FALSE, solver = "Rglpk",
                              seed = NULL, runs = 1,
                              time_limit = Inf, threads = 1, correct_sizes = TRUE,
                              low_memory = FALSE) {

  # Make sure inputs are good
  verify_inputs(X = X, importances = importances, ratio = ratio, q_s = q_s,
                st = st, z = z, treated = treated, integer = integer, solver = solver)
  multi_comp <- !is.null(q_star_s)
  z <- factor(z)
  group <- levels(z)
  k <- length(group)
  kc2 <- choose(k, 2)
  if (!(!multi_comp & is.null(ratio_star))) {
    # TODO: Add verify inputs here to check everything about these
    if (is.null(treated_star)) {
      stop("If `q_star_s` or `ratio_star` are not `NULL`, `treated_star` must also specify the treatment group for reference in the second comparison.")
    }
    if (correct_sizes) {
      correct_sizes <- FALSE
      warning("Sample sizes are only correct in expectation for multiple comparisons. `correct_sizes` has thus been switched to `FALSE`.")
    }
  }
  # Look at strata counts
  frtab <- table(z, st)
  stratios <- frtab / frtab[group == treated, ]
  st_vals <- as.character(colnames(frtab))  # stratum values
  S <- length(st_vals)  # number of strata

  # TODO: Move all this figuring to a new function that outputs the final q_s
  if (!is.null(q_s) & is.vector(q_s)) {
    q_s <- matrix(rep(q_s, k), byrow = TRUE, nrow = k, dimnames = list(NULL, names(q_s)))
    q_s[group == treated, ] <- frtab[group == treated, ]
  }
  if (!is.null(ratio) & length(ratio) == 1) {
    ratio <- rep(ratio, k)
    ratio[group == treated] <- 1
  }
  # Determine number of controls desired for each stratum
  if (is.null(q_s)) {
    if (is.null(ratio)) {
      ratio <- sapply(1:nrow(stratios), function(i) min(1, min(stratios[i, ])))
    }
    q_s <- round(ratio %*% t(frtab[group == treated, ]))
    n_s <- table(z, st)
    if (any(q_s > n_s)) {
      stop("The ratio you specified is not feasible.
            Please supply `q_s` instead of `ratio` or lower the `ratio` input.",
           call. = FALSE)
    }
  } else {
    q_s <- q_s[, st_vals]
  }

  # Prepare importances
  if (is.null(importances)) {
    importances <- setNames(rep(1, ncol(X)), colnames(X))
  } else {
    importances <- importances[colnames(X)]
  }

  nvars_per_group <- dim(X)[2]
  nvars <- nvars_per_group * kc2
  N <- length(z)

  # Run linear program to choose control units
  lp_results <- balance_LP(z = z, X = X, importances = importances,
                           st = st, st_vals = st_vals, S = S,
                           q_s = q_s, q_star_s = q_star_s, weight_star = weight_star,
                           N = N, integer = integer, solver = solver,
                           time_limit = time_limit, threads = threads)

  if (is.null(lp_results)) {
    return(NULL)
  } else {
    Q <- rowSums(q_s)

    if (!low_memory) {
      # Epsilons for the linear program solution
      if (!multi_comp) {
        lp_results$lpdetails$eps <- lp_results$o$solution[(N + 1):(N + (2 * nvars))]
        lp_results$lpdetails$eps <- matrix(lp_results$lpdetails$eps, nvars, 2)
      } else {
        lp_results$lpdetails$eps <- lp_results$o$solution[(2 * N + 1):(2 * N + (2 * nvars))]
        lp_results$lpdetails$eps <- matrix(lp_results$lpdetails$eps, nvars, 2)
        lp_results$lpdetails$eps_star <- lp_results$o$solution[(2 * N + 2 * nvars + 1):(2 * N + (4 * nvars))]
        lp_results$lpdetails$eps_star <- matrix(lp_results$lpdetails$eps_star, nvars, 2)
      }
      if (!is.null(colnames(X))) {
        if (k == 2) {
          rownames(lp_results$lpdetails$eps) <- colnames(X)
          if (multi_comp) {
            rownames(lp_results$lpdetails$eps_star) <- colnames(X)
          }
        } else {
          rownames(lp_results$lpdetails$eps) <- 1:nvars
          if (multi_comp) {
            rownames(lp_results$lpdetails$eps_star) <- 1:nvars
          }
          pairs <- combn(group, 2)
          for (pair_num in 1:kc2) {
            group1 <- pairs[1, pair_num]
            group2 <- pairs[2, pair_num]
            rownames(lp_results$lpdetails$eps)[((pair_num - 1) * nvars_per_group + 1):(pair_num * nvars_per_group)] <-
              paste0(colnames(X), "_", group1, ":", group2)
            if (multi_comp) {
              rownames(lp_results$lpdetails$eps_star)[((pair_num - 1) * nvars_per_group + 1):(pair_num * nvars_per_group)] <-
                paste0(colnames(X), "_", group1, ":", group2)
            }
          }
        }
      }

      colnames(lp_results$lpdetails$eps) <- c("positive", "negative")
      if (multi_comp) {
        colnames(lp_results$lpdetails$eps_star) <- c("positive", "negative")
      }
    }

    # Record objectives
    lp_results$lpdetails$objective <- lp_results$o$optimum
    # Obj without importances for LP
    if (!multi_comp) {
      # sum of epsilons
      lp_results$lpdetails$objective_wo_importances <- sum(lp_results$o$solution[(N + 1):(N + (2 * nvars))])
    } else {
      # sum of epsilons and epsilon_stars
      lp_results$lpdetails$objective_wo_importances <- sum(lp_results$o$solution[(2 * N + 1):(2 * N + (4 * nvars))])
    }

    best_objective <- Inf
    run_objectives <- rep(NA, runs)
    run_objectives_wo_importances <- rep(NA, runs)

    if (is.null(seed)) {
      seed <- sample(1:1000000, 1)
    }
    set.seed(seed)

    balance_matrices <- create_balance_matrices(X = X, z = z, N = N, nvars = nvars_per_group,
                                                kc2 = kc2, q_s = q_s, q_star_s = q_star_s, return = "X")

    eps <- NULL
    eps_star <- NULL

    for (run in 1:runs) {
      # Run randomized rounding
      if (correct_sizes) {
        rr_results_temp <- randomized_rounding(o = lp_results$o, N = N, st = st,
                                               st_vals = st_vals, S = S, z = z)
      } else {
        rr_results_temp <- randomized_rounding_expectation(o = lp_results$o, N = N,
                                                           multi_comp = multi_comp)
      }

      # Calculate and format results

      if (!low_memory) {
        # Epsilons for the randomized rounding (or integer) solution
        eps_temp <- matrix(0, nvars, 2)
        eps_temp_star <- matrix(0, nvars, 2)
        for (i in 1:nvars) {
          row <- as.vector(balance_matrices$x_blk[i, ])
          inbalance <- sum(row * rr_results_temp$select)
          if (inbalance < 0) {
            eps_temp[i, 1] <- abs(inbalance)
          } else if (inbalance > 0) {
            eps_temp[i, 2] <- inbalance
          }
          if (multi_comp) {
            row <- as.vector(balance_matrices$x_blk2[i, ])
            inbalance_star <- sum(row * rr_results_temp$select)
            if (inbalance_star < 0) {
              eps_temp_star[i, 1] <- abs(inbalance_star)
            } else if (inbalance > 0) {
              eps_temp_star[i, 2] <- inbalance_star
            }
          }
        }
        if (k == 2) {
          if (!is.null(colnames(X))) {
            rownames(eps_temp) <- colnames(X)
            if (multi_comp) {
              rownames(eps_temp_star) <- colnames(X)
            }
          }
        } else {
          if (!is.null(colnames(X))) {
            rownames(eps_temp) <- 1:nvars
            if (multi_comp) {
              rownames(eps_temp_star) <- 1:nvars
            }
            pairs <- combn(unique(z), 2)
            for (pair_num in 1:kc2) {
              group1 <- pairs[1, pair_num]
              group2 <- pairs[2, pair_num]
              rownames(eps_temp)[((pair_num - 1) * nvars_per_group + 1):(pair_num * nvars_per_group)] <- paste0(colnames(X), "_", group1, ":", group2)
              if (multi_comp) {
                rownames(eps_temp_star)[((pair_num - 1) * nvars_per_group + 1):(pair_num * nvars_per_group)] <- paste0(colnames(X), "_", group1, ":", group2)
              }
            }
          }
        }

        colnames(eps_temp) <- c("positive", "negative")
        if (multi_comp) {
          colnames(eps_temp_star) <- c("positive", "negative")
        }
      }

      # Objective value for the randomized rounding (or integer) solution
      run_objectives[run] <- sum(abs(importances * (matrix(balance_matrices$x_blk,
                                                           nrow = nrow(balance_matrices$x_blk),
                                                           ncol = ncol(balance_matrices$x_blk))
                                                    %*% rr_results_temp$select)))
      run_objectives_wo_importances[run] <- sum(abs(matrix(balance_matrices$x_blk,
                                                           nrow = nrow(balance_matrices$x_blk),
                                                           ncol = ncol(balance_matrices$x_blk))
                                                    %*% rr_results_temp$select))
      if (multi_comp) {
        run_objectives[run] <- run_objectives[run] +
          sum(abs(importances * weight_star * (matrix(balance_matrices$x_blk2,
                                                      nrow = nrow(balance_matrices$x_blk2),
                                                      ncol = ncol(balance_matrices$x_blk2))
                                               %*% rr_results_temp$select)))
        run_objectives_wo_importances[run] <- run_objectives_wo_importances[run] +
          sum(abs(matrix(balance_matrices$x_blk2,
                         nrow = nrow(balance_matrices$x_blk2),
                         ncol = ncol(balance_matrices$x_blk2)) %*% rr_results_temp$select))
      }

      if (run_objectives[run] < best_objective) {
        if (!low_memory) {
          eps <- eps_temp
          eps_star <- eps_temp_star
        }
        objective_wo_importances <- run_objectives_wo_importances[run]
        rr_results <- rr_results_temp
        best_objective <- run_objectives[run]
      }

    }

    selected_star <- NULL
    pr_star <- NULL
    if (multi_comp) {
      selected_star <- rr_results$select[(N+1):(2*N)]
      pr_star <- rr_results$pr[(N+1):(2*N)]
    }

    return(list(objective = best_objective,
                objective_wo_importances = objective_wo_importances,
                eps = eps,
                eps_star = eps_star,
                importances = importances,
                weight_star = weight_star,
                selected = rr_results$select[1:N],
                pr = rr_results$pr[1:N],
                selected_star = selected_star,
                pr_star = pr_star,
                rrdetails = list(
                  seed = seed,
                  run_objectives = run_objectives,
                  run_objectives_wo_importances = run_objectives_wo_importances),
                lpdetails = lp_results$lpdetails))
  }
}

#' Verify the inputs to \code{\link{optimize_controls}()}
#'
#' Makes sure that the inputs to \code{\link{optimize_controls}()} are in the correct
#' format and feasible.
#'
#' @inheritParams optimize_controls
#'
#' @return No return value. If there is a problem with the inputs to \code{\link{optimize_controls}()},
#' an error is raised.
#'
#' @keywords internal

verify_inputs <- function(X, importances, ratio, q_s, st, z, treated, integer, solver) {
  if (is.data.frame(X))
    X <- as.matrix(X)
  stopifnot(!is.null(ratio) || !is.null(q_s))
  stopifnot(is.null(ratio) | (is.vector(ratio) & all(ratio > 0)))
  stopifnot(is.vector(st) | is.factor(st))
  stopifnot(is.vector(z) | is.factor(z))
  stopifnot(is.matrix(X))
  stopifnot(length(z) == (dim(X)[1]))
  stopifnot(length(z) == length(st))
  stopifnot(is.logical(integer))
  if (!solver %in% c("gurobi", "Rglpk")) {
    stop("\"solver\" must be one of \"gurobi\" or \"Rglpk\".",
         call. = FALSE)
  }
  if (solver == "gurobi" && !requireNamespace("gurobi", quietly = TRUE)) {
    stop("Package \'gurobi\' needed if \"solver\" parameter set to \"gurobi\". Please
         install it or switch the \"solver\" parameter to \"Rglpk\".",
         call. = FALSE)
  }
  z <- factor(z)
  group <- levels(z)
  if (!treated %in% group) {
    stop("\"treated\" must be one of the values in \"z\".")
  }
  frtab <- table(z, st)
  if (min(frtab[group == treated, ]) == 0) {
    warning("Note that at least one stratum has no treated individuals.")
  }
  if (!is.null(q_s)) {
    if (is.vector(q_s)) {
      q_s <- matrix(c(q_s, frtab[group == treated, ]), byrow = TRUE, nrow = 2, dimnames = list(NULL, names(q_s)))
      if (group[1] == treated) {
        q_s <- q_s[c(2, 1), ]
      }
    }
    n_s <- table(z, st)
    if (any(q_s[, colnames(n_s)] > n_s)) {
      stop("At least one of the entries for `q_s` is greater than the number of units available in the stratum.
           Please lower `q_s` such that all entries are at most the number of available units.",
           call. = FALSE)
    }
  }
  stopifnot(sort(names(importances)) == sort(colnames(X)))
}


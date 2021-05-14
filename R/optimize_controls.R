#' Select control units that optimize covariate balance
#'
#' Select control units within strata that optimize covariate balance.
#' Uses randomized rounding of a linear program or a mixed
#' integer linear program.
#'
#' @inheritParams stand
#' @param X a matrix or data frame containing constraints in the columns. The number
#'   of rows should equal the length of \code{z}. Balance is achieved when a constraint
#'   sums to 0, such that numbers closer to 0 are better. When a constraint
#'   does not apply to a particular unit, the entry should be \code{NA}.
#'   This should typically be generated using \code{\link{generate_constraints}}.
#' @param ratio a numeric specifying the desired ratio of controls to treated in
#'   each stratum. If \code{NULL}, \code{q_s} should be specified.
#' @param q_s a named vector indicating how many control units are to be selected from each stratum.
#'   If \code{NULL}, \code{ratio} should be specified. If both are specified, \code{q_s} will take priority.
#'   Typically, if the desired ratio is not feasible for every stratum, \code{q_s} should be generated
#'   using \code{\link{generate_qs}}.
#' @param importances a vector with length equal to the number of constraints or columns
#'   in \code{X}. This can be generated using \code{\link{generate_constraints}} and each nonnegative value
#'   denotes how much to prioritize each constraint, with the default being 1
#'   for all constraints.
#' @param integer a logical stating whether to use a mixed integer programming solver
#'   instead of randomized rounding. Default is \code{FALSE}.
#' @param solver a character stating which solver to use to run the linear program.
#'   Options are "Rglpk" (default) or "gurobi". You must have gurobi installed to use "gurobi".
#'   If gurobi is available, this is the recommended solver.
#' @param seed the seed to use when doing the randomized rounding of the linear program.
#'   This will allow results to be reproduced if desired. The default is \code{NULL},
#'   which will choose a random seed to use and return.
#' @param runs the number of times to run randomized rounding of the linear solution.
#'   The objective values of all runs will be reported, but the detailed results
#'   will only be reported for the run with the lowest objective value. The default is 1.
#' @param time_limit numeric stating maximum amount of seconds for which the
#'   program is allowed to run before aborting. Default is \code{Inf} for no time limit.
#'
#' @return List containing:
#' \describe{
#'   \item{\code{objective}}{objective value of the randomized rounding or mixed integer
#'     linear program solution.}
#'   \item{\code{objective_wo_importances}}{objective value of the randomized rounding or mixed integer
#'     linear program solution not weighted by the variable importances.}
#'   \item{\code{eps}}{the amount of imbalance obtained in each constraint from the linear program.}
#'   \item{\code{importances}}{the importance of each on the balance constraints.}
#'   \item{\code{selected}}{the selected treated and control units.}
#'   \item{\code{controls}}{the linear program weight assigned to each control and
#'     whether it was selected by randomized rounding.}
#'   \item{\code{rrdetails}}{A list containing:
#'   \describe{
#'   \item{\code{seed}}{the seed used before commencing the random sampling.}
#'   \item{\code{raw_objective}}{objective value of the randomized rounding or mixed integer
#'     linear program solution before the denominator has been corrected for the number of
#'     units chosen with missing covariate values.}
#'   \item{\code{raw_objective_wo_importances}}{objective value of the randomized rounding or mixed integer
#'     linear program solution not weighted by the variable importances
#'     before the denominator has been corrected for the number of
#'     units chosen with missing covariate values.}
#'   \item{\code{raw_eps}}{the amount of imbalance obtained in each constraint from the linear program,
#'     before the denominators have been corrected for the number of
#'     units chosen with missing covariate values.}
#'   \item{\code{run_raw_objectives}}{the objective values for each run of randomized rounding,
#'   before denominators have been corrected for missingness.}
#'   \item{\code{run_raw_objectives_wo_importances}}{the objective values for each run of randomized rounding,
#'   before denominators have been corrected for missingness, not scaled by constraint importances.}
#'   \item{\code{run_objectives}}{the objective values for each run of randomized rounding.}
#'   \item{\code{run_objectives_wo_importances}}{the objective values for each run of randomized rounding,
#'   not scaled by constraint importances.}
#'   }}
#'   \item{\code{lpdetails}}{the full return of the function \code{\link[Rglpk]{Rglpk_solve_LP}}
#'     or gurobi plus information about the epsilons and objective values
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

optimize_controls <- function(z, X, st, importances = NULL,
                              ratio = NULL, q_s = NULL,
                              integer = FALSE, solver = "Rglpk",
                              seed = NULL, runs = 1,
                              time_limit = Inf) {

  # Make sure inputs are good
  verify_inputs(X = X, importances = importances, ratio = ratio, q_s = q_s, st = st, z = z, integer = integer, solver = solver)

  # Look at strata counts
  frtab <- table(z, st)
  stratios <- frtab[1, ]/frtab[2, ]
  st_vals <- as.character(colnames(frtab))  # stratum values
  S <- length(st_vals)  # number of strata

  # Determine number of controls desired for each stratum
  if (is.null(q_s)) {
    if (is.null(ratio)) {
      ratio <- min(1, min(stratios))
    }
    q_s <- round(ratio * frtab[2, ])
  } else {
    q_s <- q_s[st_vals]
  }

  # Prepare importances
  if (is.null(importances)) {
    importances <- setNames(rep(1, ncol(X)), colnames(X))
  } else {
    importances <- importances[colnames(X)]
  }

  nvars <- dim(X)[2]
  N <- sum(1 - z)
  # Run linear program to choose control units
  lp_results <- balance_LP(z = z, X = X, importances = importances, st = st, st_vals = st_vals, S = S,
                           q_s = q_s, N = N, integer = integer, solver = solver,
                           time_limit = time_limit)

  if (is.null(lp_results)) {
    return(NULL)
  } else {
    Q <- sum(q_s)

    # Epsilons for the linear program solution
    lp_results$lpdetails$raw_eps <- lp_results$o$solution[(N + 1):(N + (2 * nvars))]
    lp_results$lpdetails$raw_eps <- matrix(lp_results$lpdetails$raw_eps, nvars, 2)
    if (!is.null(colnames(X))) {
      rownames(lp_results$lpdetails$raw_eps) <- colnames(X)
    }
    colnames(lp_results$lpdetails$raw_eps) <- c("positive", "negative")

    # Record raw objectives
    lp_results$lpdetails$raw_objective <- lp_results$o$optimum

    # Obj without importances for LP
    lp_results$lpdetails$raw_objective_wo_importances <- sum(lp_results$lpdetails$raw_eps)

    # Corrected epsilons for LP

    # Locate missing values in each covariate (within any stratum)
    stripped_covs <- sapply(strsplit(row.names(lp_results$lpdetails$raw_eps), "_"),
                            function(k) {paste0(k[-length(k)], collapse = "_")})
    unique_covs <- unique(stripped_covs)
    missing_instances <- data.frame(matrix(NA, nrow = N, ncol = length(unique_covs)))
    names(missing_instances) <- unique_covs
    for (cov in unique_covs) {
      missing_instances[, cov] <- !complete.cases(X[!z, stripped_covs == cov])
    }
    # Correct denominator for amount of missingness
    lp_results$lpdetails$eps <- lp_results$lpdetails$raw_eps *
      Q / (Q - sapply(stripped_covs, function(con) {sum(lp_results$o$solution[1:N] * missing_instances[, con])}))
    # Corrected obj for LP
    lp_results$lpdetails$objective_wo_importances <- sum(lp_results$lpdetails$eps)
    lp_results$lpdetails$objective <- sum(importances * lp_results$lpdetails$eps)

    best_objective <- Inf
    run_raw_objectives <- rep(NA, runs)
    run_raw_objectives_wo_importances <- rep(NA, runs)
    run_objectives <- rep(NA, runs)
    run_objectives_wo_importances <- rep(NA, runs)

    if (is.null(seed)) {
      seed <- sample(1:1000000, 1)
    }
    set.seed(seed)

    for (run in 1:runs) {
      # Run randomized rounding
      rr_results_temp <- randomized_rounding(o = lp_results$o, N = N, z = z, st = st,
                                             st_vals = st_vals, S = S)

      # Calculate and format results

      # Epsilons for the randomized rounding (or integer) solution
      X0 <- X[z == 0, ] / sum(q_s)
      blk1 <- t(X0)
      raw_eps_temp <- matrix(0, nvars, 2)
      for (i in 1:nvars) {
        row <- blk1[i, ]
        inbalance <- sum(row * rr_results_temp$controls$select, na.rm = TRUE)
        if (inbalance < 0) {
          raw_eps_temp[i, 1] <- abs(inbalance)
        } else if (inbalance > 0) {
          raw_eps_temp[i, 2] <- inbalance
        }
      }
      if (!is.null(colnames(X))) {
        rownames(raw_eps_temp) <- colnames(X)
      }
      colnames(raw_eps_temp) <- c("positive", "negative")

      # Corrected epsilons for RR
      eps_temp <- raw_eps_temp *
        Q / (Q - sapply(stripped_covs, function(con) {sum(missing_instances[rr_results_temp$controls$select, con])}))

      # Objective value for the randomized rounding (or integer) solution
      run_raw_objectives[run] <- sum(importances * raw_eps_temp)
      run_raw_objectives_wo_importances[run] <- sum(raw_eps_temp)
      run_objectives[run] <- sum(importances * eps_temp)
      run_objectives_wo_importances[run] <- sum(eps_temp)

      if (run_raw_objectives[run] < best_objective) {
        raw_eps <- raw_eps_temp
        eps <- eps_temp
        raw_objective_wo_importances <- run_raw_objectives_wo_importances[run]
        objective_wo_importances <- run_objectives_wo_importances[run]
        objective <- run_objectives[run]
        rr_results <- rr_results_temp
        best_objective <- run_raw_objectives[run]
      }

    }

    return(list(objective = objective,
                objective_wo_importances = objective_wo_importances,
                eps = eps,
                importances = importances,
                selected = rr_results$selected,
                controls = rr_results$controls,
                rrdetails = list(
                  seed = seed,
                  raw_objective = best_objective,
                  raw_objective_wo_importances = raw_objective_wo_importances,
                  raw_eps = raw_eps,
                  run_raw_objectives = run_raw_objectives,
                  run_raw_objectives_wo_importances = run_raw_objectives_wo_importances,
                  run_objectives = run_objectives,
                  run_objectives_wo_importances = run_objectives_wo_importances),
                lpdetails = lp_results$lpdetails))
  }
}

#' Verify the inputs to \code{\link{optimize_controls}}
#'
#' Makes sure that the inputs to \code{\link{optimize_controls}} are in the correct
#' format and feasible.
#'
#' @inheritParams optimize_controls
#' @keywords internal

verify_inputs <- function(X, importances, ratio, q_s, st, z, integer, solver) {
  if (is.data.frame(X))
    X <- as.matrix(X)
  stopifnot(!is.null(ratio) || !is.null(q_s))
  stopifnot(is.null(ratio) | (is.vector(ratio) & (ratio > 0)))
  stopifnot(is.vector(st) | is.factor(st))
  stopifnot(is.vector(z))
  stopifnot(is.matrix(X))
  stopifnot(all((z == 0) | (z == 1)))
  stopifnot(length(z) == (dim(X)[1]))
  stopifnot(length(z) == length(st))
  stopifnot(is.logical(integer))
  if (!solver %in% c("gurobi", "Rglpk")) {
    stop("\"solver\" must be one of \"gurobi\" or \"Rglpk\".",
         call. = FALSE)
  }
  if (solver == "gurobi" && !requireNamespace("gurobi", quietly = TRUE)) {
    stop("Package \"gurobi\" needed if \"solver\" parameter set to \"gurobi\". Please
         install it or switch the \"solver\" parameter to \"Rglpk\".",
         call. = FALSE)
  }
  frtab <- table(z, st)
  if (min(frtab[2, ]) == 0) {
    warning("Note that at least one stratum has no treated individuals.")
  }
  stratios <- frtab[1, ]/frtab[2, ]
  st_vals <- as.character(colnames(frtab))  # stratum values
  if (!is.null(ratio) && min(stratios) < ratio) {
    stop("The ratio you specified is not feasible.
            Please supply `q_s` instead of `ratio` or lower the `ratio` input.",
         call. = FALSE)
  }
  if (!is.null(q_s)) {
    n_s <- table(z, st)[1, ]
    if (any(q_s[names(n_s)] > n_s)) {
      stop("At least one of the entries for `q_s` is greater than the number of controls available in the stratum.
           Please lower `q_s` such that all entries are at most the number of available controls.",
           call. = FALSE)
    }
  }
  stopifnot(sort(names(importances)) == sort(colnames(X)))
}


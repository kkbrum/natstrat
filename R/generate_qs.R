#' Calculate desired number of controls per stratum
#'
#' Figure out how many units to take from each stratum when some strata are deficient.
#' The result should be used as an input to \code{\link{optimize_controls}()}.
#'
#' @inheritParams stand
#' @param ratio a numeric specifying the desired ratio of controls to treated in
#'   each stratum.
#' @param max_ratio a numeric specifying the maximum ratio to allow in a stratum to achieve
#'   the overall \code{ratio} specified. If \code{NULL}, it is set by default to 1.1 times the desired \code{ratio}.
#'   To have no maximum ratio, set this to \code{Inf}.
#' @param max_extra_s single numeric or named vector with values corresponding to the maximum desired number
#'   of extra controls to be chosen from each stratum to achieve the overall
#'   \code{ratio} specified. If this is a vector, the names should correspond to the stratum
#'   values from \code{st}. The default is 5 for each stratum.
#'   To have no maximum, set this to \code{Inf}.
#'   If both \code{max_ratio} and \code{max_s} are specified, the
#'   maximum of the two will be used for each stratum.
#' @param strata_dist matrix with both row and column names with names corresponding
#'   to the stratum values from \code{st} and entries corresponding to the distance
#'   associated with taking a control from the stratum associated with the row when
#'   the desired stratum is the one associated with the column. Lower distance values
#'   are more desirable replacements. Typically the diagonal should be 0, meaning
#'   there is no penalty for choosing a unit from the correct stratum.
#' @return A named vector stating how many controls to take from each stratum.
#'
#' @export

generate_qs <- function(z, st, ratio = NULL,
                        max_ratio = NULL, max_extra_s = 5, strata_dist = NULL) {

  # Make sure inputs are good
  verify_inputs_EMD(ratio = ratio, st = st, z = z)

  # Look at strata counts
  frtab <- table(z, st)
  if (min(frtab[2, ]) == 0) {
    warning("Note that at least one stratum has no treated individuals.")
  }
  stratios <- frtab[1, ]/frtab[2, ]
  if (is.null(ratio)) {
    ratio <- min(1, min(stratios))
  }
  st_vals <- as.character(colnames(frtab))  # stratum values
  S <- length(st_vals)  # number of strata
  # Determine number of controls desired for each stratum
  desired_qs <- round(ratio * frtab[2, ])
  # Number of controls available by strata
  n_s <- frtab[1, ]

  # Make sure the maximum per strata input is okay
  if(!is.null(max_extra_s)) {
    if (length(max_extra_s) > 1) {
      if (length(max_extra_s) != S) {
        stop("Length of \"max_extra_s\" must be 1 or match number of strata occuring in \"st\"",
             call. = FALSE)
      }
      if (is.null(names(max_extra_s))) {
        stop("If \"max_extra_s\" is a vector, it must have names corresponding to the stratum values.")
      }
      max_extra_s <- max_extra_s[st_vals]
    }
  } else {
    max_extra_s <- 5
  }
  max_s <- desired_qs + max_extra_s

  if (is.null(max_ratio)) {
    max_ratio <- 1.1 * ratio
  }
  if (max_ratio < Inf) {
    max_s <- pmax(max_s, round(max_ratio * frtab[2, ]))
  } else {
    max_s <- pmax(max_s, Inf)
  }

  max_s <- pmin(max_s, n_s)

  if (sum(max_s) < sum(desired_qs)) {
    stop("The \"max_ratio\" or \"max_extra_s\" supplied is not feasible because
           it does not allow for enough extra controls in the strata where they are available.
           Please increase these maximums or reduce the desired ratio.")
  }

  # Solve for minimum work value of earth-movers distance if not enough controls in every stratum
  if (is.null(strata_dist)) {
    # No preference for which controls are selected
    strata_dist <- matrix(1, nrow = S, ncol = S)
    diag(strata_dist) <- 0
  } else {
    if (is.null(rownames(strata_dist)) | is.null(colnames(strata_dist))) {
      stop("\"strata_dist\" must have both column and row names corresponding to the stratum values.")
    }
    strata_dist <- strata_dist[st_vals, st_vals]
  }
  if (min(stratios) < ratio) {
    # Find how many controls to take from each stratum using earth-movers distance
    strata_dist_flat <- matrix(t(strata_dist), nrow = 1)
    colnames(strata_dist_flat) <- rep(st_vals, each = S)
    q_s <- presolve_EMD(S, desired_qs, max_s, strata_dist_flat)
  } else { q_s <- desired_qs } # No need for earth-movers distance if there are enough controls per strata

  return(q_s)
}

#' Verify the inputs to the earthmover's distance problem
#'
#' Check that the ratio, strata, and treated indicator provided
#' are in the correct forms and that the desired ratio is feasible
#' across the population.
#'
#' @inheritParams generate_qs
#' @keywords internal


verify_inputs_EMD <- function(ratio, st, z) {
  stopifnot(is.null(ratio) | (is.vector(ratio) & (ratio > 0)))
  stopifnot(is.vector(st) | is.factor(st))
  stopifnot(is.vector(z))
  stopifnot(all((z == 0) | (z == 1)))
  stopifnot(length(z) == length(st))
  N <- sum(1 - z)  # Total number of controls available
  Q <- ratio * sum(z)  # Total number of controls to be selected
  if (Q > N) {
    stop(paste0("There are not enough total controls to use the ratio you provided.
               The maximum ratio you may have is ", N/sum(z),
                " which will select all treated and control units."),
         call. = FALSE)
  }
}


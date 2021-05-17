#' Standardize covariate vector for balance constraint
#'
#' This function is used by \code{\link{generate_constraints}()} to standardize
#' covariate vectors to become balance constraints. This standardization
#' is done such that the balance constraint will be minimized when the treated
#' and control groups (within or across strata) have equal means.
#' The function subtracts the treated mean (within or across strata) and
#' divides by the treated or pooled standard deviation (across strata).
#'
#' @param z a treatment indicator vector with \code{i}th entry equal to 0 if
#'   unit \code{i} is a control and equal to 1 if unit \code{i} is treated.
#' @param x a covariate vector with \code{i}th entry equal to the
#'   covariate value of unit \code{i}. This should have the same order of units and
#'   length as \code{z}.
#' @param st a stratum vector with the \code{i}th entry equal to the
#'   stratum of unit \code{i}. This should have the same order of units and length
#'   as \code{z}.
#' @param ist an optional specification of the
#'   target stratum within which the balance constraint is desired. Must be one of the
#'   values within \code{st}. By default, this is \code{NULL}, meaning the generated
#'   constraint balances across strata.
#' @param denom_variance character stating what variance to use in the standardization:
#'   either the default "treated", meaning the standardization will use the
#'   treated variance (across all strata), or "pooled", meaning
#'   the standardization will use the average of the treated and control variances.
#' @param autogen_missing whether to automatically generate missingness constraint
#'   and how heavily to prioritize it. Should be a numeric
#'   or \code{NULL} value. \code{NULL} indicates that
#'   a constraint to balance the rate of missingness (denoted by \code{NA}
#'   in \code{x}) should not be automatically generated. Note that this is not
#'   recommended unless the user has already accounted for missing values.
#'   If not \code{NULL}, \code{autogen_missing} should be a numeric stating how heavily
#'   to prioritize generated missingness constraints over covariate constraint.
#'   The default is 50.
#' @return A list with two components:
#' \describe{
#' \item{covariate}{a balance constraint for the standardized covariate values
#'   of either all treated and control units,
#'   or just those in stratum \code{ist}.}
#' \item{missingness}{a corresponding balance constraint for the rate of missingness if
#' \code{autogen_missing} not \code{NULL}, otherwise \code{NULL}.}
#' }
#' @importFrom stats sd
#' @export

stand <- function(z, x, st, ist = NULL, denom_variance = "treated", autogen_missing = 50) {

  # Verify inputs ----
  if (!denom_variance %in% c("treated", "pooled")) {
    stop("* `denom_variance` must be one of `treated` or `pooled`.",
         call. = FALSE)
  }
  if (!(is.vector(z) && is.vector(x) && is.vector(st))) {
    stop("`z`, `x`, and `st` must all be vectors.",
         call. = FALSE)
  }
  if (!all((z == 0) | (z == 1))) {
    stop("`z` must contain all 0s and 1s.",
         call. = FALSE)
  }
  if (length(z) != length(x) || length(z) != length(st)) {
    stop("`z`, `x`, and `st` must all have same length.",
         call. = FALSE)
  }
  if (!is.null(ist) && (!ist %in% unique(st))) {
    stop("`ist` must be one of the values in `st`.",
         call. = FALSE)
  }

  # Define the units for which to generate the constraint ----
  if (is.null(ist)) {
    ind <- rep(TRUE, length(z))
  } else {
    ind <- (st == ist)
  }

  # Define the target mean ----
  loc <- mean(x[ind & z == 1], na.rm = TRUE)
  if (is.na(loc)) {
    stop("At least one stratum of a covariate has no treated units with nonmissing values.",
         call. = FALSE)
  }

  # Define the value by which to scale ----
  if (denom_variance == "treated") {
    scl <- sd(x[z == 1], na.rm = TRUE)
    # If there is no variance in treated group, use pooled value
    # (which is half of the control variance since treated variance = 0)
    if (is.na(scl) || scl == 0) {
      warning("There is a covariate with no variance in the treated group. Standardization will thus use the average of the treated and control variance for this covariate.")
      scl <- sd(x[z == 0], na.rm = TRUE) / sqrt(2)
    }
  } else if (denom_variance == "pooled") {
    scl <- sqrt((sd(x[z == 0], na.rm = TRUE)^2 + sd(x[z == 1], na.rm = TRUE)^2) / 2)
  }

  # Perform standardization ----
  x_stand <- (x[ind] - loc) / scl
  x_stand[x[ind] == loc] <- 0  # Catches the 0/0 case if no variance in treated and control group

  # Deal with missing data (coded as NA) ----
  # Add missingness constraint
  miss_stand <- NULL
  if (!is.null(autogen_missing) && sum(is.na(x)) > 0) {
    miss <- is.na(x)
    loc_miss <- mean(miss[ind & z == 1])
    if (denom_variance == "treated") {
      scl_miss <- sd(miss[z == 1])
      if (is.na(scl_miss) || scl_miss == 0) {
        scl_miss <- sd(miss[z == 0]) / sqrt(2)
      }
    } else if (denom_variance == "pooled") {
      scl_miss <- (sd(miss[z == 0]) + sd(miss[z == 1])) / 2
    }
    miss_stand <- (miss[ind] - loc_miss) / scl_miss
    miss_stand[is.na(miss_stand)] <- 0  # Catches the 0/0 case if no variance in treated and control group
  }

  return(list("covariate" = x_stand, "missingness" = miss_stand))
}


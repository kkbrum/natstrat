#' Sample integer solution from linear programming solution
#'
#' The linear programming solution of \code{\link{balance_LP}} that is used
#' within \code{\link{optimize_controls}} sometimes selects fractional control units.
#' Here, we select any unit the linear programming solution chose with coefficient 1.
#' Then, we select the remaining required number of units from those that have
#' fractional solutions by sampling with probabilities equal to the linear
#' programming solution. Used within \code{\link{optimize_controls}}.
#'
#' @inheritParams optimize_controls
#' @inheritParams balance_LP
#' @param o linear programming results, as found in the `o` element of the
#' returned list from \code{\link{balance_LP}}.
#'
#' @return List containing:
#' \describe{
#' \item{\code{controls}}{dataframe with two columns: \code{pr}, which contains
#' the coefficient determined for that control unit from the linear programming
#' solution, and \code{select}, a boolean vector stating whether that
#' control unit was selected for inclusion by randomized rounding.}
#' \item{\code{selected}}{boolean vector stating whether each unit (both treated
#' and control) was selected.}
#' }
#'
#' @keywords internal

randomized_rounding <- function(o, N, z, st, st_vals, S) {

  pr <- o$solution[1:N]
  pr <- round(pr, 5)
  select <- rep(FALSE, N)
  select[pr == 1] <- TRUE
  ind <- 1:N
  stc <- st[z == 0]  # Strata for controls
  for (j in 1:S) {
    w <- (pr < 1) & (pr > 0) & (stc == st_vals[j])
    needed <- sum(pr[pr != 1 & (stc == st_vals[j])])  # Sample how many? Total weight of non-1 entries
    needed <- round(needed, 3)
    if (needed > 0) {
      if (length(ind[w]) == 1) {
        sa <- ind[w]
      } else if (needed == 1) {
        sa <- sample(ind[w], needed, replace = FALSE, prob = pr[w])
      } else {
        if (sum(w) < 100) {
          sampled <- pps::sampford(pr[w], needed)
        } else {
          sampled <- sampling::UPmidzuno(pr[w])
          sampled <- sampled == 1
        }
        sa <- (ind[w])[sampled]
      }
      select[sa] <- TRUE
    }
  }
  controls <- data.frame(pr, select)
  selected <- rep(TRUE, length(z))
  selected[z == 0] <- select

  return(list(controls = controls, selected = selected))
}

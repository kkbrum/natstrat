#' Sample integer solution from linear programming solution with sample sizes correct in expectation
#'
#' The linear programming solution of \code{\link{balance_LP}()} that is used
#' within \code{\link{optimize_controls}()} sometimes selects fractional control units.
#' Here, we select any unit the linear programming solution chose with coefficient 1.
#' Then, we select sample each unit with a fractional solution with
#' probability equal to the linear programming solution. The total sample
#' size is then correct in expectation. Used within \code{\link{optimize_controls}()}
#' if \code{correct_sizes = FALSE}.
#'
#' @inheritParams optimize_controls
#' @inheritParams balance_LP
#' @param o linear programming results, as found in the `o` element of the
#' returned list from \code{\link{balance_LP}()}.
#'
#' @return Dataframe with two columns: \code{pr}, which contains
#' the coefficient determined for that unit from the linear programming
#' solution, and \code{select}, a boolean vector stating whether that
#' unit was selected for inclusion by randomized rounding.
#'
#' @keywords internal
#' @importFrom stats rbinom rmultinom

randomized_rounding_expectation <- function(o, N, multi_comp = FALSE) {
  if (multi_comp) {
    pr <- o$solution[1:(2 * N)]
  } else {
    pr <- o$solution[1:N]
  }
  pr <- round(pr, 5)
  select <- rep(FALSE, length(pr))
  for (i in 1:N) {
    if (multi_comp) {
      draw <- rmultinom(1, 1, c(pr[i], pr[N + i], round(1 - pr[i] - pr[N + i], 5)))
      if (draw[1] == 1) {
        select[i] <- TRUE
      } else if (draw[2] == 1) {
        select[N + i] <- TRUE
      }
    } else {
      draw <- rbinom(1, 1, pr[i])
      if (draw[1] == 1) {
        select[i] <- TRUE
      }
    }
  }

  units <- data.frame(pr, select)

  return(units = units)
}

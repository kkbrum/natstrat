#' @description
#' Natural strata fix a constant ratio of controls to treated units within
#' each stratum. This ratio need not be an integer. The control units are
#' chosen using randomized rounding of a linear program that balances many
#' covariates.
#' The gurobi commercial optimization software is recommended, but not required.
#' If available, the gurobi R package can be installed following the instructions
#' at the useful link about gurobi below.
#'
#' @details
#' To achieve the desired ratio of control to treated units,
#' a subset of control units are
#' chosen using by optimizing the balance of many covariates using
#' either randomized rounding of a linear program or
#' an integer program. The main function in this package is \code{\link{optimize_controls}}.
#' To create the input constraints for this function, you should use
#' \code{\link{generate_constraints}}.
#'
#' @references Brumberg, K; Small, D. S.; Rosenbaum, P. R. (2021+),
#' "Using Randomized Rounding of Linear Programs to Obtain Unweighted
#' Natural Strata That Balance Many Covariates".
#' @keywords internal
"_PACKAGE"

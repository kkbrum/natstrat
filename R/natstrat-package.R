#' @description
#' Natural strata can be used in observational studies to balance
#' the distributions of many covariates across any number of treatment
#' groups and any number of comparisons. These strata have proportional
#' amounts of units within each stratum across the treatments, allowing
#' for simple interpretation and aggregation across strata. Within each
#' stratum, the units are chosen using randomized rounding of a linear
#' program that balances many covariates. For more details, see Brumberg et al. (2022) <doi:10.1111/rssa.12848> and Brumberg et al.(2023) <doi:10.1093/jrsssc/qlad010>.
#' To solve the linear program, the 'Gurobi' commercial optimization software
#' is recommended, but not required. The 'gurobi' R package can be installed by following the instructions
#' \href{https://docs.gurobi.com/projects/optimizer/en/current/reference/r/setup.html}{here} after claiming your free academic license 
#' \href{https://www.gurobi.com/academia/academic-program-and-licenses/}{here}.
#'
#' @details
#' To achieve the desired ratio of control to treated units,
#' a subset of control units are
#' chosen using by optimizing the balance of many covariates using
#' either randomized rounding of a linear program or
#' an integer program. The main function in this package is \code{\link{optimize_controls}()}.
#' To create the input constraints for this function, you should use
#' \code{\link{generate_constraints}()}.
#'
#' @keywords internal
"_PACKAGE"

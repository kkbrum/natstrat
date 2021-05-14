z <- c(rep(0, 10), rep(1, 5))
data <- data.frame(color = c(rep("Red", 5), rep("White", 2), rep("Blue", 3), rep("White", 2), rep("Red", 3)),
                   number = 1:15,
                   category = c(rep(c("1", "2"), 5), "1", rep("2", 3), "1"))
data$number[c(1, 5, 11)] <- NA
constraints <- suppressWarnings(generate_constraints(list(color + number ~ 2 * category), z, data = data,
                                                     autogen_missing = 4))
results <- optimize_controls(z = z, X = constraints$X, st = data$category, ratio = 1.5,
                             importances = constraints$importances,
                             integer = FALSE, solver = "Rglpk", seed = 1, runs = 5,
                             time_limit = Inf)

test_that("optimization gives correct raw results", {
  expect_equal(results$lpdetails$raw_objective, 16.76667, tolerance = .000001)
  expect_equal(results$rrdetails$raw_objective, 17.27342, tolerance = .000001)
})

test_that("optimization gives correct corrected results", {
  expect_equal(results$lpdetails$objective, 20.47115, tolerance = .000001)
  expect_equal(results$objective, 22.3415, tolerance = .000001)
})

test_that("optimization chooses correct number of units", {
  expect_equal(sum(results$controls$pr[results$controls$pr < 1 & results$controls$pr > 0]),
               sum(results$controls$select[results$controls$pr < 1 & results$controls$pr > 0]))
  expect_equal(sum(results$controls$select), sum(round(table(z, data$category) * 1.5) [2,]))
})

test_that("epsilon corrected properly", {
  expect_equal(results$rrdetails$raw_eps['number_category1', ] * sum(results$controls$select) /
                 (sum(results$controls$select) - sum(is.na(constraints$X[results$selected & !z, "number_category1"]))),
               results$eps['number_category1', ])
  expect_equal(results$objective, sum(results$importances * results$eps))
  expect_equal(results$objective_wo_importances, sum(results$eps))
  expect_equal(results$lpdetails$raw_eps['number_category1', ] * sum(results$controls$pr) /
                 (sum(results$controls$pr) - sum(results$controls$pr * is.na(constraints$X[!z, "number_category1"]))),
               results$lpdetails$eps['number_category1', ])
  expect_equal(results$lpdetails$objective, sum(results$importances * results$lpdetails$eps))
  expect_equal(results$lpdetails$objective_wo_importances, sum(results$lpdetails$eps))

})

data_for_sds <- cbind(data[, 1:2], is.na(data$number))
names(data_for_sds)[3] <- "number_missing"
sds <- check_balance(z, data_for_sds, data$category, results$selected, message = FALSE)

test_that("epsilons across strata equal standardized diff in means across", {
  expect_equal(sds$sd_across[, "abs_stand_diff_after"],
               as.numeric(rowSums(results$eps[paste0(row.names(sds$sd_across), "_1"),])))
})

test_that("sum of epsilons within strata equal weighted avg of within strata standardized diff in means", {
  covs <- row.names(sds$sd_across)
  stripped_row_names <- sapply(strsplit(row.names(results$eps), "_"),
                               function(k) {paste0(k[-length(k)], collapse = "_")})
  sum_eps_within <- sapply(covs, function(cov) {
    sum(results$eps[stripped_row_names == cov &
                            grepl("_category", row.names(results$eps))])})
  expect_equal(sds$sd_strata_avg[, "abs_stand_diff_after"],
               as.numeric(sum_eps_within))

})


z <- c(rep(0, 15), rep(1, 5))
data <- data.frame(color = c(rep("Red", 5), rep("White", 2), rep("Blue", 5), rep("White", 4), rep("Red", 4)),
                   number = 1:20,
                   category = c("1", "1", "1", rep(c( "2", "3"), 6), "1", rep("2", 2), "3", "1"))
data$number[c(1, 5, 11, 16)] <- NA
constraints <- suppressWarnings(generate_constraints(list(color + number ~ 2 * category), z, data = data,
                                                     autogen_missing = 4))
strata_dist <- matrix(c(0, 2, 1, 2, 0, 1, 1, 1, 0), byrow = TRUE, ncol = 3, dimnames = list(c("1", "2", "3"), c("1", "2", "3")))
q_s <- generate_qs(z, st = data$category, ratio = 2.5, max_ratio = NULL, max_extra_s = NULL, strata_dist = strata_dist)
results_emd <- suppressWarnings(optimize_controls(z = z, X = constraints$X, st = data$category, q_s = q_s,
                                                  importances = constraints$importances,
                                                  integer = FALSE, solver = "Rglpk", seed = 1, runs = 5,
                                                  time_limit = Inf))


test_that("EMD chooses correct number of units", {
  expect_equal(sum(results_emd$controls$select), sum(round(2.5 * table(z, data$category)[2, ])))
})

test_that("EMD gives right objective", {
  expect_equal(results_emd$rrdetails$raw_objective, 29.31395, tolerance = 0.0001)
  expect_equal(results_emd$objective, 34.78992, tolerance = 0.0001)
})

test_that("Choosing from closest strata", {
  expect_equal(as.numeric(table(data$category[results_emd$selected & !z]))[3], 4)
})

test_that("Specified max ratio and max extra working", {
  q_s2 <- generate_qs(z, st = data$category, ratio = 2.5, max_ratio = 0, max_extra_s = 1, strata_dist = strata_dist)
  results_emd2 <- suppressWarnings(optimize_controls(z = z, X = constraints$X, st = data$category, q_s = q_s2,
                                                    importances = constraints$importances,
                                                    integer = FALSE, solver = "Rglpk", seed = 1, runs = 5,
                                                    time_limit = Inf))
  expect_equal(as.numeric(table(data$category[results_emd2$selected & !z]))[3], 3)
  q_s3 <- generate_qs(z, st = data$category, ratio = 2.5, max_ratio = Inf, max_extra_s = 1, strata_dist = strata_dist)
  results_emd3 <- suppressWarnings(optimize_controls(z = z, X = constraints$X, st = data$category,  q_s = q_s3,
                                                     importances = constraints$importances,
                                                     integer = FALSE, solver = "Rglpk", seed = 1, runs = 5,
                                                     time_limit = Inf))
  expect_equal(as.numeric(table(data$category[results_emd3$selected & !z]))[3], 4)
  q_s4 <- generate_qs(z, st = data$category, ratio = 2.5,  max_ratio = 3, max_extra_s = 0, strata_dist = strata_dist)
  results_emd4 <- suppressWarnings(optimize_controls(z = z, X = constraints$X, st = data$category, q_s = q_s4,
                                                     importances = constraints$importances,
                                                     integer = FALSE, solver = "Rglpk", seed = 1, runs = 5,
                                                     time_limit = Inf))
  expect_equal(as.numeric(table(data$category[results_emd4$selected & !z]))[3], 3)
})



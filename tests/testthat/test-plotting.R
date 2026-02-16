# Tests for visualization helper functions

# Shared test data
make_ors_results <- function() {
  set.seed(42)
  data.frame(
    gene_sp1 = paste0("SP1_G", 1:20),
    gene_sp2 = paste0("SP2_G", 1:20),
    CCS = runif(20, -0.5, 0.9),
    ORS = runif(20, 0, 1),
    logORS = rnorm(20, mean = 1, sd = 1),
    type = rep(c("1:1", "1:N", "N:1", "N:M"), each = 5),
    stringsAsFactors = FALSE
  )
}


test_that("plot_ors_distribution returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  ors <- make_ors_results()
  p <- plot_ors_distribution(ors)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ors_distribution validates input", {
  skip_if_not_installed("ggplot2")
  expect_error(plot_ors_distribution(data.frame(x = 1)), "logORS")
})

test_that("plot_ors_by_type returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  ors <- make_ors_results()
  p <- plot_ors_by_type(ors)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ors_by_type requires type column", {
  skip_if_not_installed("ggplot2")
  ors <- make_ors_results()
  ors$type <- NULL
  expect_error(plot_ors_by_type(ors), "type")
})

test_that("plot_ccs_vs_ors returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  ors <- make_ors_results()
  p <- plot_ccs_vs_ors(ors)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ccs_vs_ors validates input", {
  skip_if_not_installed("ggplot2")
  ors <- make_ors_results()
  ors$CCS <- NULL
  expect_error(plot_ccs_vs_ors(ors), "CCS")
})

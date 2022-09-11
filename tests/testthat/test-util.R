
test_that("renaming works", {
  df1 <- data.frame(x = 1, y = "b")
  df2 <- rename_check(df1, "y", "z")
  expect_equal(names(df2)[2], "z")
})

test_that("renaming fails because variable not found", {
  df1 <- data.frame(x = 1, y = "b")
  expect_error(rename_check(df1, "k", "z"))
})

test_that("renaming fails because more than one variable is found", {
  df1 <- data.frame(y = 1, y = "b")
  names(df1)[2] <- "y"
  expect_error(rename_check(df1, "y", "z"))
})

test_that("get_dat_ordered1 renames correctly", {
  df <- structure(list(xx = c(0.87, 0.02, 0.46, 0.57, 0.83, 0.24), 
                       yx = c(6.14, 6.11, 9.9, NA, NA, 14.17), 
                       cutx = c(NA, NA, NA, 4, 5.81, NA), 
                       uncensoredx = c(1, 1, 1, 0, 0, 1)), 
                  row.names = c(1L, 2L, 3L, 4L, 5L, 6L), class = "data.frame")
  result <- get_dat_ordered1(df, x = "xx", y = "yx", uncensored = "uncensoredx", threshold = "cutx")
  expect_equal(names(result), c("x", "y", "cut", "uncensored", "y_comb"))
})

test_that("get_dat_ordered1 orders correctly", {
  df <- structure(list(xx = c(0.87, 0.02, 0.46, 0.57, 0.83, 0.24), 
                       yx = c(6.14, 6.11, 9.9, NA, NA, 14.17), 
                       cutx = c(NA, NA, NA, 4, 5.81, NA), 
                       uncensoredx = c(1, 1, 1, 0, 0, 1)), 
                  row.names = c(1L, 2L, 3L, 4L, 5L, 6L), class = "data.frame")
  result <- get_dat_ordered1(df, x = "xx", y = "yx", uncensored = "uncensoredx", threshold = "cutx")
  expect_equal(result$uncensored, c(1,1,1,1,0,0))
})



# test_data_base <- data.frame(
#   x = rep(2009:2020, each = 3),
#   threshold = NA,
#   uncensored = 1
# )
# saveRDS(test_data_base, "tests/testthat/fixtures/test_data_clean.rds")


test_that("lc_clean1 - file doesn't need cleaning", {
  test_data_base <- readRDS(test_path("fixtures", "test_data_clean.rds"))
  test_data <- set_x_to_censored(test_data_base, c(2010:2011, 2013, 2016))
  test_result <- lc_clean1(test_data)  
  expect_equal(test_data, test_result)
})

test_that("lc_clean1 - clean 4 years", {
  test_data_base <- readRDS(test_path("fixtures", "test_data_clean.rds"))
  test_data <- set_x_to_censored(test_data_base, c(2010:2016))
  test_result <- lc_clean1(test_data) 
  expect_equal(subset(test_data, x>=2013), test_result)
})

test_that("lc_clean1 - clean almost all", {
  test_data_base <- readRDS(test_path("fixtures", "test_data_clean.rds"))
  test_data <- set_x_to_censored(test_data_base, c(2010:2011, 2013, 2016:2019))
  test_result <- lc_clean1(test_data)  
  expect_equal(subset(test_data, x>=2019), test_result)
})

test_that("lc_clean2 - file doesn't need cleaning", {
  test_data_base <- readRDS(test_path("fixtures", "test_data_clean.rds"))
  test_data <- set_x_to_censored(test_data_base, c(2010:2019))
  test_result <- lc_clean2(test_data)  
  expect_equal(test_data, test_result)
})

test_that("lc_clean2 - delete one year", {
  test_data_base <- readRDS(test_path("fixtures", "test_data_clean.rds"))
  test_data <- set_x_to_censored(test_data_base, c(2009, 2019))
  test_result <- lc_clean2(test_data)  
  expect_equal(length(unique(test_result$x)), 11)
})

test_that("lc_clean2 - delete one year", {
  test_data_base <- readRDS(test_path("fixtures", "test_data_clean.rds"))
  test_data <- set_x_to_censored(test_data_base, c(2009, 2019))
  test_result <- lc_clean2(test_data)  
  expect_equal(subset(test_data, x>=2010), test_result)
})

test_that("lc_clean2 - delete 4 years", {
  test_data_base <- readRDS(test_path("fixtures", "test_data_clean.rds"))
  test_data <- set_x_to_censored(test_data_base, c(2009:2012, 2019))
  test_result <- lc_clean2(test_data)  
  expect_equal(subset(test_data, x>=2013), test_result)
})

test_that("check_large_m behaves as expected", {
  expect_no_warning(check_large_m(warn = TRUE,  p = 10, m = 4))
  
  expect_warning(
    check_large_m(warn = TRUE, p = 10, m = 5),
    class = "TooLarge_m_Warning"
  )
  
  expect_no_warning(check_large_m(warn = FALSE, p = 10, m = 5))
  
  expect_error(check_large_m(warn = FALSE, p = 10, m = 5.2),
               class = "InvalidArgumentError")
  expect_error(check_large_m(warn = FALSE, p = 10, m = NULL),
               class = "InvalidArgumentError")
  expect_error(check_large_m(warn = FALSE, p = 10, m = NA),
               class = "InvalidArgumentError")
})


library(testthat)
library(DoseFinding)
options(testthat.progress.max_fails = 100)

Sys.unsetenv("R_TESTS")
test_check("DoseFinding")

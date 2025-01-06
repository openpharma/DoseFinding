
# Generating test data
data(biom)
models <- Mods(linear = NULL, emax=c(0.05,0.2), linInt=c(1, 1, 1, 1), doses=c(0,0.05,0.2,0.6,1))
MM <- MCPMod(dose, resp, biom, models, Delta=0.5)

test_that("MCPMod object can be printed", {
  expect_output(print(MM), "MCPMod\\n")
  expect_output(print(MM), "Multiple Contrast Test:\\n")
  expect_output(print(MM), "Estimated Dose Response Models:")
})

test_that("summary.MCPMod summarizes and prints an MCPMod object", {
  expect_output( summary(MM), "MCP part \\n")
  expect_output( summary(MM), "Mod part \\n")
  expect_output( summary(MM), "Model selection criteria \\(AIC\\):")
  expect_output( summary(MM), "Estimated TD\\, Delta=0\\.5\\n")
})

test_that("plot.MCPMod plots the fitted dose-response model", {
  expect_silent(plot(MM, plotData = "meansCI"))
  expect_silent(plot(MM, plotData = "means"))
  expect_silent(plot(MM, plotData = "raw"))
  expect_silent(plot(MM, plotData = "none"))
})

test_that("predict.MCPMod provides predictions from the fitted dose-response model", {
  pred <- predict(MM, se.fit = TRUE, doseSeq = c(0,0.2,0.4, 0.9, 1), predType="ls-means")
  expect_true(is.list(pred))
  expect_true(is.list(pred[[1]]))  # Ensure each model provides a list
})

test_that("plot.MCPMod stops with appropriate error when no models significant", {
  # Create a scenario where no models are significant
  models_no_sig <- Mods(linear = NULL, doses=c(0,0.05,0.2,0.6,1))
  MM_no_sig <- MCPMod(dose, resp, biom, models_no_sig, Delta=0.5, critV = 9999)
  expect_error(plot(MM_no_sig))
})
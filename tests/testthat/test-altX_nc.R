skip_on_cran()

################################################################################
################################################################################
### file test-altX_nc.R
################################################################################
################################################################################

################################################################################
### Check altX_nc effect ###
################################################################################

mynet <- sienaDependent(array(c(s502, s503), dim=c(50, 50, 2)))
mybeh <- sienaDependent(s50a[, 2:3], type = "behavior")
mydata <- sienaDataCreate(mynet, mybeh)

altx_model <- getEffects(mydata)
altx_model <- setEffect(altx_model, altX,
	name = "mynet", interaction1 = "mybeh");
altx_controls <- sienaAlgorithmCreate(projname = NULL, seed = 42)
altx_ans <- siena07(
	altx_controls,
	batch = TRUE,
	silent = TRUE,
	data = mydata,
	effects = altx_model,
	returnChains = FALSE);

altx_nc_model <- getEffects(mydata)
altx_nc_model <- setEffect(altx_nc_model, altX_nc,
	name = "mynet", interaction1 = "mybeh");
altx_nc_controls <- sienaAlgorithmCreate(projname = NULL, seed = 42)
altx_nc_ans <- siena07(
	altx_nc_controls,
	batch = TRUE,
	silent = TRUE,
	data = mydata,
	effects = altx_nc_model,
	returnChains = FALSE);


test_that("Target statistics are correct", {
	adj <- mynet[, , 2]
	beh <- mybeh[, , 1]
	beh_centered <- beh - mean(mybeh)

	altx_nc_target <- sum(adj * rep(beh, each = nrow(adj)))
	altx_target <- sum(adj * rep(beh_centered, each = nrow(adj)))

	expect_equal(altx_target, altx_ans$targets[4])
	expect_equal(altx_nc_target, altx_nc_ans$targets[4])
	expect_false(isTRUE(all.equal(altx_ans$targets[4], altx_nc_ans$targets[4])))
})

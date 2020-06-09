test_that("Missing group", {
    expect_error(create_standata_exposure(unit_id=1:10,
                                       conc=1:10), "The lengths")
})


test_that("Missing unit_id", {
    expect_error(create_standata_exposure(group=rep(1, 10),
                                          conc=1:10), "The lengths")
})


test_that("Missing conc", {
    expect_error(create_standata_exposure(group=rep(1, 10),
                                          unit_id=1:10), "The lengths")
})

test_that("Incorrect length for conc.", {
    expect_error(create_standata_exposure(group=rep(1, 10),
                                          unit_id=1:10,
                                          conc=1:5), "The lengths")
})

test_that("Incorrect length for unit_id", {
    expect_error(create_standata_exposure(group=rep(1, 10),
                                          unit_id=1:5,
                                          conc=1:10), "The lengths")
})

test_that("Incorrect length for unit_id", {
    expect_error(create_standata_exposure(group=1,
                                          unit_id=1:5,
                                          conc=1:10), "The lengths")
})


test_that("Incorrect length for clust_id", {
    expect_error(create_standata_exposure(group=rep(1, 10),
                                          unit_id=1:10,
                                          conc=1:10,
                                          clust_id=1), "The lengths")
})

test_that("Incorrect length for time trend", {
    expect_error(create_standata_exposure(group=rep(1, 10),
                                          unit_id=1:10,
                                          conc=1:10,
                                          time=0:1), "The lengths")
})


test_that("Missing clust_id", {
    expect_is(create_standata_exposure(group=1:10,
                                          unit_id=1:10,
                                          conc=1:10), "standata_exposure")
    expect_is(create_standata_exposure(group=1:10,
                                       unit_id=1:10,
                                       conc=1:10,
                                       clust_id=0), "standata_exposure")
})

test_that("Create exposure standata object", {
    expect_s3_class(create_standata_exposure(group=1:10,
                                       unit_id=1:10,
                                       conc=1:10), "standata_exposure")

})

test_that("Simple exposure model fit without prior", {
    exp_standata2 <- create_standata_exposure(group=rep(1, 10),
                                             unit_id=1:10,
                                             conc=1:10)
    expect_message(sample_exposure_model(exp_standata2, B=100), "failed to create the sampler")
})

test_that("Simple exposure model fit", {
    exp_standata <- create_standata_exposure(group=rep(1, 20),
                                             unit_id=rep(1:10, each=2),
                                             conc=rnorm(20))
    exp_standata <- add_priors(exp_standata)
    expect_s4_class(sample_exposure_model(exp_standata, B=500, chains=2), "stanfit")
})

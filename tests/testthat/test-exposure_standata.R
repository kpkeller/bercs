test_that("Mismatched lengths", {
    expect_error(create_standata_exposure(group=1:10,
                                          unit_id=1:10,
                                          conc=1:5), "The lengths")
    expect_error(create_standata_exposure(group=1:10,
                                          unit_id=1:5,
                                          conc=1:10), "The lengths")
    expect_error(create_standata_exposure(group=1:10,
                                          unit_id=1:10,
                                          conc=1:10,
                                          clust_id=1:5), "The lengths")
    expect_error(create_standata_exposure(group=1:10,
                                          unit_id=1:10,
                                          conc=1:10,
                                          clust_id=c(0, 0),
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

test_that("Missing time", {
    expect_is(create_standata_exposure(group=1:10,
                                       unit_id=1:10,
                                       conc=1:10), "standata_exposure")
    expect_is(create_standata_exposure(group=1:10,
                                       unit_id=1:10,
                                       conc=1:10,
                                       ), "standata_exposure")
})

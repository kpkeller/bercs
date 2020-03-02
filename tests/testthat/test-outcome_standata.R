test_that("Mismatched lengths", {
    expect_error(create_standata_outcome_singlestudy(study=1:10,
                                                     unit_id=1:10,
                                                     conc=1:5), "The lengths")
})

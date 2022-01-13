

test_that("Create simple outcome simulation object.", {

    expect_s3_class(create_outcome_simulation_skeleton(nunits=10),
                                                        "outsim")

})



test_that("Mismatched lengths in outsim.", {
    expect_error(create_outcome_simulation_skeleton(nunits=10,
                                                     x=1:3), "The length")
})

outsim_update_covariate


test_that("Outcome simulation update covariates.", {
    os <- create_outcome_simulation_skeleton(nunits=10)
    expect_error(outsim_update_covariate(os, x=2),
                    "unused")
    expect_error(outsim_update_covariate(os, covariate="x",
                                            value=1:7),
                    "same length")

})



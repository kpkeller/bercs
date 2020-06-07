test_that("Mismatched lengths", {
    expect_error(create_standata_outcome_singlestudy(unit_id=1:10,
                                                     conc=1:5,
                                                     case=rep(0:1, each=5)), "The lengths")
    expect_error(create_standata_outcome_singlestudy(unit_id=1:10,
                                                     conc=1:10,
                                                     case=rep(0, 7)), "The lengths")
    expect_error(create_standata_outcome_singlestudy(unit_id=1:7,
                                                     conc=1:10,
                                                     case=rep(0, 10)), "The lengths")
    expect_error(create_standata_outcome_singlestudy(study=1,
                                                     unit_id=1:10,
                                                     conc=1:10,
                                                     case=rep(0, 10)), "The lengths")
})

test_that("Create simple outcome standata object.", {


    expect_s3_class(create_standata_outcome_singlestudy(study=rep(1, 10),
                                                        unit_id=1:10,
                                                        conc=rnorm(10),
                                                        case=rbinom(10, 1, 0.5)),
                                                        "standata_outcome")


})






test_that("Simple outcome model runs.", {
    out_standata <- create_standata_outcome_singlestudy(study=rep(1, 10),
                                        unit_id=1:10,
                                        conc=rnorm(10),
                                        case=rbinom(10, 1, 0.5))

    expect_error(create_standata_outcome_singlestudy(study=1:10,
                                                     unit_id=1:10,
                                                     conc=1:5), "The lengths")
})

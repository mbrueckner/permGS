context("util")

test_that("logrank_trafo is correct",
    time <- rexp(100, 1)
    status <- rbinom(n, 1, 0.6)
    
    x <- coin::logrank_trafo(Surv(time, status))
    y <- logrank_trafo(time, status)

    expect_identical(x, y)
}

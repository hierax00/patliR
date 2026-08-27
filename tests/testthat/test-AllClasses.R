test_that("patliR_project() creates a valid, empty project", {
  proj <- .test_project()
  expect_s4_class(proj, "PatliRProject")
  expect_equal(nrow(compounds(proj)), 0)
  expect_equal(nrow(matrixRaw(proj)), 0)
  expect_equal(nrow(binarizedMatrix(proj)), 0)
  expect_true(dir.exists(projectDir(proj)))
  expect_true(dir.exists(cacheDir(proj)))
})

test_that("validity is enforced: compounds must have the required columns", {
  proj <- .test_project()
  expect_error(
    compounds(proj) <- data.frame(id = "C0001"),
    class = "error"
  )
})

test_that("validity is enforced: matrix_raw and binarized must share columns", {
  proj <- .test_project()
  matrixRaw(proj) <- data.frame(compound_id = "C0001", LEA_ET = 100)
  expect_error(
    binarizedMatrix(proj) <- data.frame(compound_id = "C0001", FLO_ET = 1),
    class = "error"
  )
})

test_that("patliRResults() round-trips a named entry", {
  proj <- .test_project()
  expect_equal(patliRResults(proj), list())
  patliRResults(proj, "foo") <- data.frame(x = 1)
  expect_equal(patliRResults(proj, "foo"), data.frame(x = 1))
  expect_null(patliRResults(proj, "bar"))
})

test_that("show() prints a summary without erroring", {
  proj <- .test_project()
  expect_output(print(proj), "PatliRProject")
})

test_that(".cli_escape() doubles literal curly braces (glue's own escape convention)", {
  expect_equal(patliR:::.cli_escape("no braces here"), "no braces here")
  expect_equal(patliR:::.cli_escape("{a} and {b}"), "{{a}} and {{b}}")
  expect_equal(patliR:::.cli_escape('"class_results": [{"x": 1}]'), '"class_results": [{{"x": 1}}]')
})

test_that(".fetch_external() does not crash when the underlying error message contains curly braces", {
  ## Regression: a real NPClassifier failure (via compounds_classify(),
  ## fetch_mode = "warn_and_cache") raised an error whose message echoed
  ## back a raw JSON fragment (`"class_results": [...`) containing literal
  ## `{`/`}` -- cli::cli_warn() treats every bullet as a glue template, so
  ## that fragment made it try to parse() the JSON as R code and threw a
  ## confusing meta-error that fully masked the real "could not reach the
  ## resource" warning. See DEVLOG.md.
  cache_dir <- tempfile("patliR_test_cache_")
  fetch_fun <- function() stop('bad response: {"class_results": ["a", "b"]}')

  expect_warning(
    result <- patliR:::.fetch_external(fetch_fun, cache_dir, "brace_test", mode = "warn_and_cache"),
    "class_results"
  )
  expect_null(result)

  expect_error(
    patliR:::.fetch_external(fetch_fun, cache_dir, "brace_test", mode = "abort"),
    "class_results"
  )
})

test_that(".fetch_external() still surfaces a plain (brace-free) failure message unchanged", {
  cache_dir <- tempfile("patliR_test_cache_")
  fetch_fun <- function() stop("no ChEMBL molecule found for this SMILES (flexmatch)")

  expect_warning(
    result <- patliR:::.fetch_external(fetch_fun, cache_dir, "plain_test", mode = "warn_and_cache"),
    "no ChEMBL molecule found"
  )
  expect_null(result)
})

test_that(".match_column_flexible() prefers an exact match", {
  available <- c("Probability", "probability_2", "*Probability")
  expect_equal(patliR:::.match_column_flexible(available, "Probability"), "Probability")
})

test_that(".match_column_flexible() tolerates case/whitespace/decoration differences", {
  ## Real case: a SuperPred export used "*Probability" instead of
  ## "Probability" across runs.
  expect_equal(patliR:::.match_column_flexible(c("*Probability"), "Probability"), "*Probability")
  expect_equal(patliR:::.match_column_flexible(c(" probability "), "Probability"), " probability ")
  expect_equal(patliR:::.match_column_flexible(c("PROBABILITY"), "Probability"), "PROBABILITY")
  expect_equal(patliR:::.match_column_flexible(c("_Probability_"), "Probability"), "_Probability_")
})

test_that(".match_column_flexible() returns NA when nothing matches, even loosely", {
  expect_true(is.na(patliR:::.match_column_flexible(c("Foo", "Bar"), "Probability")))
})

test_that(".match_column_flexible() does not confuse a decorated column with an unrelated one", {
  ## "probability_2" strips to "probability2" which does not normalize to
  ## "probability" -- must not be mistaken for the target column.
  expect_true(is.na(patliR:::.match_column_flexible(c("probability_2"), "Probability")))
})

test_that(".parse_percent_column() leaves a plain 0..1 fraction column unchanged", {
  expect_equal(patliR:::.parse_percent_column(c("0.92", "0.5", "0.1")), c(0.92, 0.5, 0.1))
  expect_equal(patliR:::.parse_percent_column(c(0.92, 0.5)), c(0.92, 0.5))
})

test_that(".parse_percent_column() parses a '%' string column into a 0..1 fraction", {
  expect_equal(patliR:::.parse_percent_column(c("96.55%", "50%")), c(0.9655, 0.5))
})

test_that(".parse_percent_column() rescales a whole 0..100 column even without a literal '%'", {
  expect_equal(patliR:::.parse_percent_column(c("85", "50")), c(0.85, 0.5))
})

test_that(".parse_percent_column() decides the scale once per column, not per row", {
  ## One value has a literal '%' -> the WHOLE column is divided by 100,
  ## even the value(s) that look like they could already be a fraction.
  expect_equal(patliR:::.parse_percent_column(c("96.55%", "0.5")), c(0.9655, 0.005))
  ## No literal '%' anywhere, but the column max is > 1 -> same rule.
  expect_equal(patliR:::.parse_percent_column(c("85", "0.5")), c(0.85, 0.005))
})

test_that(".parse_percent_column() returns NA for unparseable entries without breaking the rest of the column", {
  result <- patliR:::.parse_percent_column(c("96.55%", "not_a_number", NA))
  expect_equal(result[1], 0.9655)
  expect_true(is.na(result[2]))
  expect_true(is.na(result[3]))
})

test_that(".read_csv_safe() reads a file missing its trailing newline without warning", {
  path <- tempfile(fileext = ".csv")
  ## writeChar with no trailing "\n" reproduces the "incomplete final line"
  ## condition real third-party exports (SuperPred, SwissADME, ...) trigger.
  writeChar("a,b\n1,2", path, eos = NULL)
  expect_no_warning(result <- patliR:::.read_csv_safe(path))
  expect_equal(result$a, 1)
  expect_equal(result$b, 2)
  unlink(path)
})

test_that(".read_csv_safe() still lets other, real warnings through", {
  ## Mock utils::read.csv() itself to raise a warning unrelated to
  ## "incomplete final line" -- proves .read_csv_safe() only muffles that
  ## one specific warning and lets everything else (malformed quoting,
  ## embedded nulls, ...) surface unchanged, as documented.
  testthat::local_mocked_bindings(
    read.csv = function(...) {
      warning("some other, unrelated read.csv() warning")
      data.frame(a = 1, b = 2)
    },
    .package = "utils"
  )
  path <- tempfile(fileext = ".csv")
  writeLines(c("a,b", "1,2"), path)
  expect_warning(patliR:::.read_csv_safe(path), "some other, unrelated")
  unlink(path)
})

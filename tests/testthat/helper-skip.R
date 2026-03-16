# Helper: skip_slow() — call at the start of any slow test.
# Run the full test suite with:  Sys.setenv(RSENA_FULL_TESTS = "1"); devtools::test()
# or from the shell:             RSENA_FULL_TESTS=1 Rscript -e "devtools::test()"
skip_slow <- function() {
  if (!identical(Sys.getenv("RSENA_FULL_TESTS"), "1")) {
    skip("Slow test. Set RSENA_FULL_TESTS=1 to run all tests.")
  }
}

context("ENT parsing")
load(test_path("T00032_truth.RData"))

T00032_test <-
  parse_ent(
    input_file = test_path("T00032.ent"),
    output_dir = ".",
    return_list = TRUE
  )



test_that("Input Assertions", {

  expect_error(parse_ent(input_file = "Bla.ent", output_dir = "."),
               "input_file not found.")

  expect_error(parse_ent(input_file = test_path("T00032.ent"),
                         output_dir = "./Bla"),
               "output_dir not found.")

})


test_that("T00032 Output File exists", {
  expect_true(file.exists(test_path("T00032.RData")))
})



test_that("T00032 output completely identical", {

  #test overall identical
  expect_identical(T00032_test, T00032_truth)

})

test_that("T00032 subunits identical", {

  #test subunits identical
  expect_identical(T00032_test$organism, T00032_truth$organism)
  expect_identical(T00032_test$protein, T00032_truth$protein)
  expect_identical(T00032_test$pathway, T00032_truth$pathway)
  expect_identical(T00032_test$enzyme, T00032_truth$enzyme)
  expect_identical(T00032_test$module, T00032_truth$module)
  expect_identical(T00032_test$db_links, T00032_truth$db_links)
  expect_identical(T00032_test$peptides, T00032_truth$peptides)

  #test subunit number of rows identical
  expect_equal(nrow(T00032_test$protein), nrow(T00032_truth$protein))
  expect_equal(nrow(T00032_test$pathway), nrow(T00032_truth$pathway))
  expect_equal(nrow(T00032_test$enzyme), nrow(T00032_truth$enzyme))
  expect_equal(nrow(T00032_test$module), nrow(T00032_truth$module))
  expect_equal(nrow(T00032_test$db_links), nrow(T00032_truth$db_links))
  expect_equal(nrow(T00032_test$peptides), nrow(T00032_truth$peptides))
  expect_equal(nrow(T00032_test$organism), nrow(T00032_truth$organism))


})



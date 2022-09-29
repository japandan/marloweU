test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

#checking inputs
#assertthat::assert_that(file.exists(input_file),
#                        msg = "input_file not found.")
#assertthat::assert_that(dir.exists(output_dir),
#                        msg = "output_dir not found.")
test_that("Input Assertions", {
  hydrogen_mass = "ACAT"
  input_file = "../../data-raw/uniprot/castor.fasta"
  output_dir = "../../data/rdata"

  expect_error(parse_fasta( input_file, output_dir, FALSE, hydrogen_mass ),
  "hydrogen_mass needs to be a numeric number")

  expect_error( parse_fasta( input_file="foo", output_dir, FALSE, hydrogen_mass ),
  "input_file not found")

  expect_error( parse_fasta( input_file, output_dir="foo", FALSE, hydrogen_mass ),
                "output_dir not found")

})


#' Parse a UniProt FASTA file
#'
#' @param input_file path to the input file
#' @param output_dir directory to write the output file, name will be the name
#'   of the input file with the .RData extension
#' @param hydrogen_mass mass of hydrogen to use when calculating neutral mass.
#'   Default is 1.00727646627 which comes from
#'   [NIST
#'   proton molar mass](https://physics.nist.gov/cgi-bin/cuu/Category?view=html&Atomic+and+nuclear.x=123&Atomic+and+nuclear.y=21)
#' @param return_list Should the parsed FASTA file be returned as well as written
#'   to a file? Defaults to FALSE.
#'
#' @return Default is to return an RData object to the specified location.
#'
#' @export
#'
#' @author Daniel T. Vogel
#'
#' @importFrom plyr adply ldply
#' @importFrom readr read_file
#' @importFrom OrgMassSpecR Digest
#' @importFrom stringr str_match str_trim str_replace_all str_extract_all str_match_all
#' @importFrom assertthat assert_that

# CRITICAL FIELDS TO EXTRACT FROM FASTA
#protein_id
#name
#organism
#aa-count
#aa-sequence
parse_fasta <- function(input_file,
                      output_dir,
                      return_list = FALSE,
                      hydrogen_mass = 1.00727646627){


  getwd()
  #checking inputs
  assertthat::assert_that(file.exists(input_file),
                          msg = "input_file not found.")
  assertthat::assert_that(dir.exists(output_dir),
                          msg = "output_dir not found.")
  assertthat::assert_that((is.numeric(hydrogen_mass)),
                          msg = "hydrogen_mass needs to be a numeric number")

}

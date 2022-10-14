#' Digest one protein from KEGG ENT file
#'
#' This function is meant to be called as part of parse_ent() which parses the
#' KEGG database ENT files to prepare them for upload to the MySQL database.
#'
#' @param protein_info a character vector of length 2 containing the protein_id
#'   and protein sequence
#' @param hydrogen_mass mass of hydrogen to use when calculating neutral mass.
#'   Default is 1.00727646627 which comes from
#'   [NIST
#'    proton molar mass](https://physics.nist.gov/cgi-bin/cuu/Category?view=html&Atomic+and+nuclear.x=123&Atomic+and+nuclear.y=21)
#'
#' @return Returns a dataframe with the protein ID, peptide string, and peptide
#'   mass.
#'
#' @author Sarah C. Jenson
#'
#' @export
#' @importFrom OrgMassSpecR Digest
#' @importFrom assertthat assert_that
digest_ent_protein <- function(protein_info, hydrogen_mass = 1.00727646627, digest_enzyme = "trypsin")
{
  # debug statements can be used to get metrics for performance
  # print( paste0("Digesting Protein=", protein_info[1], ", AA seq length=", width( protein_info[2] ),", protease=",digest_enzyme ))

  assertthat::assert_that(length(protein_info) == 2,
                          msg = "protein_info must have length 2")

  assertthat::assert_that(is.vector(protein_info),
                          msg = "protein_info needs to be a vector")

  assertthat::assert_that(is.character(protein_info[2]) & !is.na(protein_info[2]),
                msg = "aaseq in protein_info[2] must be non NA character string")

  assertthat::assert_that(is.character(protein_info[1]) & !is.na(protein_info[1]),
                msg = "protein_id in protein_info[1] must be a non-NA character string")


  #using suppressWarnings so it won't complain about proteins with no trypsin cutsites
  peps <-
    suppressWarnings(
      OrgMassSpecR::Digest(
        protein_info[2],
        enzyme = digest_enzyme,
        missed = 0,
        IAA = FALSE,
        custom = list(
          code = c("X", "B", "Z", "J", "U", "O"),
          mass = c(0,
                   0,
                   0,
                   131.094635,
                   168.964203,
                   255.158295)
        )
      )
    )

  #filtering for length
  peps <- peps[nchar(peps$peptide) > 4,]

 #handling proteins that don't have any > 4 peptides
  if (nrow(peps) == 0) {
    print(paste0(
      "Protin ID: ",
      protein_info[1],
      " had no tryptic peptides with length > 4."
    ))

    output <-
      data.frame(
        protein_id = protein_info[[1]],
        peptide = NA,
        mass = NA,
        stringsAsFactors = FALSE
      )
  } else
  {
    #converting M+1 mass to neutral mass and making dataframe
    output <-
      data.frame(
        protein_id = protein_info[[1]],
        peptide = peps$peptide,
        mass = (peps$mz1 - hydrogen_mass),
        stringsAsFactors = FALSE
      )

    #setting peptides with X, B, or Z to have a mass of -9
    output$mass[grepl("X|B|Z", output$peptide)] <- -9
  }

  return(output)
}

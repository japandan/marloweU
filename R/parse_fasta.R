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


library( Biostrings )
library( S4Vectors )

# CRITICAL FIELDS TO EXTRACT FROM FASTA
#protein_id
#name
#organism
#aa-count
#aa-sequence

# UniProtKB Fasta headers
#
# Example:
# >db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
#
# Actual Examples from a UniProt Fasta file.  Note that there are | delimiters, positional, and named column
# >sp|O65039|CYSEP_RICCO Vignain OS=Ricinus communis OX=3988 GN=CYSEP PE=1 SV=1
# >sp|B9RK42|GPC1_RICCO Glycerophosphocholine acyltransferase 1 OS=Ricinus communis OX=3988 GN=GPC1 PE=1 SV=1
# >sp|B9RU15|ATXR5_RICCO Probable Histone-lysine N-methyltransferase ATXR5 OS=Ricinus communis OX=3988 GN=ATXR5 PE=1 SV=1
#
# Where:
#
# db is 'sp' for UniProtKB/Swiss-Prot and 'tr' for UniProtKB/TrEMBL.
# | separator
# UniqueIdentifier is the primary accession number of the UniProtKB entry.
# | separator
# EntryName is the entry name of the UniProtKB entry.
# ProteinName is the recommended name of the UniProtKB entry as annotated in the RecName field. For UniProtKB/TrEMBL entries without a RecName field, the SubName field is used. In case of multiple SubNames, the first one is used. The 'precursor' attribute is excluded, 'Fragment' is included with the name if applicable.
# OrganismName is the scientific name of the organism of the UniProtKB entry.
# OrganismIdentifier is the unique identifier of the source organism, assigned by the NCBI.
# GeneName is the first gene name of the UniProtKB entry. If there is no gene name, OrderedLocusName or ORFname, the GN field is not listed.
# ProteinExistence is the numerical value describing the evidence for the existence of the protein.
# SequenceVersion is the version number of the sequence.

parse_fasta <- function(input_file,
                      output_dir=".",
                      return_list = FALSE,
                      hydrogen_mass = 1.00727646627){


  #checking inputs
  assertthat::assert_that(file.exists(input_file), msg = "input_file not found.")
  assertthat::assert_that(dir.exists(output_dir), msg = "output_dir not found.")
  assertthat::assert_that((is.numeric(hydrogen_mass)), msg = "hydrogen_mass needs to be a numeric number")

  # Function to read FASTA aa file
  fasta_list <- Biostrings::readAAStringSet(input_file)
  print(paste0(basename(input_file), " has ", length(fasta_list), " entries."))

  # Break apart the names row into variables...splits into columns db,UniqueIdentifier, and then multiple fields
  columns <-strsplit( names(fasta_list),"\\|")
  print( columns )

  UniqueIdentifier <- sapply( columns, "[", 2 )
  print(paste0( "UniqueIdenifier has ", length(UniqueIdentifier), " entries."))
  print( UniqueIdentifier)

  protein_info <-sapply( columns, "[", 3 )
  print(paste0( "protein_info has ", length(protein_info), " entries."))
  print( protein_info)
  #
  # we need to split protein info into EntryName<space>ProteinName<space>OS=xxx OX=##### GN=xxxxxPE=# SV=#
  # The ProteinName can contain spaces so this is a little tricky



  # temporary set the info to the fasta_list so something is returned
  organism_info <- fasta_list

  # Toggle the return of the organism_info based on parameter return_list
  if (return_list) {
    return(organism_info)
  } else{
    invisible(input_file)
  }


}



# test code for the function.  Convert to assetthat test later
fasta <- parse_fasta( "data-raw/uniprot/castor.head.fasta",return_list = TRUE)
names( fasta )
#str( fasta )

fasta <- parse_fasta( "data-raw/uniprot/castor.fasta",return_list = TRUE)


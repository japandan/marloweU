#' Parse a UniProt uniref FASTA file
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

# CRITICAL FIELDS TO EXTRACT FROM FASTA and their uniref field name
#protein_id. <-- UniqueIdentifier
#name <-- ClusterName
#organism <-- Tax=TaxonName
#aa-count
#aa-sequence

# UniRef fasta fields
# >UniqueIdentifier ClusterName n=Members Tax=TaxonName TaxID=TaxonIdentifier RepID=RepresentativeMember
#
# Where:
#
# UniqueIdentifier is the primary accession number of the UniRef cluster.
# ClusterName is the name of the UniRef cluster.
# Members is the number of UniRef cluster members.
# TaxonName is the scientific name of the lowest common taxon shared by all UniRef cluster members.
# TaxonIdentifier is the NCBI taxonomy identifier of the lowest common taxon shared by all UniRef cluster members.
# RepresentativeMember is the entry name of the representative member of the UniRef cluster.
# e.g.
#"UniRef50_A0A5A9P0L4 Peptidylprolyl isomerase n=1 Tax=Triplophysa tibetana TaxID=1572043 RepID=A0A5A9P0L4_9TELE"
#"UniRef50_A0A410P257 Glycogen synthase n=2 Tax=Candidatus Velamenicoccus archaeovorus TaxID=1930593 RepID=A0A410P257_9BACT"
#"UniRef50_A0A8J3NBY6 Uncharacterized protein n=2 Tax=Actinocatenispora rupis TaxID=519421 RepID=A0A8J3NBY6_9ACTN"
#"UniRef50_Q8WZ42 Titin n=2871 Tax=Vertebrata TaxID=7742 RepID=TITIN_HUMAN"
#"UniRef50_A0A401TRQ8 Uncharacterized protein (Fragment) n=2 Tax=Chiloscyllium TaxID=34767 RepID=A0A401TRQ8_CHIPU"
#"UniRef50_A0A6J2WDG0 titin n=196 Tax=cellular organisms TaxID=131567 RepID=A0A6J2WDG0_CHACN"

parse_fasta <- function(input_file,
                      output_dir=".",
                      return_list = FALSE,
                      hydrogen_mass = 1.00727646627){


  #checking inputs
  assertthat::assert_that( file.exists(input_file), msg = "input_file not found.")
  assertthat::assert_that( dir.exists(output_dir),  msg = "output_dir not found.")
  assertthat::assert_that( ( is.numeric(hydrogen_mass) ), msg = "hydrogen_mass needs to be a numeric number")

  # Function to read FASTA aa file
  fasta_list <- Biostrings::readAAStringSet( input_file )
  print(paste0("File ", basename(input_file), " has ", length( fasta_list), " entries."))

  #protein_id. <- UniqueIdentifier
  entry_id <- stringr::word( names(fasta_list), 1 )
  print( paste0("entry_id: ", entry_id) )

  #name <- ClusterName
  name <-  stringr::str_match( names(fasta_list), "\\s(.+)\\sn=")
  print( paste0("name: ",name[2]) )

  #organism <-- Tax=TaxonName
  organism <- stringr::str_match( names(fasta_list), "Tax=(.+)\\sTaxID=")
  print( paste0("organism: ",organism[2]))

  #aa_count
  aa_count <- width( fasta_list )
  #print( paste0("aa_count:", aa_count ))
  #aa_sequence
  aaseq <- toString( fasta_list )
  print( paste0("aaseq", aaseq ))


  #binding the results from each regular expression into one big matrix
  all_info <-
    cbind(
      entry_id,
      name[,2],
      organism[,2],
      aa_count,
      aaseq
    )

  #adding column names
  colnames(all_info) <-
    c(
      "protein_id",
      "name",
      "organism",
      "aa_count",
      "aaseq"
    )

  #rownames(all_info) <- 1:nrow(all_info)
  organism_info <- all_info

  # Toggle the return of the organism_info based on parameter return_list
  if (return_list) {
    return(organism_info)
  } else{
    invisible(input_file)
  }


}



# test code for the function.  Convert to assetthat test later
matrix <- parse_fasta( "data-raw/uniprot/castor.bean.protein.fasta", return_list = TRUE )
matrix

castor_matrix <- parse_fasta( "data-raw/uniprot/castor.bean.taxid.3988.uniref.fasta", return_list = TRUE)
dim(castor_matrix)
# show the first 5 entries, without the aaseq.
# we need to examine the aaseq returned to ensure multiple sequences are correct aa_count, etc
castor_matrix[1:5,1:4]






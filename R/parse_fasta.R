#' parse_fasta.R

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
#' @importFrom Biostring readAAStringSet
#
parse_fasta <- function(input_file,
                      output_dir=".",
                      return_list = FALSE,
                      hydrogen_mass = 1.00727646627){


  #checking inputs
  assertthat::assert_that( file.exists(input_file), msg = "input_file not found.")
  assertthat::assert_that( dir.exists(output_dir),  msg = "output_dir not found.")
  assertthat::assert_that( ( is.numeric(hydrogen_mass) ), msg = "hydrogen_mass needs to be a numeric number")


  #library( Biostrings )
  #library( S4Vectors )

  # CRITICAL FIELDS TO EXTRACT FROM FASTA and their uniref field name
  #
  # protein_id. <-- UniqueIdentifier
  # name <-- ClusterName
  # organism <-- Tax=TaxonName
  # aa-count <- from the AAString after reading with readAAStringSet
  # aaseq <- <- from the AAString after reading with readAAStringSet
  #
  # UniRef fasta fields
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
  # >UniqueIdentifier ClusterName n=Members Tax=TaxonName TaxID=TaxonIdentifier RepID=RepresentativeMember
  # >UniRef50_A0A5A9P0L4 Peptidylprolyl isomerase n=1 Tax=Triplophysa tibetana TaxID=1572043 RepID=A0A5A9P0L4_9TELE


  # Function to read FASTA AA file
  fasta_list <- Biostrings::readAAStringSet( input_file )
  print(paste0("File ", basename(input_file), " has ", length( fasta_list), " entries."))

  #protein_id. <- UniqueIdentifier
  entry_id <- stringr::word( names(fasta_list), 1 )
  #print( paste0("entry_id: ", entry_id) )

  #name <- ClusterName
  name <-  stringr::str_match( names(fasta_list), "\\s(.+)\\sn=")
  #print( paste0("name: ",name[2]) )

  #organism <-- Tax=TaxonName
  organism <- stringr::str_match( names(fasta_list), "Tax=(.+)\\sTaxID=")
  #print( paste0("organism: ",organism[2]))

  #aa_count
  aa_count <- width( fasta_list )#
  #print( paste0("aa_count:", aa_count ))
  #aa_sequence
  aaseq <- toString( fasta_list )
  #print( paste0("aaseq", aaseq ))


  # binding the results from each regular expression into one big matrix
  # Adding NA placeholders for missing info from FASTA in case we want to get it from another source
  all_info <-
    cbind(
      entry_id,
      NA, #cds
      name[,2],
      NA, #definition
      NA, #orthology
      organism[,2],
      NA, #pathway
      NA, #module
      NA, #brite
      NA, #position
      NA, #motif
      NA, #dblinks
      aa_count,
      aaseq
    )

# These are the items we are able to pull from the FASTA in case we want to create a smaller all_info matrix
#  colnames(all_info) <-
#    c(
#      "protein_id",
#      "name",
#      "organism",
#      "aa_count",
#      "aaseq"
#    )

  colnames(all_info) <-
    c(
      "protein_id",
      "cds",
      "name",
      "definition",
      "orthology",
      "organism",
      "pathway",
      "module",
      "brite",
      "position",
      "motif",
      "dblinks",
      "aa_count",
      "aaseq"
    )



  #removing unnecessary objects because they will be very large for large ENT files
  #The are already inserted into all_info
  rm(
    entry_id,
    name,
  #  organism,  # don't remove until we know it is correct
    aa_count,
    aaseq
  )



  #Making Output Dataframes#####

  #peptides####
  raw_seq <- all_info[, c("protein_id", "aaseq")]

  peptides <-
    plyr::adply(raw_seq,
                1,
                digest_ent_protein,
                .id = NULL)

  #removing proteins with no tryptic peptides of length > 4
  if (sum(is.na(peptides$peptide)) > 0) {
    no_pep_prots <- unique(peptides$protein_id[is.na(peptides$peptide)])

    all_info <-
      all_info[!(all_info[, "protein_id"] %in% no_pep_prots),]

    peptides <- peptides[!is.na(peptides$peptide), ]

    print(paste0(
      length(no_pep_prots),
      " proteins had no tryptic peptides with length > 4."
    ))
  }

  ##  uncomment as you check each section
  #
  # Protein
  #
  # These fields are not in the FASTA
  # definition <- NA
  # orthology <- NA
  # position <- NA
  # motif <- NA
  #
  protein <- as.data.frame(all_info[, c("protein_id",
                                         "name",
                                         "definition",
                                         "orthology",
                                         "position",
                                         "motif",
                                         "aaseq")], stringsAsFactors = FALSE)

  #
  # #removing extra space in orthology to match output from previous code
  # protein$orthology <- gsub("\\s\\s", " ", protein$orthology)
  #
  #
  # #Organism####
  #organism <-
  #   as.data.frame(unique(all_info[, c("organism", "cds")]),
  #                 stringsAsFactors = FALSE)
  #
  # #making sure there is only one organism
  # assertthat::assert_that(
  #   nrow(organism) == 1,
  #   msg = paste0(
  #     nrow(organism),
  #     " organisms were found. Only one distinct organism is allowed."
  #   )
  # )
  # #there must be at least 1 non-NA
  # assertthat::assert_that(!is.na(organism$organism[1]),
  #                         msg = "No organism fields detected.")
  #
  # #extracting info
  # organism_split <- stringr::str_extract_all(organism[1, 1], "\\S+")
  # organism$letter_code <- organism_split[[1]][1]
  # organism$genus <- organism_split[[1]][2]
  # organism$species <- organism_split[[1]][3]
  #
  # #removing the 3 letter code from the organism field to match output from
  # #previous code
  # organism$organism <-
  #   gsub("^[a-z]{3}\\s\\s", "", organism$organism)
  #
  #
  # colnames(organism)[2] <- "kegg_id"
  # rm(organism_split)
  #

  # Pathway
  # Since the FASTA has no Pathway entries, return an empty dataframe

  pathway <- data.frame(
                protein_id = NA,
                short_path = NA,
                description = NA)


  # Enzyme
  # Since FASTA contains no enzyme entries, return an empty dataframe

  enzyme <- data.frame(
                protein_id = NA,
                Enzyme = NA,
                description = NA)

  # Module
  # Since FASTA file contains no Module entries, return an empty dataframe

  module <- data.frame(
                protein_id = NA,
                module_code = NA,
                description = NA)


  #
  # db_links
  # Since FASTA has noDBlinks, return an empty dataframe

  db_links <- data.frame(
                protein_id = NA,
                database = NA,
                id = NA)




  # #Saving output####
  # output_name <- basename(input_file)
  # output_name <- gsub("ent", "RData", output_name)
  # output_path <-
  #   paste(normalizePath(output_dir), output_name, sep = "\\")
  #
  organism_info <- list(
     organism = organism,
     protein = protein,
     pathway = pathway,
     enzyme = enzyme,
     module = module,
     db_links = db_links,
     peptides = peptides
   )
  #
  # save(organism_info, file = output_path)
  #
  # print(paste0(basename(input_file), " contained ", nrow(protein), " proteins and ", nrow(peptides), " peptides."))
  # print(paste0("Parsing took ", difftime(Sys.time(), start_time, units = "secs"), " seconds."))
  # print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")



  # Toggle the return of the organism_info based on parameter return_list
  if (return_list) {
    return(organism_info)
  } else{
    invisible(input_file)
  }


}




  # test code for the function.  Convert to assetthat test later
#  matrix <- parse_fasta( "data-raw/uniprot/castor.bean.protein.fasta", return_list = TRUE )
#  matrix

  # This file has 14625 protein entries
#castor_matrix <- parse_fasta( "/nbacc/uniprot/castor.bean.taxid.3988.uniref.fasta", return_list = TRUE)
#dim(castor_matrix)
  # show the first 5 entries, without the aaseq.
  # we need to examine the aaseq returned to ensure multiple sequences are correct aa_count, etc
#castor_matrix[1:5,1:4]






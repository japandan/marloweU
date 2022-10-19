#' parse_uniref.R

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
#'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==
#' @export
#'
#' @author Daniel T. Vogel
#'
#' @importFrom plyr adply ldply
#' @importFrom readr read_file
#' @importFrom OrgMassSpecR Digest
#' @importFrom stringr str_match str_trim str_replace_all str_extract_all str_match_all
#' @importFrom assertthat assert_that
#' @importFrom Biostrings readAAStringSet
#
parse_uniref <- function(fasta_input_file,
                      output_dir=".",
                      return_list = FALSE,
                      hydrogen_mass = 1.00727646627){


  #checking inputs
  assertthat::assert_that( file.exists(fasta_input_file), msg = "input_file not found.")
  assertthat::assert_that( dir.exists(output_dir),  msg = "output_dir not found.")
  assertthat::assert_that( ( is.numeric(hydrogen_mass) ), msg = "hydrogen_mass needs to be a numeric number")

  #library( S4Vectors )

  # CRITICAL FIELDS TO EXTRACT FROM FASTA and their uniref field name
  #
  # protein_id. <-- UniqueIdentifier
  # name <-- ClusterName
  # organism <-- Tax=TaxonName
  # aa-count <- from the AAString after reading with readAAStringSet
  # aaseq    <- from the AAString after reading with readAAStringSet
  # kegg-id  <- TaxID
  #
  # UniRef fasta fields
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

  # start the timer to calculate how long this takes
  print(paste0("Starting file: ", basename(fasta_input_file), " at ", Sys.time()))
  start_time <- Sys.time()

  # Function to read FASTA AA file
  fasta_list <- Biostrings::readAAStringSet( fasta_input_file )

  no_entries <- length( fasta_list)
  print(paste0("File ", basename(fasta_input_file), " has ", no_entries , " entries."))

  #Use regular expresions to pull out the needed fields.  Need to be returned as.matrix
  #protein_id. <- UniqueIdentifier (first field, can contain underlines)
  #entry_id <- as.matrix( stringr::word( names(fasta_list), 1 ))
  entry_id <- stringr::str_match( names(fasta_list), "(?:^|(?:[.!?]\\s))(\\w+)")
  #print( paste0("entry_id: ", entry_id) )

  #name <- ClusterName
  name <-  stringr::str_match( names(fasta_list), "\\s(.+)\\sn=")
  #print( paste0("name: ",name[2]) )

  #organism <-- Tax=TaxonName
  # match the genus + species after Tax= and then the word after TaxID
  organism <- stringr::str_match( names(fasta_list), "Tax=(.+)\\sTaxID=(.+)\\s")
  print( paste0("organism: ",organism))

  #aa_count
  aa_count <- as.matrix( width( fasta_list ))
  #print( paste0("aa_count:", aa_count ))

  #aa_sequence
  # convert the sequences from AAString to conventional string
  names( fasta_list ) <- NULL
  aaseq <- as.matrix( sapply( fasta_list, toString  ))
  #print( paste0("aaseq", aaseq ))




  # binding the results from each regular expression into one big matrix
  # Adding NA placeholders for missing info from FASTA in case we want to get it from another source
  blanks <- matrix( data = NA, nrow = no_entries )
  all_info <-
    cbind(
      entry_id[,2],
      organism[,3], #TaxID replaces(kegg-id)
      name[,2],
      blanks, #definition
      blanks, #orthology
      organism[,2],
      blanks, #pathway
      blanks, #module
      blanks, #brite
      blanks, #position
      blanks, #motif
      blanks, #dblinks
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

# add row and column names
  rownames(all_info) <- 1:nrow(all_info)
  colnames(all_info) <-
    c(
      "protein_id",
      "taxid",
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
    organism,
    aa_count,
    aaseq
  )


  # Making Output Dataframes #####

  #peptides####
  raw_seq <- all_info[ , c("protein_id", "aaseq"), drop = FALSE]

  peptides <-
    plyr::adply(raw_seq,
                1,                   # iterate over rows
                digest_ent_protein,  # call the digest function for each row
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
  protein <- as.data.frame(all_info[,
                        c("protein_id",
                          "name",
                          "definition",
                          "orthology",
                          "position",
                          "motif",
                          "aaseq"), drop = FALSE] ,stringsAsFactors = FALSE)

  #
  # #removing extra space in orthology to match output from previous code
  # protein$orthology <- gsub("\\s\\s", " ", protein$orthology)
  #
  #
  # we need to drop all the organisms that are not matching the filename or taxid
  # #Organism####
  organism <-
     as.data.frame(unique(all_info[, c("organism","taxid"), drop = FALSE]),
                   stringsAsFactors = FALSE)
  #
  # #making sure there is only one organism
    assertthat::assert_that(
        nrow(organism) == 1,
        msg = paste0(
        nrow(organism)," organisms were found. Only one distinct organism is allowed."
      )
    )

  # #there must be at least 1 non-NA
    assertthat::assert_that(!is.na(organism$organism[1]),
                           msg = "No organism fields detected.")
  #
  # #extracting info
   organism_split <- stringr::str_extract_all(organism[1, 1], "\\S+")
   organism$letter_code <- NA
   organism$genus <- organism_split[[1]][1]
   organism$species <- organism_split[[1]][2]
   colnames(organism)[2] <- "kegg_id"
   rm(organism_split)
   print( organism )
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


  #Saving output####
  output_name <- basename(fasta_input_file)
  output_name <- gsub("fasta", "RData", output_name)
  output_path <- paste(normalizePath(output_dir), output_name, sep = "/")

  # output_path <-
  #   paste(normalizePath(output_dir), output_name, sep = "\\"). # windows version

  # building the organism data structure
  organism_info <- list(
     organism = organism,
     protein  = protein,
     pathway  = pathway,
     enzyme   = enzyme,
     module   = module,
     db_links = db_links,
     peptides = peptides
   )

  # uncomment when ready to create files
  save(organism_info, file = output_path)
  #
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print(paste0( "Output saved in: ", output_path ))
  print(paste0( basename(fasta_input_file), " contained ", nrow(protein), " proteins and ", nrow(peptides), " peptides." ))
  print(paste0( "Parsing took ", difftime(Sys.time(), start_time, units = "secs"), " seconds." ))
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

  # Toggle the return of the organism_info based on parameter return_list
  if (return_list) {
    return(organism_info)
  } else{
    invisible(input_file)
  }

}


#test_all<-function() {
      # test code for the function. Test edge cases, 1 entry multiple taxid in a file.
      output_dir <- "data/RData"

      # "castor.bean.protein.fasta contained 1 proteins and 30 peptides."
      fasta_input_file <-"data-raw/uniprot/uniref50_castor_bean.3988.head.fasta"
      castor1 <- parse_uniref( fasta_input_file,
                             output_dir,
                             return_list = TRUE )

      #print( str( castor1 ))

      # "uniref50_castor_bean.3988.head.fasta contained 1 proteins and 5 peptides."
      castor2 <- parse_uniref( "/nbacc/uniprot/uniref50_castor_bean.3988.fasta",
                                     output_dir,
                                     return_list = TRUE)
      print( str( castor2 ))
      #
      # # show the first 5 entries, without the aaseq.
      # # we need to examine the aaseq returned to ensure multiple sequences are correct aa_count, etc
      #
      # # truncated entry with 4 proteins and 2 taxid and 1 entry with no taxid
      # # "uniref50_chlamydia_pneumoniae.head.fasta contained 4 proteins and 38 peptides."
      # fasta_input_file <- "data-raw/uniprot/uniref50_chlamydia_pneumoniae.head.fasta"
      # clap1 <- parse_uniref( fasta_input_file,
      #                      output_dir,
      #                      return_list = TRUE )
      #
      # print( str( clap1 ))
      #
      # # uniref50_chlamydia_pneumoniae.fasta contained 872 proteins and 14513 peptides.
      # "uniref50_chlamydia_pneumoniae.fasta contained 872 proteins and 14513 peptides."
      # clap2 <- parse_uniref( "data-raw/uniprot/uniref50_chlamydia_pneumoniae.fasta",
      #                      output_dir,
      #                      return_list = TRUE )
      #
      # print( str( clap2 ))
#}

#test_all()








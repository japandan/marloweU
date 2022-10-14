#' parse_ent.R
#'
#' Parse a KEGG ENT file
#'
#' @param input_file path to the input file
#' @param output_dir directory to write the output file, name will be the name
#'   of the input file with the .RData extension
#' @param hydrogen_mass mass of hydrogen to use when calculating neutral mass.
#'   Default is 1.00727646627 which comes from
#'   [NIST
#'   proton molar mass](https://physics.nist.gov/cgi-bin/cuu/Category?view=html&Atomic+and+nuclear.x=123&Atomic+and+nuclear.y=21)
#' @param return_list Should the parsed ENT file be returned as well as written
#'   to a file? Defaults to FALSE.
#'
#' @return Default is to return inpuwrites an RData object to
#'   the specified location.
#' @export
#'
#' @author Sarah C. Jenson
#'
#' @importFrom plyr adply ldply
#' @importFrom readr read_file
#' @importFrom OrgMassSpecR Digest
#' @importFrom stringr str_match str_trim str_replace_all str_extract_all str_match_all
#' @importFrom assertthat assert_that
parse_ent <- function(input_file,
                      output_dir,
                      return_list = FALSE,
                      hydrogen_mass = 1.00727646627){


  #checking inputs
  assertthat::assert_that(file.exists(input_file),
                          msg = "input_file not found.")
  assertthat::assert_that(dir.exists(output_dir),
                          msg = "output_dir not found.")
  assertthat::assert_that((is.numeric(hydrogen_mass)),
                          msg = "hydrogen_mass needs to be a numeric number")

  print(paste0("Starting file: ", basename(input_file), " at ", Sys.time()))

    start_time <- Sys.time()

  raw <- read_file(input_file)

  entries <- unlist(strsplit(raw, "\\r?\\n\\/\\/\\/\\r?\\n"))

  rm(raw)

  print(paste0(basename(input_file), " has ", length(entries), " entries."))

  #Extracting raw values from ent file using regular expressions####
  #str_match is vectorized so it will look for the regular expression in every
  #entry in the list and return a character matrix, one row for each entry in
  #the original list with a column for every capturing group in the regular
  #expression
  entry_id <- stringr::str_match(entries, "(?:^|\\s)ENTRY\\s+(\\S+)\\s+")

  cds <- stringr::str_match(entries, "(?:^|\\s)CDS\\s+(\\S+)\\s+")
  name <-
    stringr::str_match(entries, "(?:^|\\s)NAME\\s+((?:.|\\n|\\r)+?)(?=\\r?\\n\\S)")

  definition <-
    stringr::str_match(entries, "(?:^|\\s)DEFINITION\\s+((?:.|\\n|\\r)+?)(?=\\r?\\n\\S)")

  orthology <-
    stringr::str_match(entries, "(?:^|\\s)ORTHOLOGY\\s+((?:.|\\n|\\r)+?)(?=\\r?\\n\\S)")

  #only matches until the first newline
  organism <-
    stringr::str_match(entries, "(?:^|\\s)ORGANISM\\s+(.*)(?=\\r?\\n)")

  pathway <-
    stringr::str_match(entries, "(?:^|\\s)PATHWAY\\s+((?:.|\\n|\\r)+?)(?=\\r?\\n\\S)")

  module <-
    stringr::str_match(entries, "(?:^|\\s)MODULE\\s+((?:.|\\n|\\r)+?)(?=\\r?\\n\\S)")

  brite <-
    stringr::str_match(entries, "(?:^|\\s)BRITE\\s+((?:.|\\n|\\r)+?)(?=\\r?\\n\\S)")

  position <-
    stringr::str_match(entries, "(?:^|\\s)POSITION\\s+((?:.|\\n|\\r)+?)(?=\\r?\\n\\S)")

  motif <-
    stringr::str_match(entries, "(?:^|\\s)MOTIF\\s+((?:.|\\n|\\r)+?)(?=\\r?\\n\\S)")

  dblinks <-
    stringr::str_match(entries, "(?:^|\\s)DBLINKS\\s+((?:.|\\n|\\r)+?)(?=\\r?\\n\\S)")

  aa_info <-
    stringr::str_match(entries,
                       "(?:^|\\s)AASEQ\\s+(\\d+)\\r?\\n((?:.|\\n|\\r)+?)(?=\\r?\\n\\S)")

  aa_info[, 2] <- stringr::str_trim(aa_info[, 2])
  aa_info[, 3] <- stringr::str_replace_all(aa_info[, 3], "\\s", "")



  #binding the results from each regular expression into one big matrix
  all_info <-
    cbind(
      entry_id[, 2],
      cds[, 2],
      name[, 2],
      definition[, 2],
      orthology[, 2],
      organism[, 2],
      pathway[, 2],
      module[, 2],
      brite[, 2],
      position[, 2],
      motif[, 2],
      dblinks[, 2],
      aa_info[, 2:3]
    )

  #adding column names
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


  rownames(all_info) <- 1:nrow(all_info)

  #removing rows with no amino acid sequence
  all_info <- all_info[!is.na(all_info[, "aaseq"]), ]

  #making sure all entries have a protein_id

  #if protein_id is missing use protein name
  na_pos <- is.na(all_info[, "protein_id"])
  if (sum(na_pos) > 0) {
    all_info[na_pos, "protein_id"] <- all_info[na_pos, "name"]
    print(paste0(
      sum(na_pos),
      " missing entry fields were replaced with the entry name."
    ))
  }


  #if protein_id is still NA use filename followed by entry position
  na_pos <- is.na(all_info[, "protein_id"])
  if (sum(na_pos) > 0) {
    input_name <- basename(input_file)
    input_name <- gsub(".ent", "", input_name, ignore.case = TRUE)

    all_info[na_pos, "protein_id"] <-
      paste0(input_name,
             ": ",
             rownames(all_info)[na_pos])

    print(
      paste0(
        sum(na_pos),
        " missing entry fields were replaced with the Kegg ID followed by the entry position."
      )
    )
  }


  #removing unnecessary objects because they will be very large for large ENT files
  rm(
    entry_id,
    cds,
    name,
    definition,
    orthology,
    organism,
    pathway,
    module,
    brite,
    position,
    motif,
    dblinks,
    aa_info,
    na_pos
  )

  print( str( all_info ) )
  print(paste0(colnames(all_info) ))
  print(paste0( all_info[1]))

  #Making Output Dataframes#####

  #peptides####
  raw_seq <- all_info[, c("protein_id", "aaseq")]
  #print( raw_seq )
  print(paste0("raw_seq has ", length( raw_seq ), " length."))
  print(paste0("raw_seq is a vector ", is.vector( raw_seq )))

  #print( raw_seq[1] )
  print(paste0("raw_seq[1] has ", length( raw_seq[1] ), " length."))
  print(paste0("raw_seq[1] is a vector ", is.vector( raw_seq[1] )))

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



  #Protein####
  protein <- as.data.frame(all_info[, c("protein_id",
                               "name",
                               "definition",
                               "orthology",
                               "position",
                               "motif",
                               "aaseq")], stringsAsFactors = FALSE)


  #removing extra space in orthology to match output from previous code
  protein$orthology <- gsub("\\s\\s", " ", protein$orthology)


  #Organism####
  organism <-
    as.data.frame(unique(all_info[, c("organism", "cds")]),
                  stringsAsFactors = FALSE)

  #making sure there is only one organism
  assertthat::assert_that(
    nrow(organism) == 1,
    msg = paste0(
      nrow(organism),
      " organisms were found. Only one distinct organism is allowed."
    )
  )
  #there must be at least 1 non-NA
  assertthat::assert_that(!is.na(organism$organism[1]),
                          msg = "No organism fields detected.")

  #extracting info
  organism_split <- stringr::str_extract_all(organism[1, 1], "\\S+")
  organism$letter_code <- organism_split[[1]][1]
  organism$genus <- organism_split[[1]][2]
  organism$species <- organism_split[[1]][3]

  #removing the 3 letter code from the organism field to match output from
  #previous code
  organism$organism <-
    gsub("^[a-z]{3}\\s\\s", "", organism$organism)


  colnames(organism)[2] <- "kegg_id"
  rm(organism_split)


  #Pathway####

  #if the file has no Pathway entries return an empty dataframe
  if (sum(!is.na(all_info[, 7])) == 0) {
    pathway <-
      data.frame(protein_id = NA,
                 short_path = NA,
                 description = NA)

    print(paste0(basename(input_file), " contained no pathway information."))
  } else
  {
    pathway_raw <-
      unique(all_info[!is.na(all_info[, 7]), c("protein_id", "pathway")])

    #using regular expression to extract every line in the pathway section as a
    #separate row and the short path and description as separate columns
    pathway_split <-
      stringr::str_match_all(pathway_raw[, 2], "(\\S+)\\s+(.+)(?=(?:\\r?\\n)|$)")

    names(pathway_split) <- pathway_raw[, 1]

    #converting the list of matrices into a dataframe
    pathway <-
      plyr::ldply(pathway_split, function(x) {
        as.data.frame(x, stringsAsFactors = FALSE)
      })

    pathway$V1 <- NULL

    colnames(pathway) <-
      c("protein_id", "short_path", "description")

    rm(pathway_split, pathway_raw)
  }


  #Enzyme(Brite)####

  #if the file contains no enzyme entries then return an empty dataframe
  if (sum(!is.na(all_info[, 9])) == 0) {
    enzyme <- data.frame(protein_id = NA,
                         Enzyme = NA,
                         description = NA)
    print(paste0(basename(input_file), " contained no enzyme information."))

  } else{
    enzyme_raw <-
      unique(all_info[!is.na(all_info[, 9]), c("protein_id", "brite")])

    enzyme_split <-
      str_match_all(enzyme_raw[, 2], "(\\d\\.(?:\\d+|-)\\.(?:\\d+|-)\\.(?:\\d+|-))\\s+?(.*)\\r?\\n")

    names(enzyme_split) <- enzyme_raw[, 1]

    enzyme <-
      ldply(enzyme_split, function(x) {
        as.data.frame(x, stringsAsFactors = FALSE)
      })

    enzyme$V1 <- NULL

    colnames(enzyme) <- c("protein_id", "Enzyme", "description")

    #changing empty strings to NA
    enzyme$description[!stringr::str_detect(enzyme$description, "\\S")] <- NA

    rm(enzyme_raw, enzyme_split)

  }

  #Module####

  #if the file contains no Module entries return an empty dataframe
  if (sum(!is.na(all_info[, 8])) == 0) {
    module <-
      data.frame(
        protein_id = NA,
        module_code = NA,
        description = NA
      )
    print(paste0(basename(input_file), " contained no module information."))

  } else
  {
    module_raw <-
      unique(all_info[!is.na(all_info[, 8]), c("protein_id", "module")])

    module_split <-
      stringr::str_match_all(module_raw[, 2], "(\\S+)\\s+(.+)(?=(?:\\r?\\n)|$)")

    names(module_split) <- module_raw[, 1]

    module <-
      plyr::ldply(module_split, function(x) {
        as.data.frame(x, stringsAsFactors = FALSE)
      })

    module$V1 <- NULL

    colnames(module) <-
      c("protein_id", "module_code", "description")

    rm(module_raw, module_split)
  }

  #db_links####
  #if a file contains no DBlinks return an empty dataframe
  if (sum(!is.na(all_info[, 12])) == 0) {
    db_links <- data.frame(protein_id = NA,
                           database = NA,
                           id = NA)
    print(paste0(basename(input_file), " contained no db link information."))

  } else{
    db_links_raw <-
      unique(all_info[!is.na(all_info[, 12]), c("protein_id", "dblinks")])


    #old regex which gets all accessions "(\\S+):\\s+(.+)(?=(?:\\r?\\n)|$)"
    #new regex only gets first accession
    db_links_split <-
      stringr::str_match_all(db_links_raw[, 2],
                             "(\\S+):\\s+(\\S+)")

    names(db_links_split) <- db_links_raw[, 1]

    db_links <-
      plyr::ldply(db_links_split, function(x) {
        as.data.frame(x, stringsAsFactors = FALSE)
      })

    db_links$V1 <- NULL

    colnames(db_links) <-
      c("protein_id", "database", "id")

    rm(db_links_split, db_links_raw)
  }





  #Saving output####
  output_name <- basename(input_file)
  output_name <- gsub("ent", "RData", output_name)
  output_path <-
    paste(normalizePath(output_dir), output_name, sep = "\\")

  organism_info <- list(
    organism = organism,
    protein = protein,
    pathway = pathway,
    enzyme = enzyme,
    module = module,
    db_links = db_links,
    peptides = peptides
  )

  save(organism_info, file = output_path)

  print(paste0(basename(input_file), " contained ", nrow(protein), " proteins and ", nrow(peptides), " peptides."))
  print(paste0("Parsing took ", difftime(Sys.time(), start_time, units = "secs"), " seconds."))
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

  if (return_list) {
    return(organism_info)
  } else{
    invisible(input_file)
  }
}

# test file
#parse_ent( "data-raw/T00032.ent",".")

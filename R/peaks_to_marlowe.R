#' Generate tags from PEAKS
#'
#' @description
#' This functions will take in a dataframe read in from a de novo peptides.csv
#' file exported from Peaks De Novo CMD.
#' `r lifecycle::badge("maturing")`
#'
#' @param all_peaks_denovo dataframe of peaks de novo results
#' @param LC_threshold LC threshold used to define tags see \code{\link{generate_tags}}
#' @param prolineBlock should P after R or K prevent splitting at that site
#' @param output_dir desired output directory
#' @param min_length minimum tag length
#' @param ppm_error +/- allowed mass error in ppm, defaults to 15
#' @param return_tags Default `TRUE`, Should the tags dataframe be returned? If
#'   FALSE then just the path to the tags dataframe is returned.
#' @param new_mods (Default NULL) A named vector containing the masses of any
#'   additional modifications other then methionine oxidation or cysteine
#'   carbamidomethylation/alkylation. Vector names must be the rounded masses in
#'   quotes that appear within the `()` in the `Peptide` column. Ex:
#'   ABC(+57.02)HK => c("57.02" = 57.021464). It is recommended you get the
#'   masses from unimod.org which is the source for the default PTMs methionine
#'   oxidation (mass: 15.994915) and cysteine carbamidomethylation/alkylation
#'   (mass: 57.021464)).
#' @param hydrogen_mass mass of hydrogen, default 1.00727646627 comes from
#'  [NIST
#'   proton molar mass](https://physics.nist.gov/cgi-bin/cuu/Category?view=html&Atomic+and+nuclear.x=123&Atomic+and+nuclear.y=21)
#'
#' @return If `return_tags = TRUE` a dataframe of de novo tags is returned as
#'   well as saved to the specified `output_dir`. If `return_tags = FALSE` then
#'   only the path to the output RData file is returned.
#' @export
#'
#' @author Isabelle O'Bryon, Sarah C. Jenson
#'
#'
#' @examples \dontrun{\donttest{
#' # See vignette "Generate De novo Tags" using
#' # `browseVignettes("CandidateSearchDatabase")`}}
#'
#' @importFrom plyr ddply
#' @importFrom assertthat assert_that
generate_PEAKS_tags <- function(all_peaks_denovo,
                                     output_dir,
                                     LC_threshold,
                                     prolineBlock,
                                     min_length,
                                     ppm_error = 15,
                                     return_tags = TRUE,
                                     new_mods = NULL,
                                     hydrogen_mass = 1.00727646627) {
  #checking inputs
  assertthat::assert_that(dir.exists(output_dir),
                          msg = "output_dir not found.")

  assertthat::assert_that("PTM" %in% colnames(all_peaks_denovo),
                          msg = cat("Column `PTM` not found in all_peaks_denovo which has columns:",
                                    colnames(all_peaks_denovo)))

  assertthat::assert_that("Scan" %in% colnames(all_peaks_denovo),
                          msg = cat("Column `Scan` not found in all_peaks_denovo which has columns:",
                                    colnames(all_peaks_denovo)))


  assertthat::assert_that("local.confidence...." %in% colnames(all_peaks_denovo),
                          msg = cat("Column `local.confidence....` not found in all_peaks_denovo which has columns:",
                                    colnames(all_peaks_denovo)))

  assertthat::assert_that("Peptide" %in% colnames(all_peaks_denovo),
                          msg = cat("Column `Peptide` not found in all_peaks_denovo which has columns:",
                                    colnames(all_peaks_denovo)))


  assertthat::assert_that("Source.File" %in% colnames(all_peaks_denovo),
                          msg = cat("Column `Source.File` not found in all_peaks_denovo which has columns:",
                                    colnames(all_peaks_denovo)))

  assertthat::assert_that(length(unique(all_peaks_denovo$Source.File)) == 1,
                          msg = "all_peaks_denovo$Source.File contains more than one distinct file.")

  ## Call create_tags to make a dataframe of all the tags with ALC above the
    ## LC_threshold
    all_peaks_denovo_tag <-
      plyr::ddply(
        all_peaks_denovo,
        "Scan",
        create_peaks_tags,
        min_conf = LC_threshold,
        min_length = min_length,
        pro_block = prolineBlock
      )

    ## Convert all_peaks_denovo_tag into a dataframe with mass windows that can
    ## be used as input to MARLOWE database query
    all_peaks_tags <-
      plyr::adply(all_peaks_denovo_tag,
                  1,
                  calc_mass_window,
                  LC_threshold = LC_threshold,
                  ppm_error = ppm_error,
                  new_mods = new_mods,
                  hydrogen_mass = hydrogen_mass,
                  .expand = FALSE,
                  .id = NULL)


    ## Add Weight and sample_tag_id columns

    all_peaks_tags$Weight <- 1
    all_peaks_tags$sample_tag_id <- 1:nrow(all_peaks_tags)

    #Making output filename and saving results
    file_name <- unique(all_peaks_denovo$Source.File)
    file_nameclean <- gsub(".raw", "", file_name, ignore.case = T)
    file_nameclean <- gsub(".mzML", "", file_nameclean, ignore.case = T)
    file_nameclean <- gsub(".mgf", "", file_nameclean, ignore.case = T)
    file_nameclean <- paste0(file_nameclean, "_Peaks_tags.RData")


    output_filepath <- file.path(normalizePath(output_dir),
                                 file_nameclean, fsep = "/")

    print(c("Saving to ",output_filepath ))
    save(all_peaks_tags, file = output_filepath)

    if(return_tags == TRUE){
      return(all_peaks_tags)
    } else {
      return(output_filepath)
    }

  }


#' Calculate mass window for PEAKS results
#'
#' This function takes in a dataframe of Peaks De Novo tags and
#' scan information and converts it into the format that Marlowe can use.
#' `r lifecycle::badge("maturing")`
#'
#' @param peaks_tag_df Dataframe of Peaks tags
#' @param LC_threshold LC threshold used to define tags
#' @param ppm_error +/- allowed mass error in ppm, defaults to 15
#' @param new_mods (Default NULL) A named vector containing the masses of any
#'   additional modifications other then methionine oxidation or cysteine
#'   carbamidomethylation/alkylation. Vector names must be the rounded masses in
#'   quotes that appear within the `()` in the `Peptide` column. Ex:
#'   ABC(+57.02)HK => c("57.02" = 57.021464). It is recommended you get the
#'   masses from unimod.org which is the source for the default PTMs methionine
#'   oxidation (mass: 15.994915) and cysteine carbamidomethylation/alkylation
#'   (mass: 57.021464)).
#' @param hydrogen_mass mass of hydrogen, default 1.00727646627 comes from
#'  [NIST
#'   proton molar mass](https://physics.nist.gov/cgi-bin/cuu/Category?view=html&Atomic+and+nuclear.x=123&Atomic+and+nuclear.y=21)
#'
#' @return Returns a dataframe in MARLOWE format
#'
#' @author Isabelle O'Bryon, Sarah C. Jenson
#'
#' @export
#'
#' @importFrom assertthat assert_that
#' @importFrom stringr str_detect str_extract_all
calc_mass_window <- function(peaks_tag_df,
                            LC_threshold,
                            ppm_error = 15,
                            new_mods = NULL,
                            hydrogen_mass = 1.00727646627) {
  #checking inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  assertthat::assert_that("tryptic_tag" %in% colnames(peaks_tag_df),
                          msg = "Column `tryptic_tag` not found in peaks_tag_df")

  assertthat::assert_that("Peptide" %in% colnames(peaks_tag_df),
                          msg = "Column `Peptide` not found in peaks_tag_df")

  assertthat::assert_that("m.z" %in% colnames(peaks_tag_df),
                          msg = "Column `m.z` not found in peaks_tag_df")

  assertthat::assert_that("z" %in% colnames(peaks_tag_df),
                          msg = "Column `z` not found in peaks_tag_df")

  assertthat::assert_that("Scan" %in% colnames(peaks_tag_df),
                          msg = "Column `Scan` not found in peaks_tag_df")

  assertthat::assert_that(length(unique(peaks_tag_df$Scan)) == 1,
                          msg = "peaks_tag_df$Scan contains more than one distinct value.")



  if (is.null(new_mods) == FALSE) {
    assertthat::assert_that(
      (is.vector(new_mods) &
         is.numeric(new_mods)),
      msg = paste0(
        "new_mods is a ",
        class(new_mods),
        " but must be a numeric vector"
      )
    )
    assertthat::assert_that(is.null(names(new_mods)) == FALSE,
                            msg = "new_mods must be a named vector")
  }

  #adjusting mass for modifications~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #(?<=\()([\+|\-]\d+\.\d+)+(?=\))

  if(peaks_tag_df$PTM != ""){
    # Default named vector of modification names (as they appear in PEAKS results) and
    # their unimod masses
    mod_vector <- c("15.99" = 15.994915, "57.02" = 57.021464)

    #adding new_mods to default mods if present
    if(is.null(new_mods) == FALSE){
      mod_vector <- c(mod_vector, new_mods)
    }

    #extracting short modification masses from Peptide column
    short_mod_mass <- unlist(stringr::str_extract_all(peaks_tag_df$Peptide,
                                               pattern = "(?<=\\()([\\+|\\-]\\d+\\.\\d+)+(?=\\))"))

    #removing + sign
    short_mod_mass <- as.character(as.numeric(short_mod_mass))

    mod_mass <- sum(mod_vector[short_mod_mass], na.rm = TRUE)

    #calculating neutral mass
    peptide_total_mass <- calc_neutral_mass(mz = peaks_tag_df$m.z,
                                            charge = peaks_tag_df$z,
                                            hydrogen_mass =  hydrogen_mass)

    #adjusting neutral mass
    peptide_total_mass <- peptide_total_mass - mod_mass

  } else {
    #calculating neutral mass
    peptide_total_mass <- calc_neutral_mass(mz = peaks_tag_df$m.z,
                                            charge = peaks_tag_df$z,
                                            hydrogen_mass =  hydrogen_mass)
  }

  #calculating min and max mass in Da corresponding to ppm_min and ppm_max
  min_mass <- ppm_to_da(-ppm_error, peptide_total_mass)
  max_mass <- ppm_to_da(ppm_error, peptide_total_mass)

  original_peptide <- gsub("[^A-Z]", "", peaks_tag_df$Peptide)

  output_df <- data.frame(
    tag = peaks_tag_df$tryptic_tag,
    pep_mass = peptide_total_mass,
    min_mass = min_mass,
    max_mass = max_mass,
    min_mass_ppm = -ppm_error,
    max_mass_ppm = ppm_error,
    original_peptide = original_peptide,
    scan_num = peaks_tag_df$Scan,
    LC_threshold = LC_threshold,
    stringsAsFactors = FALSE
  )

  return(output_df)
}


#' create_peaks_tags
#'
#' Generate high confidence tags from PEAKS denovo results.
#' `r lifecycle::badge("maturing")`
#'
#' @param peaks_denovo_df A dataframe of PEAKS results corresponding to one scan.
#' @param min_conf Minimum confidence score
#' @param min_length Minimum length
#' @param pro_block Should P after R or K prevent digestion at that site?
#'
#' @return Returns a dataframe of peaks tags
#'
#' @export
#'
#' @author Isabelle O'Bryon, Sarah C. Jenson
#' @importFrom assertthat assert_that
create_peaks_tags <- function(peaks_denovo_df,
                              min_conf,
                              min_length,
                              pro_block){

  #checking inputs
  assertthat::assert_that("local.confidence...." %in% colnames(peaks_denovo_df),
                          msg = "Column `local.confidence....` not found in input dataframe.")

  assertthat::assert_that("Peptide" %in% colnames(peaks_denovo_df),
                          msg = "Column `Peptide` not found in input dataframe.")

  assertthat::assert_that(length(unique(peaks_denovo_df$Scan)) == 1,
                          msg = "`Scan` column contains more than one distinct value.")



  local.confidence <-
    as.numeric(unlist(strsplit(
      peaks_denovo_df$local.confidence...., split = " "
    )))

  peptide <- gsub("[^A-Z]", "", peaks_denovo_df$Peptide)
  conf_tag <- ""

  for (j in 1:length(local.confidence)) {
    aa <- substr(peptide, j, j)
    if (local.confidence[j] > min_conf) {
      conf_tag <- paste(conf_tag, aa, sep = "")
    } else{
      conf_tag <- paste(conf_tag, ".", sep = "")
    }
  }

  tag_list <- unlist(strsplit(conf_tag, ".", fixed = T))

  # removing empty and short strings
  tag_list <- tag_list[tag_list != ""]

  # if tags is empty return NULL
  if (length(tag_list) == 0) {
    return(NULL)
  }

  #digesting tags
  tryp_tag_list <- MakeSearchSim::digestPeptides(
    peptides = tag_list,
    prolineBlock = pro_block,
    discardLessThan = min_length
  )

  # if tags is empty return NULL
  if (length(tryp_tag_list) == 0) {
    return(NULL)
  }

  #duplicating input row if more than one tag in PSM
  tag_peaks_denovo_df <- peaks_denovo_df[rep(1, length(tryp_tag_list)),]
  tag_peaks_denovo_df$tryptic_tag <- tryp_tag_list

  return(tag_peaks_denovo_df)
}

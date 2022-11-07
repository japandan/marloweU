


#' PEAKS tags from multiple source files
#'
#'
#' @param PEAKS_output_csv_path Path to the de novo peptides.csv file generated
#'   by PEAKS which contains the de novo results.
#' @param output_dir Relative or absolute path to the directory where tags
#'   should be written as RData files. Defaults to current working directory.
#' @param min_ALC Minimum ALC (average local confidence) required for the PSM to
#'   be considered for tag generation. Default 50.
#' @param min_peptide_length Minimum length of of the overall peptide required
#'   for the PSM to be considered for tag generation. Default 6.
#' @param min_tag_length Minimum tag length. Default 5.
#'
#' @inheritParams generate_PEAKS_tags
#'
#' @return If `return_tags = TRUE` then returns a list of dataframes (or a
#'   single dataframe of only one Source.File), each entry being for each
#'   distinct value found within the Source.File column of the PEAKS output
#'   file. If `return_tags = FALSE` then a list of paths to RData files holding
#'   the tags will be returned instead.
#' @export
#'
#' @examples \dontrun{\donttest{
#' # See vignette "Generate De novo Tags" using
#' # `browseVignettes("CandidateSearchDatabase")`}}
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#' @importFrom assertthat assert_that
#' @importFrom purrr map
#' @importFrom rlang .data
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
batch_PEAKS_tags <- function(PEAKS_output_csv_path,
                             output_dir = ".",
                             min_ALC = 50,
                             LC_threshold = 80,
                             min_peptide_length = 6,
                             min_tag_length = 5,
                             ppm_error = 15,
                             return_tags = TRUE,
                             new_mods = NULL,
                             prolineBlock = TRUE,
                             hydrogen_mass = 1.00727646627) {



  #checking inputs
  assertthat::assert_that(file.exists(PEAKS_output_csv_path),
                          msg = "PEAKS_output_csv_path not found.")


  assertthat::assert_that(dir.exists(output_dir),
                          msg = "output_dir not found.")


  #loading PEAKS results
  all_peaks <- read.csv(PEAKS_output_csv_path,
                        header = TRUE,
                        stringsAsFactors = FALSE)


  #checking column names
  assertthat::assert_that(
    "PTM" %in% colnames(all_peaks),
    msg = cat(
      "Column `PTM` not found in input csv which has columns:",
      colnames(all_peaks)
    )
  )

  assertthat::assert_that(
    "Scan" %in% colnames(all_peaks),
    msg = cat(
      "Column `Scan` not found in input csv which has columns:",
      colnames(all_peaks)
    )
  )


  assertthat::assert_that(
    "local.confidence...." %in% colnames(all_peaks),
    msg = cat(
      "Column `local.confidence....` not found in input csv which has columns:",
      colnames(all_peaks)
    )
  )

  assertthat::assert_that(
    "Peptide" %in% colnames(all_peaks),
    msg = cat(
      "Column `Peptide` not found in input csv which has columns:",
      colnames(all_peaks)
    )
  )


  assertthat::assert_that(
    "Source.File" %in% colnames(all_peaks),
    msg = cat(
      "Column `Source.File` not found in input csv which has columns:",
      colnames(all_peaks)
    )
  )

  assertthat::assert_that(
    "m.z" %in% colnames(all_peaks),
    msg = cat(
      "Column `m.z` not found in input csv which has columns:",
      colnames(all_peaks)
    )
  )

  assertthat::assert_that(
    "z" %in% colnames(all_peaks),
    msg = cat(
      "Column `z` not found in input csv which has columns:",
      colnames(all_peaks)
    )
  )

  #filtering by peptide length and ALC
  all_peaks <- all_peaks %>%
    dplyr::filter(.data$length >= min_peptide_length &
                    .data$ALC.... >= min_ALC)

  #checking to make sure we have PSMs left
  if(nrow(all_peaks) == 0){

    stop("No PSMs passed min_peptide_length and min_ALC filters.")

  }


  #splitting by source file if necessary and running generate_PEAKS_tags
  if(length(unique(all_peaks$Source.File)) == 1){

    print("Only 1 source file detected in Source.File column.")

    #removing extension from Source.File
    source_dataset <- gsub("\\.\\w+$", "", unique(all_peaks$Source.File))
    print(paste("Source File without extension:", source_dataset))


    tag_df <- generate_PEAKS_tags(all_peaks_denovo = all_peaks,
                                  output_dir = output_dir,
                                  LC_threshold = LC_threshold,
                                  prolineBlock = prolineBlock,
                                  min_length = min_tag_length,
                                  ppm_error = ppm_error,
                                  return_tags = return_tags,
                                  new_mods = new_mods,
                                  hydrogen_mass = hydrogen_mass)

    print("Tag generation complete, returning a single dataframe.")
    return(tag_df)

  } else {

    print("Multiple source files detected in Source.File column.")
    all_peaks$Source.File <- factor(all_peaks$Source.File)
    peaks_result_list <- split(all_peaks, all_peaks$Source.File)

    #removing extension from list names which are the source files
    list_names <- names(peaks_result_list)
    new_names <- gsub("\\.\\w+$", "", list_names)
    #adding names to peaks_result_list so they will be carried through to output
    names(peaks_result_list) <- new_names

    print("Starting tag generation ...")

    #generating tags for each file
    peaks_tags_list <- purrr::map(
      peaks_result_list,
      generate_PEAKS_tags,
      output_dir = output_dir,
      LC_threshold = LC_threshold,
      prolineBlock = prolineBlock,
      min_length = min_tag_length,
      ppm_error = ppm_error,
      return_tags = return_tags,
      new_mods = new_mods,
      hydrogen_mass = hydrogen_mass
    )
    print("Tag generation complete, returning a list of dataframes.")
    return(peaks_tags_list)
  }



}








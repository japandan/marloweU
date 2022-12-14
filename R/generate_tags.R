#'Calculate precursor mass error in PPM
#'
#'This function calculates the mass error in PPM between the assigned peptide
#'mass and the precursor mass of the spectrum. **Note:** The mass of hydrogen
#'can made a significant difference in the final precursor error in ppm. Ex:
#'Using 1 vs 1.00**78246** can change the ppm error by up to 26 ppm. (The
#'standard filtering is +/- 15 ppm).
#'
#'@param peptide_neutral_mass A single value or numeric vector containing peptide neutral masses
#'@param spectrum_mz A single value or numeric vector containing the spectrum m/z value
#'@param charge A single value or numeric vector containing the spectrum charge
#'@param hydrogen_mass mass of hydrogen, default 1.00727646627 comes from
#'  [NIST
#'   proton molar mass](https://physics.nist.gov/cgi-bin/cuu/Category?view=html&Atomic+and+nuclear.x=123&Atomic+and+nuclear.y=21)
#'
#'@return A single value or numeric vector containing mass error in PPM. It uses the formula:
#'  \preformatted{ mz_peptide <- (peptide_neutral_mass + (hydrogen_mass *
#'  charge)) / charge mz_diff <- spectrum_mz - mz_peptide ppm_error <- (mz_diff
#'  / mz_peptide) * 1e+06}
#'@export
#'
#'@author Sarah C. Jenson
#'
#' @examples
#' \donttest{
#' \dontrun{
#'
#' Biod_novor$calc_ppm <-
#'      calc_ppm(peptide_neutral_mass = Biod_novor$pepMass.denovo.,
#'      spectrum_mz = Biod_novor$mz.data.,
#'      charge = Biod_novor$z)
#'}}
calc_ppm <- function(peptide_neutral_mass,
                          spectrum_mz,
                          charge,
                          hydrogen_mass = 1.00727646627){
  assertthat::assert_that(length(peptide_neutral_mass) == length(spectrum_mz),
                          msg = "peptide_neutral_mass and spectrum_mz need to be the same length")
  assertthat::assert_that(length(peptide_neutral_mass) == length(charge),
                          msg = "peptide_neutral_mass and charge need to be the same length")
  assertthat::assert_that(is.numeric(peptide_neutral_mass),
                          msg = "peptide_neutral_mass needs to be a numeric vector")
  assertthat::assert_that(is.numeric(spectrum_mz),
                          msg = "spectrum_mz needs to be a numeric vector")
  assertthat::assert_that((is.numeric(hydrogen_mass)),
                          msg = "hydrogen_mass needs to be a numeric number")
  assertthat::assert_that(is.numeric(charge),
                          msg = "charge needs to be a numeric vector")

    mz_peptide <-
        (peptide_neutral_mass + (hydrogen_mass * charge)) / charge

    mz_diff <- spectrum_mz - mz_peptide

    ppm_error <- (mz_diff / mz_peptide) * 1e+06

    return(ppm_error)
}


#' Convert from ppm error to mass in Da
#'
#'
#'
#' @param ppm A single numeric value representing mass error in ppm.
#' @param base_mass A single or vector of numeric values representing the mass
#'   your ppm error is based on.
#'
#' @return Returns a mass or vector of masses in Da that is different from the
#'   input `base_mass` by the specified `ppm`.
#'
#' @author Sarah C. Jenson
#'
#' @export
#' @importFrom assertthat assert_that
ppm_to_da <- function(ppm, base_mass)
{

  #checking inputs
  assertthat::assert_that(length(ppm) == 1, msg = "ppm must be a single value" )
  assertthat::assert_that(is.numeric(ppm), msg = "ppm must be numeric")
  assertthat::assert_that(is.numeric(base_mass), msg = "base_mass must be numeric")

  if (ppm == 0)
  {
    return(base_mass)
  } else {
    error_mass <- ((ppm * base_mass) / (1e+06)) + base_mass

    return(error_mass)
  }
}



#' Calculate the neutral mass from m/z and charge
#'
#' This function calculates the neutral mass of a peptide using the spectrum m/z
#' and charge. It uses the formula `neutral_mass <- (mz * z) - (z *
#' hydrogen_mass)`. **Note:** The mass of hydrogen can made a significant
#' difference when the neutral mass is used to calculate the precursor error in
#' ppm. Ex: Using 1 vs 1.00**78246** can change the ppm error by up to 26
#' ppm. (The standard filtering is +/- 15 ppm).
#'
#' @param mz vector of numeric m/z values
#' @param charge vector of charges corrisponding to the m/z values
#' @param hydrogen_mass mass of hydrogen, default is 1.00727646627 and comes
#'   from
#'   [NIST
#'   proton molar mass](https://physics.nist.gov/cgi-bin/cuu/Category?view=html&Atomic+and+nuclear.x=123&Atomic+and+nuclear.y=21)
#'
#' @return A numeric vector containing the neutral masses.
#'
#' @export
#'
#' @author Sarah C. Jenson
#'
#' @examples
#' \donttest{
#' \dontrun{
#' resultdf$neutral_mass <-
#'               calc_neutral_mass{mz = resultdf$mz.data,
#'                                charge = resultdf$z}
#'
#' calc_neutral_mass(mz = 301.1857, z = 3) #900.5336
#'
#' }}
calc_neutral_mass <- function(mz, charge, hydrogen_mass = 1.00727646627){

    #checking inputs
    assertthat::assert_that(length(mz) == length(charge),
                            msg = "mz and z must be same length")
    assertthat::assert_that(is.numeric(mz), msg = "mz must be numeric")
    assertthat::assert_that(is.numeric(charge), msg = "z must be numeric")

    #calculating neutral mass
    neutral_mass <- (mz * charge) - (charge * hydrogen_mass)
    return(neutral_mass)
}



#' Generate high confidence tags from Novor results
#'
#' This function is called by generate_tags and is not meant to be used by the user.
#'
#' @param denovo_df_row a dataframe with a single row. The peptide column should
#'   contain clean sequences with no modifications.
#' @param LC_threshold minimum Local Confidence score to include an amino acid
#'   in a tag.
#' @param pep_col Name of the column containing the peptide sequence. Should
#'   have no modificaiton info.
#' @param scan_col name of the column which contains the scan numbers.
#' @param min_length Minimum tag length to return.
#' @param min_mass_ppm ppm error to use when calculating the minimum mass,
#'   default -15 and must be <=0.
#' @param max_mass_ppm ppm error to use when calculating the maximum  mass,
#'   default +15 and must be >= 0.
#' @param trypsin_digest Should the tags be digested with trypsin?
#' @param prolineBlock should P after R or K prevent splitting at that site.

#'
#' @return A dataframe with columns 'tag' (containing the tag string),
#'   'pep_mass' (containing the mass of the original souce peptide),
#'   'original_peptide' (containing the source peptide),
#'   'scan_num' (containing the scan number of the source PSM), and
#'   'LC_threshold' (containing the LC_threshold).
#' @author Sarah C. Jenson
#' @importFrom stringr str_c
#' @importFrom MakeSearchSim digestPeptides
generate_tag_single_row <- function(denovo_df_row,
                                    LC_threshold,
                                    pep_col,
                                    scan_col,
                                    min_length,
                                    min_mass_ppm,
                                    max_mass_ppm,
                                    trypsin_digest,
                                    prolineBlock) {
    #checking inputs
    assertthat::assert_that(nrow(denovo_df_row) == 1,
                msg = paste("generate_tag_single_row got a dataframe with",
                                        nrow(denovo_df_row), "rows."))


    #other input checking done in generate_tags()

    # extracting local confidence scores from the
    # aaScore column
    LC_scores <- as.integer(unlist(strsplit(denovo_df_row$aaScore,
                                   split = "-"),
                                   use.names = FALSE))

    # making LC_scores has the right length
    assertthat::assert_that(
        length(LC_scores) == nchar(denovo_df_row[, pep_col]),
        msg = "LC_scores and peptide have different length.")

    # making a dataframe with one row per amino acid.
    # One column holds the amino acid and the second
    # column holds its associated local confidence
    # score
    pepdf <- data.frame(LC_score = LC_scores,
                        amino_acids = unlist(strsplit(denovo_df_row[,
                         pep_col], split = ""), use.names = FALSE),
                        stringsAsFactors = FALSE)

    # replacing all amino acids with LC score <
    # LC_threshold with X
    pepdf$amino_acids[pepdf$LC_score < LC_threshold] <- "X"

    # collapsing the amino acids into a single string,
    # then splitting by X, which returns a vector of
    # strings that contains the tags. X's are returned
    # as empty strings Ex: XXXABCXXDEF would return '',
    # 'ABC', '', 'DEF'
    tags <- unlist(strsplit(stringr::str_c(pepdf$amino_acids,
                                     collapse = ""),
                            split = "X"), use.names = FALSE)
    # removing empty strings
    tags <- tags[tags != ""]

    # if tags is empty return NULL
    if (length(tags) == 0) {
        return(NULL)
    }

    # digesting tags and filtering for length
    if (trypsin_digest == TRUE) {
        tags <- MakeSearchSim::digestPeptides(tags,prolineBlock = prolineBlock,
                                              discardLessThan = min_length)
    } else {
        tags <- tags[nchar(tags) >= min_length]
    }


    #calculating min and max mass in Da corrsiponsing to ppm_min and ppm_max
    min_mass <- ppm_to_da(ppm = min_mass_ppm, base_mass = denovo_df_row$adj_spectrum_mass)
    max_mass <- ppm_to_da(ppm = max_mass_ppm, base_mass = denovo_df_row$adj_spectrum_mass)



    #making sure adj_spectrum_mass is in between min_mass and max_mass
    assertthat::assert_that((min_mass < denovo_df_row$adj_spectrum_mass) &
                              (denovo_df_row$adj_spectrum_mass < max_mass),
                            msg = paste0(
                              "spectrum mass is not between min_mass and max_mass calculated using min_mass_ppm: ",
                              min_mass_ppm,
                              " and max_mass_ppm: ",
                              max_mass_ppm
                            )
    )


    # if tags is empty return NULL
    if (length(tags) == 0) {
        return(NULL)
    } else {
        # creating output dataframe with a row for each
        # distinct tag. Has columns for the tag string and
        # the original peptide mass
        output_df <- data.frame(tag = tags,
                                pep_mass = denovo_df_row$adj_spectrum_mass,
                                min_mass = min_mass,
                                max_mass = max_mass,
                                min_mass_ppm = min_mass_ppm,
                                max_mass_ppm = max_mass_ppm,
                                original_peptide = denovo_df_row[, pep_col],
                                scan_num = denovo_df_row[, scan_col],
                                LC_threshold = LC_threshold,
                                stringsAsFactors = FALSE)

        return(output_df)
    }

}


#' Generate tags from Novor result dataframe
#'
#' Takes a dataframe of Novor results and returns a dataframe containing all the
#' high confidence tags that results from applying the `LC_threshold` to
#' the Local Confidence scores in each peptide. Tags are digested with trypsin
#' if `trypsin_digest = TRUE`. Default column names corrispond to the
#' output of [MakeSearchSim::importSingleNovor()] from the
#' MakeSearchSim package.
#'
#' @param denovo_results a dataframe of Novor results with one row per PSM. Must
#'   contain a column named 'aaScore' which contains the Local Confidence
#'   scores. Names of columns containing other necessary information can be
#'   specified by the user below.
#' @param LC_threshold minimum Local Confidence score to include an amino acid
#'   in a tag. Default 80.
#' @param pep_col name of the column containing the modification free peptide
#'   sequence. Default 'cleanseq'.
#' @param mz_col name of the column containing the m/z values for each spectrum.
#'   Default 'mzPrecursor'
#' @param charge_col name of the column containing the charge for each spectrum.
#'   Default 'ChargePrecursor'.
#' @param scan_col name of the column which contains the scan numbers. Default
#'   'ScanNum'.
#' @param mod_col name of the column containing the peptide sequence with
#'   modification information. Default 'peptide'.
#' @param min_length Minimum tag length to return. Default 5.
#' @param min_mass_ppm ppm error to use when calculating the minimum mass,
#'   default -15 and must be <= 0.
#' @param max_mass_ppm ppm error to use when calculating the maximum  mass,
#'   default +15 and must be >= 0.
#' @param trypsin_digest Should the tags be digested with trypsin? Default =
#'   TRUE
#' @param prolineBlock should P after R or K prevent splitting at that site,
#'   Default TRUE
#' @param hydrogen_mass mass of hydrogen to use when calculating neutral mass.
#'   Default is 1.00727646627 which comes from
#'   [NIST
#'    proton molar mass](https://physics.nist.gov/cgi-bin/cuu/Category?view=html&Atomic+and+nuclear.x=123&Atomic+and+nuclear.y=21)
#' @param new_mods A named vector containing the masses of any additional
#'   modifications other then methionine oxidation or cysteine
#'   carbamidomethylation/alkylation. Vector names must be the short names that
#'   appear within the `()` in the peptides in the `mod_col`. Ex: ABC(Acetyl)HK
#'   => c(Acetyl = 42.010565). It is recommended you get the masses from
#'   unimod.org. Default is NA.
#'
#' @return A dataframe with columns:
#'   \describe{
#'     \item{tag}{containing the tag string}
#'     \item{pep_mass}{containing the neutral mass of
#'     the scan associated with the PSM as calculated
#'     by [calc_neutral_mass()] which is then adjusted to account for
#'     any modifications. Ex: A peptide with an oxidized methionine would have
#'     15.994915 subtracted from the precursor mass so it is comparable to an in
#'     silico calculated mass for the peptide.}
#'     \item{max_mass}{mass which has a `max_mass_ppm` ppm difference from
#'     the `pep_mass`}
#'     \item{min_mass}{mass which has a `min_mass_ppm` ppm
#'     difference from the `pep_mass`}
#'     \item{min_mass_ppm}{`min_mass_ppm` value supplied by user}
#'     \item{max_mass_ppm}{`max_mass_ppm` value supplied by user}
#'     \item{original_peptide}{containing the source peptide}
#'     \item{scan_num}{containing the scan number of the source PSM}
#'     \item{LC_threshold}{LC_threshold used}
#'     \item{sample_tag_id}{an incremental integer 1 -> number of rows}
#'     \item{Weight}{currently all 1, placeholder for later features}
#'   }
#'
#'   Output contains only tags with the specified
#'   `min_length`. These tags have been digested with trypsin if
#'   `trypsin_digest = TRUE`
#'
#' @author Sarah C. Jenson
#' @export
#' @importFrom plyr adply
#' @importFrom assertthat assert_that
#' @importFrom stringr str_detect str_extract_all
#' @examples \donttest{\dontrun{
#' #call showing custom mass_col and pep_col
#' denovo_tags <- generate_tags(denovo_results,
#'                   LC_threshold = 80,
#'                   pep_col = 'cleanseq',
#'                   mz_col = 'mzPrecursor',
#'                   charge_col = 'ChargePrecursor',
#'                   scan_col = 'ScanNum',
#'                   mod_col = 'peptide',
#'                   min_length = 5,
#'                   trypsin_digest = TRUE)
#' }}
generate_tags <- function(denovo_results,
                          LC_threshold = 80,
                          pep_col = "cleanseq",
                          mz_col = "mzPrecursor",
                          charge_col = "ChargePrecursor",
                          scan_col = "ScanNum",
                          mod_col = "peptide",
                          min_length = 5,
                          min_mass_ppm = -15,
                          max_mass_ppm = 15,
                          trypsin_digest = TRUE,
                          prolineBlock = TRUE,
                          hydrogen_mass = 1.00727646627,
                          new_mods = NULL) {
    # checking inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    assertthat::assert_that(mz_col %in% colnames(denovo_results),
        msg = paste("mz_col", mz_col, "not found in column names."))

    assertthat::assert_that(pep_col %in% colnames(denovo_results),
        msg = paste("pep_col", pep_col, "not found in column names."))

    assertthat::assert_that(scan_col %in% colnames(denovo_results),
        msg = paste("scan_col", scan_col, "not found in column names."))

    assertthat::assert_that(
        mod_col %in% colnames(denovo_results),
        msg = paste("mod_col", mod_col, "not found in column names.")
    )

    assertthat::assert_that("aaScore" %in% colnames(denovo_results),
        msg = paste("aaScore not found in column names."))

    assertthat::assert_that(is.character(denovo_results[, pep_col]),
                            msg = paste("pep_col", pep_col,
                                        "must be a character vector."))

    assertthat::assert_that(is.logical(trypsin_digest),
        msg = paste("trypsin_digest must be a logical. Currently it is a",
                    class(trypsin_digest)))

    assertthat::assert_that(is.logical(prolineBlock),
        msg = paste("prolineBlock must be a logical. Currently it is a",
                    class(prolineBlock)))

    #checking format of new_mods vector if it isn't NA
    if (!is.null(new_mods)) {
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


    #verifying min_mass_ppm is numeric
    assertthat::assert_that(
      is.numeric(min_mass_ppm),
      msg = paste0(
        "Min_mass_ppm must be numeric, currently it is a ",
        class(min_mass_ppm)
      )
    )

    #verifying max_mass_ppm is numeric
    assertthat::assert_that(
      is.numeric(max_mass_ppm),
      msg = paste0(
        "max_mass_ppm must be numeric, currently it is a ",
        class(max_mass_ppm)
      )
    )

    #verifying min_mass_ppm < max_mass_ppm
    assertthat::assert_that(
      min_mass_ppm < max_mass_ppm,
      msg = paste0(
        "min_mass_ppm must be less than max_mass_ppm. Currently min_mass_ppm is ",
        min_mass_ppm,
        " and max_mass_ppm is ",
        max_mass_ppm
      )
    )

    #verifying min_mass_ppm is <= 0
    assertthat::assert_that(min_mass_ppm <= 0,
                            msg = paste0("min_mass_ppm must be <= 0. Currently it is: ",
                                         min_mass_ppm))


    #verifying max_mass_ppm is >=0
    assertthat::assert_that(max_mass_ppm >= 0,
                            msg = paste0("max_mass_ppm must be >= 0. Currently it is: ",
                                         max_mass_ppm))

    #adjusting spectrum mass~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #calculating neutral mass of spectrum
    denovo_results$spectrum_neutral_mass <-
        calc_neutral_mass(mz = denovo_results[, mz_col],
                          charge = denovo_results[, charge_col],
                          hydrogen_mass = hydrogen_mass)

    #detecting PSMs with modifications
    denovo_results$has_mods <-
        stringr::str_detect(denovo_results[, mod_col], "\\(")

    denovo_results$mods <-
        stringr::str_extract_all(denovo_results[,mod_col],
                        pattern = "(?<=\\()([A-Z]|[a-z])+(?=\\))")

    denovo_results$mods[denovo_results$has_mods == FALSE] <- NA

    # named vector of modification names (as they appear in novor results) and
    # their unimod masses#
    mod_vector <- c(O = 15.994915, Cam = 57.021464)


    #adding new_mods to default mods
    if (!is.null(new_mods)) {
        mod_vector <- c(mod_vector, new_mods)
    }

    #helper function which looks up modification masses and returns sum
    sum_mod_mass <- function(mods, mod_masses){
        if (anyNA(mods)) {
            return(0)
        } else{
            return(sum(mod_masses[mods], na.rm = TRUE))
        }
    }

    denovo_results$mod_mass <- unlist(lapply(denovo_results$mods,
                                             sum_mod_mass,
                                             mod_masses = mod_vector))

    #subtracting modification mass from spectrum neutral mass to get a mass without modifications
    denovo_results$adj_spectrum_mass <-
        denovo_results$spectrum_neutral_mass - denovo_results$mod_mass

    print(paste0(sum(denovo_results$has_mods), " PSMs had modifications."))


    #generating tags~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # using adply to call generate_tage_single_row on
    # each row of the input dataframe
    output <-
        plyr::adply(
            denovo_results,
            1,
            generate_tag_single_row,
            LC_threshold = LC_threshold,
            pep_col = pep_col,
            min_length = min_length,
            min_mass_ppm = min_mass_ppm,
            max_mass_ppm = max_mass_ppm,
            trypsin_digest = trypsin_digest,
            prolineBlock = prolineBlock,
            scan_col = scan_col,
            .expand = FALSE,
            .id = NULL
        )

    assertthat::assert_that(nrow(output) > 0,
                msg = "No tags were generated, is LC_threshold too high?")

    # removing duplicates
    output <- unique(output)

    # adding id number
    output$sample_tag_id <- rownames(output)

    # adding weights (placeholder for later)
    output$Weight <- 1


    return(output)

}


#' Import Novor and generate tags
#'
#' This function is meant to be used within a funcitonal like mdply or map2 to
#' import Novor and generate tags in one step to avoid having to hold the Novor
#' results for a large batch of files in memmory. It calls
#' [MakeSearchSim::importSingleNovor()] then
#' [MakeSearchSim::filterNovor()] then [generate_tags()] on
#' a single dataset.
#'
#'
#'
#' @param Novor_Output_File path to a raw novor output file
#' @param MGF_file path to the MGF file corresponding to the 'NovorOutputFile',
#'   default NA.
#' @param novor_suffix The portion of the input file name that should be
#'   replaced with '_tags.RData' to make the output file name. Defaults to
#'   '_Novor.csv'. Ex: File1_Novor.csv => File1_tags.RData
#' @param clean_novor If TRUE(default), the metadata header will be removed and columns
#'   renamed from Novor_Output_File by [MakeSearchSim::importSingleNovor()].
#'   If FALSE, the Novor results file is read in directly as a csv,
#'   without calling [MakeSearchSim::importSingleNovor()].
#' @param output_tag_dir Directory where RData files containing tags should be
#'   written. Filename is File name is `MGF_file` with '.mgf' replaced with
#'   '_tags.RData'. Defaults to working directory.
#' @param return_tags Should the tags dataframe be returned? If FALSE (default)
#'   then the path to the tags dataframe is returned.
#' @param novor_ppm_min minimum allowed precursor error when filtering Novor results, default -15
#' @param novor_ppm_max maximum allowed precursor error when filtering Novor results, default +15
#' @param minNovorScore minimum allowed Novor score, default 20 for tag
#'   generation.
#' @param min_pep_length minimum allowed peptide length, default 6
#' @param LC_threshold minimum Local Confidence score to include an amino acid
#'   in a tag. Default 80.
#' @param min_tag_length Minimum tag length to return. Default 5.
#' @param trypsin_digest Should the tags be digested with trypsin? Default =
#'   TRUE
#' @param prolineBlock should P after R or K prevent splitting at that site,
#'   Default TRUE
#' @inheritParams generate_tags
#'
#' @return Returns either a dataframe with the tags or the path to the output file.
#'
#' @author Sarah C. Jenson
#'
#' @export
#' @importFrom assertthat assert_that
#' @importFrom MakeSearchSim importSingleNovor filterNovor
#' @importFrom utils read.csv
batch_tags <- function(Novor_Output_File,
                       MGF_file = NULL,
                       novor_suffix = "_Novor.csv",
                       clean_novor = TRUE,
                       output_tag_dir = ".",
                       return_tags = FALSE,
                       novor_ppm_min = -15,
                       novor_ppm_max = 15,
                       minNovorScore = 20,
                       min_pep_length = 6,
                       min_mass_ppm = -15,
                       max_mass_ppm = 15,
                       LC_threshold = 80,
                       min_tag_length = 5,
                       trypsin_digest = TRUE,
                       prolineBlock = TRUE,
                       hydrogen_mass = 1.00727646627,
                       new_mods = NULL)
{

  #input assertions
  assertthat::assert_that(file.exists(Novor_Output_File),
                          msg = "Novor_Output_File not found.")

  assertthat::assert_that(dir.exists(output_tag_dir),
                          msg = "output_tag_dir not found.")

  #importing Novor results
  if (clean_novor){
    raw_novor <-
      MakeSearchSim::importSingleNovor(NovorOutputFile = Novor_Output_File,
                                       MGF_file = MGF_file)
  } else {
    raw_novor <-
      read.csv(Novor_Output_File, stringsAsFactors = FALSE)
  }

 #filtering Novor results
 filt_novor <-
   MakeSearchSim::filterNovor(
     raw_novor,
     minNovorScore = minNovorScore,
     ppmErrorMin = novor_ppm_min,
     ppmErrorMax = novor_ppm_max,
     minLength = min_pep_length
   )

 #generating tags
 tags_mass_factor_df <-
   generate_tags(
     filt_novor,
     LC_threshold = LC_threshold,
     min_length = min_tag_length,
     min_mass_ppm = min_mass_ppm,
     max_mass_ppm = max_mass_ppm,
     trypsin_digest = trypsin_digest,
     prolineBlock = prolineBlock,
     hydrogen_mass = hydrogen_mass,
     new_mods = new_mods
   )

 #writing tags dataframe
 tags_filenname <-
   gsub(novor_suffix, "_tags.RData", basename(Novor_Output_File), fixed = TRUE)
 tags_path <- paste(normalizePath(output_tag_dir), tags_filenname,
                    sep = "/")

 save(tags_mass_factor_df, file = tags_path)

 if (return_tags == TRUE)
 {
   return(tags_mass_factor_df)
 } else{
   return(tags_path)
 }

}






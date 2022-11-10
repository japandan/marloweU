

#' Make heatmap by organism
#'
#' @param score_file_path path to csv with organism scores
#' @param overlap_matrix dataframe containing the overlap matrix to use for plotting
#' @param output_dir Path to the output directory, defaults to the current
#'   working directory.
#' @param output_type Extension associated with format to save the plot in. Must
#'   be one of "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", "bmp",
#'   "svg" or "wmf".
#' @param min_pos_score minimim pos_score to filter groups by, Default 1.
#' @param min_org_hits minimum number of organism hits required to display organism. Default 1.
#' @param score_suffix suffix to remove from the score_file_path to leave the
#'   dataset name, defaults to "_organism_scores.csv"
#' @param lineage_info_df An optional dataframe containing the NCBI lineage of
#'   all KEGG organisms. The `node_tax_name` column is used for shorter organism
#'   labels. Must have at minimum columns `kegg_id` and `node_tax_name`. Default
#'   `NULL`
#' @param ... additional arguments passed to \code{\link[ggplot2]{ggsave}}
#'
#' @return returns a ggplot2 object
#'
#' @author Sarah C. Jenson
#'
#' @export
#'
#' @importFrom readr read_csv cols col_character col_double
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr select mutate filter distinct inner_join full_join
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_sub
#' @import ggplot2
organism_heatmap <- function(score_file_path,
                          overlap_matrix,
                          output_dir = ".",
                          output_type = "jpeg",
                          min_org_hits = 1,
                          min_pos_score = 1,
                          score_suffix = "_organism_scores.csv",
                          lineage_info_df = NULL,
                          ...) {

  #checking inputs
  assertthat::assert_that(file.exists(score_file_path),
                          msg = "score_file_path not found.")

  assertthat::assert_that(dir.exists(output_dir), msg = "output_dir not found.")


  assertthat::assert_that(any(output_type %in% c("eps", "ps", "tex",
                                                 "pdf", "jpeg", "tiff",
                                                 "png", "bmp", "svg", "wmf")),
                          msg = "invalid output_type, must be one of: eps, ps, tex,
                                                 pdf, jpeg, tiff,
                                                 png, bmp, svg, or wmf")
  assertthat::assert_that(is.data.frame(overlap_matrix) | is.matrix(overlap_matrix),
                          msg = "overlap_matrix must be a dataframe or a matrix.")

  assertthat::assert_that(is.numeric(min_org_hits) & (min_org_hits >= 0),
                          msg = "min_org_hits must be numeric and >= 0")

  assertthat::assert_that(is.numeric(min_pos_score) & (min_pos_score >= 0),
                          msg = "min_pos_score must be numeric and >= 0")

  assertthat::assert_that(grepl(score_suffix, score_file_path, fixed = TRUE),
                          msg = "score_suffix not found in score_file_path.")


  if (!is.null(lineage_info_df)) {

    assertthat::assert_that(is.data.frame(lineage_info_df),
                            msg = "lineage_info_df must be a dataframe.")

    assertthat::assert_that("kegg_id" %in% colnames(lineage_info_df),
                            msg = "kegg_id column not found in lineage_info_df")

    assertthat::assert_that("node_tax_name" %in% colnames(lineage_info_df),
                            msg = "node_tax_name column not found in lineage_info_df")

  }


  #getting dataset name
  dataset_name <- gsub(score_suffix, "", basename(score_file_path), fixed = TRUE)

  print(paste0("Starting ", dataset_name))

  #Loading inputs
  scores <- readr::read_csv(score_file_path,
    col_types = readr::cols(
    kegg_id = readr::col_character(),
    organism_id = readr::col_double(),
    kegg_org_code = readr::col_character(),
    organism_name = readr::col_character(),
    species_name = readr::col_character(),
    genus_name = readr::col_character(),
    peptide_count = readr::col_double(),
    protein_count = readr::col_double(),
    strong_peptide_count = readr::col_double(),
    taxonomic_group = readr::col_double(),
    organism_hits = readr::col_double(),
    group_hits = readr::col_double(),
    orgscores = readr::col_double(),
    pos_scores = readr::col_double()
  ))


  #filtering scores
  filt_scores <- scores %>%
    dplyr::filter((.data$organism_hits >= min_org_hits)
                  & (.data$pos_scores >= min_pos_score))

  #making sure there are enough organisms
  if(nrow(filt_scores) <= 1){

    print(paste0("Only ",
                 nrow(filt_scores),
                 " organisms remain after filtering. Returning NA"))
    return(NA)
  }

  #adding column for which group has the top score and
  #which organism has the most hits
  filt_scores <- filt_scores %>%
    dplyr::mutate(is_top_group = dplyr::if_else(.data$pos_scores ==
                                                  max(.data$pos_scores),
                                                TRUE, FALSE)) %>%
    dplyr::mutate(is_top_org = dplyr::if_else(.data$organism_hits ==
                                                max(.data$organism_hits),
                                              TRUE, FALSE)) %>%
    dplyr::arrange(.data$taxonomic_group)


  #extracting kegg_ids
  kegg_ids <- filt_scores %>%
    dplyr::pull(.data$kegg_id) %>%
    unique()


  #verifying kegg_ids exist in matrix
  assertthat::assert_that(all(kegg_ids %in% colnames(overlap_matrix)) &
                            all(kegg_ids %in% rownames(overlap_matrix)),
                          msg = "All kegg_ids not found in overlap_matrix")

  #subsetting matrix to contain only filtered kegg_ids
  filt_matrix <- overlap_matrix[kegg_ids,kegg_ids]

  #making long version of score matrix
  long_overlap <- filt_matrix %>%
    tibble::rownames_to_column(var = "kegg_y") %>%
    tidyr::pivot_longer(cols = -.data$kegg_y,
                        names_to = "kegg_x",
                        values_to = "overlap")


  #adding taxonomic groups and pos scores to y axis
  long_overlap <- filt_scores %>%
    dplyr::select(.data$taxonomic_group,
                  .data$pos_scores,
                  .data$is_top_group,
                  .data$is_top_org,
                  .data$kegg_id) %>%
    dplyr::distinct() %>%
    dplyr::full_join(long_overlap, by = c("kegg_id" = "kegg_y")) %>%
    dplyr::rename(kegg_y = .data$kegg_id)

  #adding organism hits
  long_overlap <- filt_scores %>%
    dplyr::select(.data$kegg_id, .data$organism_hits) %>%
    dplyr::mutate(kegg_x = .data$kegg_id,
                  kegg_y = .data$kegg_id) %>%
    dplyr::full_join(long_overlap, by = c("kegg_x", "kegg_y"))


  #making organisms ordered factors so order will be retained in plot
  long_overlap <- long_overlap %>%
    dplyr::mutate(kegg_y = factor(.data$kegg_y, levels = kegg_ids, ordered = TRUE),
                  kegg_x = factor(.data$kegg_x, levels = kegg_ids, ordered = TRUE))


  #making new taxonomic group label that includes pos score
  long_overlap <- long_overlap %>%
    dplyr::mutate(pos_score_string  = format(.data$pos_scores, digits=0)) %>%
    dplyr::mutate(group_label = paste0("TG",
                                       .data$taxonomic_group,
                                       ": ",
                                       .data$pos_score_string)) %>%
    dplyr::mutate(group_label= factor(.data$group_label,
                                      levels = unique(.data$group_label),
                                      ordered = TRUE))



  #making abbreviated organism labels
  org_labels <- filt_scores %>%
    dplyr::select(.data$kegg_id,
                  .data$organism_name,
                  .data$genus_name,
                  .data$species_name) %>%
    dplyr::distinct() %>%
    dplyr::mutate(short_name = paste0(stringr::str_sub(.data$genus_name, 1, 1),
                                      ". ",
                                      .data$species_name))

  #adding node_taxon_name from lineage_info_df if it supplied
  if(!is.null(lineage_info_df)){
    org_labels <- lineage_info_df %>%
      dplyr::select(.data$kegg_id, .data$node_tax_name) %>%
      dplyr::right_join(org_labels, by = "kegg_id")
  }


  #making graph
  graph <- ggplot2::ggplot(long_overlap, aes(x = .data$kegg_x,
                                             y = .data$kegg_y,
                                             fill = .data$overlap)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(aes(label = .data$organism_hits,
                           color = .data$is_top_group,
                           fontface = dplyr::if_else(.data$is_top_org, 2, 1)), na.rm = TRUE) +
    ggplot2::scale_color_manual(breaks = c(TRUE, FALSE), values = c("red", "black")) +
    ggplot2::guides(color = FALSE) +
    ggplot2::scale_x_discrete(breaks = org_labels$kegg_id,
                              labels = org_labels$short_name,
                              expand = c(0,0))


  #using node_tax_name for org labels if lineage_info_df was supplied
  if(!is.null(lineage_info_df)){

    graph <- graph + ggplot2::scale_y_discrete(breaks = org_labels$kegg_id,
                                               labels = org_labels$node_tax_name,
                                               expand = c(0,0))
  } else {

   graph <- graph + ggplot2::scale_y_discrete(breaks = org_labels$kegg_id,
                              labels = org_labels$organism_name,
                              expand = c(0,0))

  }

  #setting up colors and organism group splits and labels
    graph <- graph + ggplot2::theme_bw() +
    ggplot2::scale_fill_distiller(palette = "Blues", direction = -1) +
    ggplot2::facet_grid(.data$group_label ~ .,
                        scales = "free",
                        space = "free",
                        as.table = FALSE,
                        switch = "y") +
    ggplot2::ggtitle(dataset_name,
                     subtitle = paste0("Organism level with minimum hits: ",
                                       min_org_hits,
                                       " and min pos_score: ",
                                       min_pos_score)) +
    ggplot2::xlab("Organism") + ggplot2::ylab("Organism") +
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                   strip.background = element_rect(fill = NA),
                   strip.text.y.left = element_text(angle = 0),
                   panel.spacing=unit(0.1,"cm"))



  #making filename for output
  output_filename <- paste0(dataset_name,
           "_org_heatmap_min_o_",
           min_org_hits,
           "_min_p_",
           min_pos_score,
           ".",
           output_type)

  #saving graph to output_dir
  #adjusted paremeters for higher resolution
  ggplot2::ggsave(filename = output_filename,
                  path = normalizePath(output_dir),
                  plot = graph,
                  device = output_type,
                  dpi = 600,
                  ...)

  return(graph)

}

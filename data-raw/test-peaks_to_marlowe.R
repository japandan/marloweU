context("PEAKS Output")


#loading input csv
test_peaks_output <- read.csv(file = test_path("Rcom_GCH4_1_2_PEAKS.csv"),
                              header = TRUE,
                              stringsAsFactors = FALSE)

#loading expected function output
load(test_path("Rcom_GCH4_1_2_marlowe_input_truth_6_7.RData"))



test_that("convert_peaks_to_marlowe", {

  #running function
  marlowe_input_test <- generate_PEAKS_tags(all_peaks_denovo = test_peaks_output,
                                            output_dir = test_path(),
                                            LC_threshold = 80,
                                            prolineBlock = TRUE,
                                            min_length = 5,
                                            ppm_error = 15,
                                            return_tags = TRUE)

  expect_identical(marlowe_input_test, marlowe_input_truth)

  expect_true(file.exists(test_path("Rcom_GCH4_1_2_Peaks_tags.RData")))

})

#cleaning up created file
if(file.exists(test_path("Rcom_GCH4_1_2_Peaks_tags.RData"))){
  file.remove(test_path("Rcom_GCH4_1_2_Peaks_tags.RData"))
}


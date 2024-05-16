#' Convert Metazoogene database to DADA2 trainset
#' 
#' This function takes any metazoogene database and converts it into a DADA2
#' trainset. Important is to use the CSV file format of MZGdb. A file path or a
#' URL can e used as input file.

MZGdb_to_DADA2 <- function(input_file, output_file) {

  library(magrittr)
  
  # Download
  readr::read_csv(input_file, col_names = FALSE, trim_ws = FALSE) %>%
    # Select ID, Sequence and taxa columns
    dplyr::select(X34, X31, X33, X2) %>%
    # Separate taxa column into separate columns
    tidyr::separate(X34, into = as.character(1:21), sep = ";") %>% 
    # Merge wanted taxa columns and add ">"
    tidyr::unite(col = "Taxa", 1,4,5,7,8,9,10,11,12,16,18,X2, sep = ";" ) %>%
    dplyr::mutate(Taxa = paste(">", Taxa, sep = "")) %>%
    # Restructure into one single column vector
    dplyr::select(X33, Taxa, X31) %>% 
    tidyr::pivot_longer(cols = c("Taxa", "X31")) %>%
    dplyr::pull(value) %>%
    # Write vector to fasta
    base::write(output_file)
}

# Execute the function for COI and 18S respectively
# These URLs are for the data set including all organisms (microbes, plants, and animals) with recorded pressence in the North Pacific.
# Visit https://www.st.nmfs.noaa.gov/copepod/collaboration/metazoogene/atlas/html-src/data__MZGdbALL__o00.html for options.
# Be awhare that the function only accepts the csv format of the database.

MZGdb_to_DADA2(input_file = "https://www.st.nmfs.noaa.gov/copepod/collaboration/metazoogene/atlas/data-src/MZGdata-coi__MZGdbALL__o07__A.csv.gz",
               output_file = "Data/Metazoogene/MZGdb_COI_NPac_ALL_mode-A_v3.0.fasta")

MZGdb_to_DADA2(input_file = "https://www.st.nmfs.noaa.gov/copepod/collaboration/metazoogene/atlas/data-src/MZGdata-18s__MZGdbALL__o07__A.csv.gz",
               output_file = "Data/Metazoogene/MZGdb_18S_NPac_ALL_mode-A_v3.0.fasta")
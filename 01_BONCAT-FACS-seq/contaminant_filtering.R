# Contaminant filtering of BONCAT-FACS-seq taxonomic profiles
#
# Desmarais et al.
# "Piezophilic activity dominates the Puerto Rico Trench
# deep subseafloor biosphere"
#
# Genera identified in the reagent-only negative control were
# removed from all biological samples prior to downstream analyses.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
})

# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

contam_path <- "../metadata/contaminants_genus.xlsx"

# Directory containing genus-level Bracken output files
bracken_dir <- "kraken_bracken/bracken"

# ------------------------------------------------------------------
# Helper function
# ------------------------------------------------------------------

norm_genus <- function(x) {
  x <- as.character(x)
  x <- gsub("^Candidatus\\s+", "", x, ignore.case = TRUE)
  x <- sub("^[dpcofgs]__", "", x)
  x <- sub("_[A-Za-z]+$", "", x)
  trimws(x)
}

# ------------------------------------------------------------------
# Load contaminant genera
# ------------------------------------------------------------------

contaminants <- read_excel(contam_path)

# The contaminant spreadsheet contains one genus column.
contaminant_genera <- contaminants[[1]] %>%
  as.character() %>%
  norm_genus() %>%
  unique()

contaminant_genera <- contaminant_genera[
  !is.na(contaminant_genera) &
  contaminant_genera != ""
]

message(
  "Loaded ",
  length(contaminant_genera),
  " unique contaminant genera."
)

# ------------------------------------------------------------------
# Load Bracken genus tables
# ------------------------------------------------------------------

bracken_files <- list.files(
  bracken_dir,
  pattern = "bracken_genus.*\\.(txt|tsv|csv)$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(bracken_files) == 0) {
  stop("No genus-level Bracken files found in: ", bracken_dir)
}

read_bracken <- function(file) {

  dat <- read.delim(
    file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  sample_id <- basename(file) %>%
    sub("\\.bracken_genus.*$", "", ., ignore.case = TRUE)

  dat %>%
    transmute(
      sample = sample_id,
      genus = norm_genus(name),
      est_reads = as.numeric(new_est_reads)
    ) %>%
    filter(
      !is.na(genus),
      genus != "",
      !is.na(est_reads)
    )
}

genus_data <- lapply(bracken_files, read_bracken) %>%
  bind_rows()

# ------------------------------------------------------------------
# Remove contaminant genera
# ------------------------------------------------------------------

genus_filtered <- genus_data %>%
  filter(!genus %in% contaminant_genera)

# Recalculate relative abundance after contaminant removal
genus_filtered <- genus_filtered %>%
  group_by(sample) %>%
  mutate(
    total_filtered_reads = sum(est_reads, na.rm = TRUE),
    relative_abundance = est_reads / total_filtered_reads
  ) %>%
  ungroup()

# ------------------------------------------------------------------
# Save cleaned table
# ------------------------------------------------------------------

write.csv(
  genus_filtered,
  file = file.path(bracken_dir, "bracken_genus_FINAL_clean.csv"),
  row.names = FALSE
)

message("Contaminant filtering complete.")

# JS1/Atribacterota predicted-proteome physicochemical properties
#
# Desmarais et al.
# "Piezophilic activity dominates the Puerto Rico Trench
# deep subseafloor biosphere"
#
# Calculates physicochemical properties from Prodigal-predicted
# protein sequences used in comparative genomic analyses.

suppressPackageStartupMessages({
  library(Biostrings)
  library(Peptides)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(readr)
  library(tibble)
  library(fs)
})

# ------------------------------------------------------------
# Input/output directories
# ------------------------------------------------------------

faa_dir <- "proteins"
results_dir <- "results"

dir_create(results_dir)

# ------------------------------------------------------------
# Amino-acid groups
# ------------------------------------------------------------

aa_letters <- c(
  "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
  "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
)

acidic_aa <- c("D", "E")
basic_aa <- c("K", "R", "H")
charged_aa <- c("D", "E", "K", "R", "H")
hydrophobic_aa <- c("A", "V", "I", "L", "M", "F", "W", "Y")

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

clean_protein_sequence <- function(x) {
  x <- toupper(as.character(x))
  str_replace_all(x, "[^ARNDCQEGHILKMFPSTWYV]", "")
}

safe_pI <- function(seq) {
  tryCatch(
    Peptides::pI(seq),
    error = function(e) NA_real_,
    warning = function(w) NA_real_
  )
}

safe_gravy <- function(seq) {
  tryCatch(
    Peptides::hydrophobicity(seq, scale = "KyteDoolittle"),
    error = function(e) NA_real_,
    warning = function(w) NA_real_
  )
}

# ------------------------------------------------------------
# Calculate properties for one predicted proteome
# ------------------------------------------------------------

summarize_one_proteome <- function(faa_file) {

  genome_id <- path_ext_remove(path_file(faa_file))

  message("Processing: ", genome_id)

  aa <- readAAStringSet(faa_file)

  seqs <- map_chr(aa, clean_protein_sequence)
  seqs <- seqs[nchar(seqs) > 0]

  if (length(seqs) == 0) {
    warning("No valid protein sequences found for: ", genome_id)
    return(tibble(genome_id = genome_id))
  }

  protein_lengths <- nchar(seqs)

  all_seq <- paste0(seqs, collapse = "")
  total_aa <- nchar(all_seq)

  aa_split <- str_split(all_seq, "", simplify = FALSE)[[1]]
  aa_counts_raw <- table(aa_split)

  aa_count_named <- setNames(rep(0, length(aa_letters)), aa_letters)

  shared_aa <- intersect(names(aa_counts_raw), aa_letters)
  aa_count_named[shared_aa] <- as.numeric(aa_counts_raw[shared_aa])

  aa_fraction_named <- aa_count_named / sum(aa_count_named)

  protein_pI <- map_dbl(seqs, safe_pI)
  protein_gravy <- map_dbl(seqs, safe_gravy)

  acidic_fraction <- sum(
    aa_fraction_named[acidic_aa],
    na.rm = TRUE
  )

  basic_fraction <- sum(
    aa_fraction_named[basic_aa],
    na.rm = TRUE
  )

  charged_fraction <- sum(
    aa_fraction_named[charged_aa],
    na.rm = TRUE
  )

  hydrophobic_fraction <- sum(
    aa_fraction_named[hydrophobic_aa],
    na.rm = TRUE
  )

  tibble(
    genome_id = genome_id,

    protein_count = length(seqs),
    total_aa = total_aa,

    protein_length_mean = mean(protein_lengths, na.rm = TRUE),
    protein_length_median = median(protein_lengths, na.rm = TRUE),

    mean_pI = mean(protein_pI, na.rm = TRUE),
    median_pI = median(protein_pI, na.rm = TRUE),

    mean_gravy = mean(protein_gravy, na.rm = TRUE),
    median_gravy = median(protein_gravy, na.rm = TRUE),

    acidic_fraction = acidic_fraction,
    basic_fraction = basic_fraction,

    acidic_basic_ratio =
      acidic_fraction / basic_fraction,

    charged_fraction = charged_fraction,
    hydrophobic_fraction = hydrophobic_fraction,

    fraction_A = aa_fraction_named["A"],
    fraction_R = aa_fraction_named["R"],
    fraction_N = aa_fraction_named["N"],
    fraction_D = aa_fraction_named["D"],
    fraction_C = aa_fraction_named["C"],
    fraction_Q = aa_fraction_named["Q"],
    fraction_E = aa_fraction_named["E"],
    fraction_G = aa_fraction_named["G"],
    fraction_H = aa_fraction_named["H"],
    fraction_I = aa_fraction_named["I"],
    fraction_L = aa_fraction_named["L"],
    fraction_K = aa_fraction_named["K"],
    fraction_M = aa_fraction_named["M"],
    fraction_F = aa_fraction_named["F"],
    fraction_P = aa_fraction_named["P"],
    fraction_S = aa_fraction_named["S"],
    fraction_T = aa_fraction_named["T"],
    fraction_W = aa_fraction_named["W"],
    fraction_Y = aa_fraction_named["Y"],
    fraction_V = aa_fraction_named["V"]
  )
}

# ------------------------------------------------------------
# Process all predicted proteomes
# ------------------------------------------------------------

faa_files <- dir_ls(
  faa_dir,
  regexp = "\\.faa$",
  type = "file"
)

if (length(faa_files) == 0) {
  stop("No .faa files found in: ", faa_dir)
}

proteome_summary <- map_dfr(
  faa_files,
  summarize_one_proteome
)

# Remove failed/extremely small predictions, as in the
# original analysis.
proteome_summary <- proteome_summary %>%
  filter(
    !is.na(protein_count),
    protein_count >= 50
  )

# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

write_csv(
  proteome_summary,
  file.path(
    results_dir,
    "js1_proteome_physicochemistry_summary.csv"
  )
)

message(
  "Proteome physicochemistry analysis complete: ",
  nrow(proteome_summary),
  " genomes retained."
)

# JS1/Atribacterota KO composition: PCoA and PERMANOVA
#
# Desmarais et al.
# "Piezophilic activity dominates the Puerto Rico Trench
# deep subseafloor biosphere"
#
# Bray-Curtis PCoA of KO copy density normalized per 1,000
# predicted proteins, with PERMANOVA and dispersion tests.
#
# Input tables should contain:
#   genome_id
#   ko_id
#   copies_per_1000
# and genome metadata containing:
#   genome_id
#   clade3_plot
#   carbon

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(vegan)
})

# ------------------------------------------------------------
# Settings used in the final analysis
# ------------------------------------------------------------

group_levels <- c("Non-Hadal", "Clade A", "Clade B")

permanova_permutations <- 9999
set.seed(123)

# Singleton and core KOs were retained in the final analysis.
drop_singleton_KOs <- FALSE
drop_core_KOs <- FALSE

# ------------------------------------------------------------
# Input/output
# ------------------------------------------------------------

ko_file <- "results/KO_copy_density_by_genome.csv"
metadata_file <- "../metadata/JS1_genome_metadata.csv"

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

ko <- read_csv(ko_file, show_col_types = FALSE)
metadata <- read_csv(metadata_file, show_col_types = FALSE)

required_ko <- c("genome_id", "ko_id", "copies_per_1000")
required_meta <- c("genome_id", "clade3_plot", "carbon")

if (!all(required_ko %in% names(ko))) {
  stop(
    "KO table must contain: ",
    paste(required_ko, collapse = ", ")
  )
}

if (!all(required_meta %in% names(metadata))) {
  stop(
    "Metadata must contain: ",
    paste(required_meta, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Prepare metadata
# ------------------------------------------------------------

metadata <- metadata %>%
  mutate(
    genome_id = as.character(genome_id),

    clade3_plot = case_when(
      clade3_plot %in%
        c("Other", "Non-hadal", "Non-Hadal") ~ "Non-Hadal",

      clade3_plot %in%
        c("A", "Clade A", "trench_A") ~ "Clade A",

      clade3_plot %in%
        c("B", "Clade B", "trench_B") ~ "Clade B",

      TRUE ~ as.character(clade3_plot)
    ),

    clade3_plot = factor(
      clade3_plot,
      levels = group_levels
    ),

    carbon = trimws(as.character(carbon)),
    carbon = ifelse(
      is.na(carbon) | carbon == "",
      "Other / unknown",
      carbon
    ),
    carbon = factor(carbon)
  ) %>%
  filter(
    !is.na(clade3_plot)
  )

# ------------------------------------------------------------
# Prepare KO copy-density table
# ------------------------------------------------------------

ko_clean <- ko %>%
  mutate(
    genome_id = as.character(genome_id),
    ko_id = as.character(ko_id),
    copies_per_1000 = as.numeric(copies_per_1000)
  ) %>%
  filter(
    !is.na(genome_id),
    !is.na(ko_id),
    !is.na(copies_per_1000),
    copies_per_1000 > 0
  ) %>%
  semi_join(metadata, by = "genome_id") %>%
  group_by(genome_id, ko_id) %>%
  summarise(
    copies_per_1000 = sum(copies_per_1000),
    .groups = "drop"
  )

# ------------------------------------------------------------
# KO prevalence
# ------------------------------------------------------------

ko_prevalence <- ko_clean %>%
  distinct(genome_id, ko_id) %>%
  count(ko_id, name = "n_genomes_with_KO") %>%
  mutate(
    n_total_genomes = n_distinct(metadata$genome_id)
  )

if (drop_singleton_KOs) {
  ko_prevalence <- ko_prevalence %>%
    filter(n_genomes_with_KO > 1)
}

if (drop_core_KOs) {
  ko_prevalence <- ko_prevalence %>%
    filter(n_genomes_with_KO < n_total_genomes)
}

ko_clean <- ko_clean %>%
  semi_join(
    ko_prevalence,
    by = "ko_id"
  )

# ------------------------------------------------------------
# Construct genome × KO matrix
# ------------------------------------------------------------

ko_matrix_df <- ko_clean %>%
  select(
    genome_id,
    ko_id,
    copies_per_1000
  ) %>%
  pivot_wider(
    id_cols = genome_id,
    names_from = ko_id,
    values_from = copies_per_1000,
    values_fill = 0
  )

# Add zero rows for retained genomes without detected KOs
missing_genomes <- setdiff(
  metadata$genome_id,
  ko_matrix_df$genome_id
)

if (length(missing_genomes) > 0) {

  zero_rows <- tibble(
    genome_id = missing_genomes
  )

  for (x in setdiff(names(ko_matrix_df), "genome_id")) {
    zero_rows[[x]] <- 0
  }

  ko_matrix_df <- bind_rows(
    ko_matrix_df,
    zero_rows
  )
}

ko_matrix_df <- ko_matrix_df %>%
  arrange(genome_id)

ko_mat <- ko_matrix_df %>%
  select(-genome_id) %>%
  as.matrix()

rownames(ko_mat) <- ko_matrix_df$genome_id

# Remove invariant KO columns
ko_variance <- apply(
  ko_mat,
  2,
  var,
  na.rm = TRUE
)

ko_mat <- ko_mat[
  ,
  ko_variance > 0,
  drop = FALSE
]

if (nrow(ko_mat) < 3) {
  stop("At least three genomes are required.")
}

if (ncol(ko_mat) < 2) {
  stop("At least two variable KOs are required.")
}

# ------------------------------------------------------------
# Match metadata order to KO matrix
# ------------------------------------------------------------

metadata_ordered <- metadata %>%
  filter(genome_id %in% rownames(ko_mat)) %>%
  mutate(
    genome_id = factor(
      genome_id,
      levels = rownames(ko_mat)
    )
  ) %>%
  arrange(genome_id) %>%
  mutate(
    genome_id = as.character(genome_id)
  )

if (!identical(
  metadata_ordered$genome_id,
  rownames(ko_mat)
)) {
  stop("Metadata and KO matrix orders do not match.")
}

# ------------------------------------------------------------
# Bray-Curtis distance
# ------------------------------------------------------------

bray_dist <- vegdist(
  ko_mat,
  method = "bray"
)

# ------------------------------------------------------------
# PCoA
# ------------------------------------------------------------

pcoa <- cmdscale(
  bray_dist,
  eig = TRUE,
  k = 2,
  add = TRUE
)

scores <- as.data.frame(
  pcoa$points[, 1:2, drop = FALSE]
)

colnames(scores) <- c(
  "PCoA1",
  "PCoA2"
)

scores$genome_id <- rownames(scores)

positive_eigenvalues <- pcoa$eig[
  pcoa$eig > 0
]

variance_explained <- (
  positive_eigenvalues /
    sum(positive_eigenvalues)
)

scores <- scores %>%
  left_join(
    metadata_ordered,
    by = "genome_id"
  )

write_csv(
  scores,
  file.path(
    results_dir,
    "KO_BrayCurtis_PCoA_scores.csv"
  )
)

write_csv(
  tibble(
    axis = c("PCoA1", "PCoA2"),
    variance_explained = variance_explained[1:2]
  ),
  file.path(
    results_dir,
    "KO_BrayCurtis_PCoA_variance.csv"
  )
)

# ------------------------------------------------------------
# PERMANOVA: JS1 clade
# ------------------------------------------------------------

permanova_clade <- adonis2(
  bray_dist ~ clade3_plot,
  data = metadata_ordered,
  permutations = permanova_permutations
)

# ------------------------------------------------------------
# PERMANOVA: carbon/substrate category
# ------------------------------------------------------------

permanova_carbon <- adonis2(
  bray_dist ~ carbon,
  data = metadata_ordered,
  permutations = permanova_permutations
)

# ------------------------------------------------------------
# Marginal PERMANOVA:
# clade and carbon considered simultaneously
# ------------------------------------------------------------

permanova_clade_carbon <- adonis2(
  bray_dist ~ clade3_plot + carbon,
  data = metadata_ordered,
  permutations = permanova_permutations,
  by = "margin"
)

# ------------------------------------------------------------
# Pairwise clade PERMANOVA
# ------------------------------------------------------------

pairwise_permanova <- function(
  group_a,
  group_b
) {

  keep <- metadata_ordered %>%
    filter(
      clade3_plot %in%
        c(group_a, group_b)
    )

  sub_mat <- ko_mat[
    keep$genome_id,
    ,
    drop = FALSE
  ]

  sub_dist <- vegdist(
    sub_mat,
    method = "bray"
  )

  keep <- keep %>%
    mutate(
      pair_group = droplevels(
        factor(clade3_plot)
      )
    )

  result <- adonis2(
    sub_dist ~ pair_group,
    data = keep,
    permutations = permanova_permutations
  )

  tibble(
    comparison = paste(
      group_a,
      "vs",
      group_b
    ),
    R2 = result$R2[1],
    F = result$F[1],
    p = result$`Pr(>F)`[1]
  )
}

pairwise_results <- bind_rows(
  pairwise_permanova(
    "Non-Hadal",
    "Clade A"
  ),
  pairwise_permanova(
    "Non-Hadal",
    "Clade B"
  ),
  pairwise_permanova(
    "Clade A",
    "Clade B"
  )
) %>%
  mutate(
    p_BH = p.adjust(
      p,
      method = "BH"
    )
  )

# ------------------------------------------------------------
# Homogeneity of multivariate dispersion
# ------------------------------------------------------------

dispersion_clade <- betadisper(
  bray_dist,
  metadata_ordered$clade3_plot
)

dispersion_clade_test <- permutest(
  dispersion_clade,
  permutations = permanova_permutations
)

dispersion_carbon <- betadisper(
  bray_dist,
  metadata_ordered$carbon
)

dispersion_carbon_test <- permutest(
  dispersion_carbon,
  permutations = permanova_permutations
)

# ------------------------------------------------------------
# Save statistical results
# ------------------------------------------------------------

capture.output(
  permanova_clade,
  file = file.path(
    results_dir,
    "PERMANOVA_clade.txt"
  )
)

capture.output(
  permanova_carbon,
  file = file.path(
    results_dir,
    "PERMANOVA_carbon.txt"
  )
)

capture.output(
  permanova_clade_carbon,
  file = file.path(
    results_dir,
    "PERMANOVA_clade_plus_carbon.txt"
  )
)

write_csv(
  pairwise_results,
  file.path(
    results_dir,
    "PERMANOVA_pairwise_clades.csv"
  )
)

capture.output(
  dispersion_clade_test,
  file = file.path(
    results_dir,
    "betadisper_clade.txt"
  )
)

capture.output(
  dispersion_carbon_test,
  file = file.path(
    results_dir,
    "betadisper_carbon.txt"
  )
)

message(
  "KO Bray-Curtis PCoA/PERMANOVA analysis complete."
)

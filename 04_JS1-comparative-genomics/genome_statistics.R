# JS1/Atribacterota genome and proteome statistics
#
# Desmarais et al.
# "Piezophilic activity dominates the Puerto Rico Trench
# deep subseafloor biosphere"
#
# Pairwise Wilcoxon rank-sum tests among Non-Hadal, Clade A,
# and Clade B genomes. P values are adjusted within each metric
# across the three pairwise comparisons using Benjamini-Hochberg FDR.
#
# Biological/statistical unit: one genome.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
})

# ------------------------------------------------------------
# Input/output
# ------------------------------------------------------------

genome_file <- "results/genome_architecture_metrics.csv"
proteome_file <- "results/js1_proteome_physicochemistry_summary.csv"
metadata_file <- "../metadata/JS1_genome_metadata.csv"

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

genome <- read_csv(genome_file, show_col_types = FALSE)
proteome <- read_csv(proteome_file, show_col_types = FALSE)
metadata <- read_csv(metadata_file, show_col_types = FALSE)

# ------------------------------------------------------------
# Standardize JS1 clade names
# ------------------------------------------------------------

standardize_clade <- function(x) {

  x <- as.character(x)

  case_when(
    x %in% c(
      "other", "Other",
      "Non-hadal", "Non-Hadal", "Non_Hadal",
      "Nonhadal", "Non_hadal"
    ) ~ "Non-Hadal",

    x %in% c(
      "trench_A", "trench A", "Trench A",
      "Clade A", "clade A", "clade_A", "A",
      "hadal_A", "Hadal A"
    ) ~ "Clade A",

    x %in% c(
      "trench_B", "trench B", "Trench B",
      "Clade B", "clade B", "clade_B", "B",
      "hadal_B", "Hadal B"
    ) ~ "Clade B",

    TRUE ~ x
  )
}

metadata <- metadata %>%
  mutate(
    genome_id = as.character(genome_id),
    js1_clade = standardize_clade(js1_clade)
  ) %>%
  filter(
    js1_clade %in%
      c("Non-Hadal", "Clade A", "Clade B")
  )

# ------------------------------------------------------------
# Genome architecture / KO metrics
# ------------------------------------------------------------

genome_metrics <- genome %>%
  mutate(genome_id = as.character(genome_id)) %>%
  select(
    genome_id,
    corrected_genome_size_mbp,
    predicted_orfs_per_corrected_mbp,
    ko_assigned_orfs_percent
  )

# ------------------------------------------------------------
# Proteome physicochemical metrics
# ------------------------------------------------------------

proteome_metrics <- proteome %>%
  mutate(genome_id = as.character(genome_id)) %>%
  select(
    genome_id,
    median_pI,
    acidic_basic_ratio,
    charged_fraction,
    hydrophobic_fraction
  )

# ------------------------------------------------------------
# Combine manuscript metrics
# ------------------------------------------------------------

metric_df <- metadata %>%
  select(
    genome_id,
    js1_clade
  ) %>%
  left_join(
    genome_metrics,
    by = "genome_id"
  ) %>%
  left_join(
    proteome_metrics,
    by = "genome_id"
  )

metric_long <- metric_df %>%
  pivot_longer(
    cols = c(
      corrected_genome_size_mbp,
      predicted_orfs_per_corrected_mbp,
      ko_assigned_orfs_percent,
      median_pI,
      acidic_basic_ratio,
      charged_fraction,
      hydrophobic_fraction
    ),
    names_to = "metric_key",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric_key,

      corrected_genome_size_mbp =
        "Corrected genome size (Mbp)",

      predicted_orfs_per_corrected_mbp =
        "Predicted ORFs per corrected Mbp",

      ko_assigned_orfs_percent =
        "KO-assigned ORFs (%)",

      median_pI =
        "Median protein pI",

      acidic_basic_ratio =
        "Acidic:basic residue ratio",

      charged_fraction =
        "Charged residue fraction",

      hydrophobic_fraction =
        "Hydrophobic residue fraction"
    ),

    js1_clade = factor(
      js1_clade,
      levels = c(
        "Non-Hadal",
        "Clade A",
        "Clade B"
      )
    )
  ) %>%
  filter(
    !is.na(value),
    is.finite(value),
    !is.na(js1_clade)
  )

# ------------------------------------------------------------
# Pairwise comparisons used in manuscript
# ------------------------------------------------------------

pairwise_comparisons <- tibble::tribble(
  ~comparison, ~group1, ~group2,

  "Clade A vs Clade B",
  "Clade A",
  "Clade B",

  "Non-Hadal vs Clade A",
  "Non-Hadal",
  "Clade A",

  "Non-Hadal vs Clade B",
  "Non-Hadal",
  "Clade B"
)

# ------------------------------------------------------------
# Pairwise Wilcoxon function
# ------------------------------------------------------------

run_pairwise_wilcox <- function(df) {

  pmap_dfr(
    pairwise_comparisons,

    function(comparison, group1, group2) {

      x <- df %>%
        filter(js1_clade == group1) %>%
        pull(value)

      y <- df %>%
        filter(js1_clade == group2) %>%
        pull(value)

      x <- x[is.finite(x)]
      y <- y[is.finite(y)]

      if (length(x) < 2 ||
          length(y) < 2) {

        return(
          tibble(
            comparison = comparison,
            group1 = group1,
            group2 = group2,
            n_group1 = length(x),
            n_group2 = length(y),
            median_group1 = median(
              x,
              na.rm = TRUE
            ),
            median_group2 = median(
              y,
              na.rm = TRUE
            ),
            difference_group1_minus_group2 =
              NA_real_,
            statistic_W = NA_real_,
            p_value = NA_real_
          )
        )
      }

      wt <- suppressWarnings(
        wilcox.test(
          x,
          y,
          exact = FALSE
        )
      )

      median_x <- median(
        x,
        na.rm = TRUE
      )

      median_y <- median(
        y,
        na.rm = TRUE
      )

      tibble(
        comparison = comparison,
        group1 = group1,
        group2 = group2,
        n_group1 = length(x),
        n_group2 = length(y),
        median_group1 = median_x,
        median_group2 = median_y,

        difference_group1_minus_group2 =
          median_x - median_y,

        statistic_W =
          unname(wt$statistic),

        p_value =
          wt$p.value
      )
    }
  )
}

# ------------------------------------------------------------
# Run tests
# ------------------------------------------------------------

pairwise_results <- metric_long %>%
  group_by(metric) %>%
  group_modify(
    ~ run_pairwise_wilcox(.x)
  ) %>%
  ungroup() %>%
  group_by(metric) %>%
  mutate(
    BH_adjusted_p = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  mutate(
    significance = case_when(
      is.na(BH_adjusted_p) ~ NA_character_,
      BH_adjusted_p < 0.001 ~ "***",
      BH_adjusted_p < 0.01 ~ "**",
      BH_adjusted_p < 0.05 ~ "*",
      TRUE ~ "ns"
    ),

    test =
      "Wilcoxon rank-sum test",

    p_adjustment_method =
      paste(
        "Benjamini-Hochberg FDR",
        "within each metric"
      )
  )

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------

write_csv(
  pairwise_results,
  file.path(
    out_dir,
    "JS1_pairwise_Wilcoxon_metrics_BH_FDR.csv"
  )
)

write_csv(
  metric_long,
  file.path(
    out_dir,
    "JS1_metrics_for_statistics.csv"
  )
)

message(
  "JS1 pairwise genome/proteome statistics complete."
)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Input table should contain, per genome:
# genome_id
# genome_size_observed_mbp
# checkm_completeness
# checkm_contamination
# predicted_orfs
# ko_assigned_orfs

input_file <- "../metadata/JS1_genome_metrics_input.csv"
output_file <- "results/genome_architecture_metrics.csv"

dir.create("results", showWarnings = FALSE)

df <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    genome_size_contam_adjusted_corrected_mbp =
      genome_size_observed_mbp * 100 /
      (checkm_completeness - checkm_contamination),

    ORF_density_per_Mbp =
      predicted_orfs /
      genome_size_contam_adjusted_corrected_mbp,

    KO_assigned_ORF_percent =
      100 * ko_assigned_orfs / predicted_orfs
  )

write_csv(df, output_file)

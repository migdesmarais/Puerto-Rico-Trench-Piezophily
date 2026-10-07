# JS1 16S rRNA gene clade abundance
#
# Calculates total JS1 relative abundance and abundance of the
# three dominant 99%-identity JS1 sequence clusters.
#
# Inputs:
#   PRT_16S_zotus_tax_FINAL.txt = finalized ASV/zOTU count table
#   vsearch_99/JS1_clusters_99.uc = VSEARCH clustering output

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

# ------------------------------------------------------------
# Input
# ------------------------------------------------------------

asv_file <- "PRT_16S_zotus_tax_FINAL.txt"
uc_file  <- "vsearch_99/JS1_clusters_99.uc"

outdir <- "results"
dir.create(outdir, showWarnings = FALSE)

# ------------------------------------------------------------
# Read finalized ASV table
# ------------------------------------------------------------

asv <- read.delim(
  asv_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# First column = ASV/zOTU identifier
# Last column = taxonomy
id_col  <- names(asv)[1]
tax_col <- names(asv)[ncol(asv)]

sample_cols <- setdiff(
  names(asv),
  c(id_col, tax_col)
)

# Total reads in each sample BEFORE extracting JS1
sample_totals <- colSums(
  asv[, sample_cols],
  na.rm = TRUE
)

# ------------------------------------------------------------
# Extract JS1-affiliated sequence variants
# ------------------------------------------------------------

js1 <- asv %>%
  filter(
    str_detect(
      .data[[tax_col]],
      "(^|;)JS1(;|$)"
    )
  )

message(
  "JS1 sequence variants detected: ",
  nrow(js1)
)

# ------------------------------------------------------------
# Read VSEARCH 99%-identity clustering
# ------------------------------------------------------------

uc <- read.delim(
  uc_file,
  header = FALSE,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE,
  fill = TRUE
)

# VSEARCH .uc:
# V1 = record type
# V2 = cluster number
# V9 = query sequence label

cluster_membership <- uc %>%
  filter(V1 %in% c("S", "H")) %>%
  transmute(
    cluster = as.integer(V2),
    zotu = V9
  ) %>%
  distinct()

# ------------------------------------------------------------
# Associate JS1 zOTUs with 99% clusters
# ------------------------------------------------------------

js1_clusters <- js1 %>%
  rename(zotu = all_of(id_col)) %>%
  inner_join(
    cluster_membership,
    by = "zotu"
  )

# ------------------------------------------------------------
# Determine the three dominant clusters
#
# Dominance is defined by summed read abundance across the
# complete PRT 16S dataset.
# ------------------------------------------------------------

cluster_totals <- js1_clusters %>%
  group_by(cluster) %>%
  summarise(
    total_reads =
      sum(
        across(
          all_of(sample_cols)
        ),
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  arrange(desc(total_reads)) %>%
  mutate(
    dominance_rank = row_number()
  )

dominant_clusters <- cluster_totals %>%
  slice_head(n = 3) %>%
  mutate(
    js1_16S_clade =
      paste0(
        "JS1-16S clade ",
        dominance_rank
      )
  )

write_csv(
  dominant_clusters,
  file.path(
    outdir,
    "JS1_16S_dominant_clusters.csv"
  )
)

# ------------------------------------------------------------
# Calculate clade abundance per sample
# ------------------------------------------------------------

clade_counts <- js1_clusters %>%
  inner_join(
    dominant_clusters %>%
      select(
        cluster,
        js1_16S_clade
      ),
    by = "cluster"
  ) %>%
  select(
    zotu,
    js1_16S_clade,
    all_of(sample_cols)
  ) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "reads"
  ) %>%
  group_by(
    sample,
    js1_16S_clade
  ) %>%
  summarise(
    reads = sum(
      reads,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Total JS1 abundance per sample
# ------------------------------------------------------------

total_js1 <- js1 %>%
  select(
    all_of(id_col),
    all_of(sample_cols)
  ) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "reads"
  ) %>%
  group_by(sample) %>%
  summarise(
    JS1_reads = sum(
      reads,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    total_reads =
      sample_totals[sample],

    JS1_relative_abundance =
      100 * JS1_reads / total_reads
  )

# ------------------------------------------------------------
# Convert clade counts to relative abundance of total
# community reads
# ------------------------------------------------------------

clade_abundance <- clade_counts %>%
  mutate(
    total_reads =
      sample_totals[sample],

    relative_abundance =
      100 * reads / total_reads
  ) %>%
  select(
    sample,
    js1_16S_clade,
    relative_abundance
  ) %>%
  pivot_wider(
    names_from = js1_16S_clade,
    values_from = relative_abundance,
    values_fill = 0
  )

final <- total_js1 %>%
  select(
    sample,
    JS1_relative_abundance
  ) %>%
  left_join(
    clade_abundance,
    by = "sample"
  )

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write_csv(
  final,
  file.path(
    outdir,
    "JS1_16S_clade_relative_abundance.csv"
  )
)

message(
  "JS1 16S abundance analysis complete."
)

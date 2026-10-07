# JS1 16S rRNA gene clade abundance
#
# Calculates total JS1 relative abundance and the abundance of
# the three dominant 99%-identity JS1 sequence clusters.
#
# The three dominant clusters were manually designated
# JS1-16S clades 1, 2 and 3 for the manuscript.
#
# Inputs:
#   PRT_16S_zotus_tax_FINAL.txt
#   vsearch_99/JS1_clusters_99.uc

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
# Read finalized ASV/zOTU table
# ------------------------------------------------------------

asv <- read.delim(
  asv_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# First column = ASV/zOTU identifier
# Last column  = taxonomy

id_col  <- names(asv)[1]
tax_col <- names(asv)[ncol(asv)]

sample_cols <- setdiff(
  names(asv),
  c(id_col, tax_col)
)

# Total reads in each sample before extracting JS1

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

# VSEARCH .uc format:
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
# Associate JS1 zOTUs with 99%-identity clusters
# ------------------------------------------------------------

js1_clusters <- js1 %>%
  rename(zotu = all_of(id_col)) %>%
  inner_join(
    cluster_membership,
    by = "zotu"
  )

# ------------------------------------------------------------
# Identify the three dominant JS1 clusters
#
# Dominance is determined from summed read abundance across
# the complete PRT 16S dataset.
# ------------------------------------------------------------

cluster_totals <- js1_clusters %>%
  group_by(cluster) %>%
  summarise(
    total_reads = sum(
      as.matrix(
        pick(all_of(sample_cols))
      ),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(total_reads))

dominant_clusters <- cluster_totals %>%
  slice_head(n = 3)

write_csv(
  cluster_totals,
  file.path(
    outdir,
    "JS1_16S_cluster_total_abundance.csv"
  )
)

message(
  "Three dominant 99%-identity clusters: ",
  paste(dominant_clusters$cluster, collapse = ", ")
)

# ------------------------------------------------------------
# Assign manuscript clade labels
# ------------------------------------------------------------
#
# IMPORTANT:
#
# Clades 1, 2 and 3 are manuscript labels assigned to the
# three dominant 99%-identity groups. They are NOT assigned
# automatically according to abundance rank.
#
# After VSEARCH clustering, replace the values below with the
# cluster IDs corresponding to the manuscript clades.
#
# Example:
#
# clade_mapping <- tibble(
#   cluster = c(4, 0, 2),
#   js1_16S_clade = c(
#     "JS1-16S clade 1",
#     "JS1-16S clade 2",
#     "JS1-16S clade 3"
#   )
# )

clade_mapping <- tibble(
  cluster = c(
    NA_integer_,
    NA_integer_,
    NA_integer_
  ),
  js1_16S_clade = c(
    "JS1-16S clade 1",
    "JS1-16S clade 2",
    "JS1-16S clade 3"
  )
)

# Stop rather than silently assigning incorrect clade labels

if (any(is.na(clade_mapping$cluster))) {
  stop(
    paste0(
      "Enter the VSEARCH cluster IDs corresponding to ",
      "JS1-16S clades 1, 2 and 3 in clade_mapping."
    )
  )
}

# Check that the manually assigned clades are the three
# dominant clusters

if (!setequal(
  clade_mapping$cluster,
  dominant_clusters$cluster
)) {
  warning(
    paste0(
      "The manually assigned clades do not correspond exactly ",
      "to the three most abundant JS1 clusters."
    )
  )
}

write_csv(
  clade_mapping,
  file.path(
    outdir,
    "JS1_16S_clade_mapping.csv"
  )
)

# ------------------------------------------------------------
# Calculate clade abundance per sample
# ------------------------------------------------------------

clade_counts <- js1_clusters %>%
  inner_join(
    clade_mapping,
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
# Calculate total JS1 abundance per sample
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
# Calculate abundance of each JS1 clade relative to the
# complete 16S community
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

# ------------------------------------------------------------
# Final table
# ------------------------------------------------------------

final <- total_js1 %>%
  select(
    sample,
    JS1_relative_abundance
  ) %>%
  left_join(
    clade_abundance,
    by = "sample"
  )

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

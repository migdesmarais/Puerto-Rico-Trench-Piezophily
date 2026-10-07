suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(broom)
})

genome_file <- "results/genome_architecture_metrics.csv"
proteome_file <- "results/js1_proteome_physicochemistry_summary.csv"
metadata_file <- "../metadata/JS1_genome_metadata.csv"

outdir <- "results"
dir.create(outdir, showWarnings = FALSE)

genome_df <- read_csv(genome_file, show_col_types = FALSE)
proteome_df <- read_csv(proteome_file, show_col_types = FALSE)
metadata <- read_csv(metadata_file, show_col_types = FALSE)

group_levels <- c(
  "Non-Hadal",
  "Clade A",
  "Clade B"
)

standardize_clade <- function(x) {

  x <- as.character(x)

  case_when(
    x %in% c(
      "Other",
      "other",
      "Non-Hadal",
      "Non-hadal",
      "Non_hadal",
      "Nonhadal",
      "Non-Hadal / other",
      "Non-hadal / other"
    ) ~ "Non-Hadal",

    x %in% c(
      "Clade A",
      "clade A",
      "clade_A",
      "A",
      "hadal_A",
      "trench_A",
      "trench A",
      "Trench A"
    ) ~ "Clade A",

    x %in% c(
      "Clade B",
      "clade B",
      "clade_B",
      "B",
      "hadal_B",
      "trench_B",
      "trench B",
      "Trench B"
    ) ~ "Clade B",

    TRUE ~ x
  )
}

clean_type <- function(x) {

  x <- as.character(x)

  case_when(
    str_to_upper(x) == "MAG" ~ "MAG",
    str_to_upper(x) == "SAG" ~ "SAG",
    TRUE ~ x
  )
}

model_input <- metadata %>%
  mutate(
    genome_id = as.character(genome_id),
    clade3_plot = standardize_clade(clade3_plot),
    clade3_plot = factor(
      clade3_plot,
      levels = group_levels
    ),
    type = factor(
      clean_type(type),
      levels = c("MAG", "SAG")
    ),
    checkm_completeness =
      as.numeric(checkm_completeness),
    checkm_contamination =
      as.numeric(checkm_contamination)
  ) %>%
  left_join(
    genome_df %>%
      select(
        genome_id,
        genome_size_contam_adjusted_corrected_mbp,
        ORF_density_per_Mbp,
        KO_assigned_ORF_percent
      ),
    by = "genome_id"
  ) %>%
  left_join(
    proteome_df %>%
      select(
        genome_id,
        charged_fraction
      ),
    by = "genome_id"
  )

model_long <- model_input %>%
  pivot_longer(
    cols = c(
      genome_size_contam_adjusted_corrected_mbp,
      ORF_density_per_Mbp,
      KO_assigned_ORF_percent,
      charged_fraction
    ),
    names_to = "response",
    values_to = "response_value"
  )

candidate_covariates <- c(
  "checkm_completeness",
  "checkm_contamination",
  "type"
)

get_usable_covariates <- function(df) {

  usable <- character()

  for (covar in candidate_covariates) {

    x <- df[[covar]]

    if (is.numeric(x)) {

      vals <- x[!is.na(x)]

      if (
        length(vals) >= 2 &&
        sd(vals) > 0
      ) {
        usable <- c(
          usable,
          covar
        )
      }

    } else {

      vals <- droplevels(
        factor(
          x[!is.na(x)]
        )
      )

      if (nlevels(vals) >= 2) {
        usable <- c(
          usable,
          covar
        )
      }
    }
  }

  usable
}

run_adjusted_lm <- function(
    response_name,
    ref_group
) {

  df <- model_long %>%
    filter(
      response == response_name,
      !is.na(response_value),
      !is.na(clade3_plot)
    ) %>%
    mutate(
      clade3_plot =
        relevel(
          factor(
            clade3_plot,
            levels = group_levels
          ),
          ref = ref_group
        )
    )

  usable_covariates <-
    get_usable_covariates(df)

  formula_text <- paste(
    "response_value ~",
    paste(
      c(
        "clade3_plot",
        usable_covariates
      ),
      collapse = " + "
    )
  )

  model_df <- df %>%
    select(
      genome_id,
      response_value,
      clade3_plot,
      all_of(
        usable_covariates
      )
    ) %>%
    drop_na()

  if (nrow(model_df) < 8) {
    return(
      tibble(
        response = response_name,
        reference_group = ref_group,
        n = nrow(model_df),
        term = NA_character_,
        estimate = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_
      )
    )
  }

  fit <- lm(
    as.formula(
      formula_text
    ),
    data = model_df
  )

  broom::tidy(fit) %>%
    transmute(
      response = response_name,
      reference_group = ref_group,
      n = nrow(model_df),
      term,
      estimate,
      std_error = std.error,
      statistic,
      p_value = p.value
    )
}

responses <- c(
  "genome_size_contam_adjusted_corrected_mbp",
  "ORF_density_per_Mbp",
  "KO_assigned_ORF_percent",
  "charged_fraction"
)

results <- bind_rows(
  lapply(
    responses,
    function(resp) {
      bind_rows(
        run_adjusted_lm(
          resp,
          "Non-Hadal"
        ),
        run_adjusted_lm(
          resp,
          "Clade A"
        )
      )
    }
  )
)

write_csv(
  results,
  file.path(
    outdir,
    "JS1_adjusted_linear_models.csv"
  )
)

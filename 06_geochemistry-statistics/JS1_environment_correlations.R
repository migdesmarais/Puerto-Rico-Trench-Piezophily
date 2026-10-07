# JS1 relative abundance versus environmental variables
#
# Reproduces Supplementary Table 13.
#
# JS1 observations are matched to TOC, delta13Corg and CaCO3 using
# maximum-cardinality, minimum-total-depth-offset one-to-one matching.
# O2 is linearly interpolated to JS1 horizons.
#
# Spearman correlations are calculated for all JS1-containing sites
# together and separately for PR05, PR07 and PR08.
# BH correction is applied separately within each site-analysis family.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(readr)
})

# Reuse matching and interpolation functions
source("transect_correlations.R", local = TRUE)

input_file <- "../metadata/PRT_geochem_long_format.xlsx"

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)

param_js1   <- "JS1"
param_toc   <- "TOC"
param_d13c  <- "d13C"
param_caco3 <- "CaCO3"

tol_js1_toc   <- 60
tol_js1_d13c  <- 60
tol_js1_caco3 <- 60


source_df <- read_excel(
  input_file,
  sheet = 1
) %>%
  rename_with(
    ~ gsub(
      "\\s+",
      "_",
      trimws(.x)
    )
  ) %>%
  mutate(
    Station =
      as.character(Station),

    water_depth =
      suppressWarnings(
        as.numeric(water_depth)
      ),

    sediment_depth =
      suppressWarnings(
        as.numeric(sediment_depth)
      ),

    Parameter =
      as.character(Parameter),

    Value =
      suppressWarnings(
        as.numeric(Value)
      )
  ) %>%
  fill(
    Station,
    water_depth
  )


build_js1_var <- function(
    param_name
) {

  source_df %>%
    filter(
      Parameter == param_name,
      is.finite(sediment_depth),
      is.finite(Value)
    ) %>%
    group_by(
      Station,
      water_depth,
      sediment_depth
    ) %>%
    summarise(
      value =
        mean(
          Value,
          na.rm = TRUE
        ),
      .groups = "drop"
    ) %>%
    arrange(
      Station,
      sediment_depth
    )
}


js1_df <-
  build_js1_var(
    param_js1
  )

toc_df <-
  build_js1_var(
    param_toc
  )

d13c_df <-
  build_js1_var(
    param_d13c
  )

caco3_df <-
  build_js1_var(
    param_caco3
  )


param_o2 <-
  detect_o2_param(
    source_df
  )

o2_df <-
  build_js1_var(
    param_o2
  )


# -------------------------------------------------------------------------
# Match JS1 to environmental variables
# -------------------------------------------------------------------------

matched_js1_toc <-
  one_to_one_match_optimal(
    js1_df,
    toc_df,
    tol_js1_toc
  )

matched_js1_d13c <-
  one_to_one_match_optimal(
    js1_df,
    d13c_df,
    tol_js1_d13c
  )

matched_js1_caco3 <-
  one_to_one_match_optimal(
    js1_df,
    caco3_df,
    tol_js1_caco3
  )


# Our interpolation helper returns O2 as x and JS1 as y.
# Spearman rho is symmetric, so this does not affect rho or P.

matched_js1_o2 <-
  interpolate_o2_to_target(
    o2_df,
    js1_df
  )


matched_js1_depth <-
  js1_df %>%
  transmute(
    Station,
    anchor_depth =
      sediment_depth,

    matched_depth =
      sediment_depth,

    depth_diff_cmbsf = 0,

    x_value = value,

    y_value =
      sediment_depth
  )


# -------------------------------------------------------------------------
# Analysis registry
# -------------------------------------------------------------------------

registry <- tibble(
  Predictor = c(
    "O2",
    "TOC",
    "d13C",
    "CaCO3",
    "Depth"
  ),

  matched_data = list(
    matched_js1_o2,
    matched_js1_toc,
    matched_js1_d13c,
    matched_js1_caco3,
    matched_js1_depth
  )
)


safe_spearman_js1 <- function(d) {

  d <- d %>%
    filter(
      is.finite(x_value),
      is.finite(y_value)
    )

  if (
    nrow(d) < 3 ||
    n_distinct(d$x_value) < 2 ||
    n_distinct(d$y_value) < 2
  ) {

    return(
      tibble(
        rho = NA_real_,
        raw_P = NA_real_,
        n = nrow(d)
      )
    )
  }


  test <- suppressWarnings(
    cor.test(
      d$x_value,
      d$y_value,
      method = "spearman",
      exact = FALSE,
      alternative = "two.sided"
    )
  )


  tibble(
    rho =
      unname(
        test$estimate
      ),

    raw_P =
      test$p.value,

    n =
      nrow(d)
  )
}


# -------------------------------------------------------------------------
# Run pooled and site-specific correlations
# -------------------------------------------------------------------------

sites <- c(
  "All",
  "PR05",
  "PR07",
  "PR08"
)


results <- map_dfr(
  sites,
  function(site_name) {

    map_dfr(
      seq_len(
        nrow(registry)
      ),
      function(i) {

        d <-
          registry$matched_data[[i]]

        if (site_name != "All") {

          d <- d %>%
            filter(
              Station ==
                site_name
            )
        }


        safe_spearman_js1(
          d
        ) %>%
          mutate(
            Site =
              site_name,

            Predictor =
              registry$Predictor[i]
          )
      }
    )
  }
)


# Separate BH family for each pooled/site analysis,
# as implemented in the final Supplementary Table 13 code.

results <- results %>%
  group_by(
    Site
  ) %>%
  mutate(
    BH_adjusted_P =
      p.adjust(
        raw_P,
        method = "BH"
      )
  ) %>%
  ungroup() %>%
  select(
    Site,
    Predictor,
    n,
    rho,
    raw_P,
    BH_adjusted_P
  )


write_csv(
  results,
  file.path(
    out_dir,
    "Supplementary_Table_13_JS1_correlations.csv"
  )
)

print(
  results,
  n = Inf
)

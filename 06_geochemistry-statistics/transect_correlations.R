# Puerto Rico Trench geochemistry and microbial-proxy correlations
#
# Reproduces the 15 transect-wide Spearman correlations reported
# in Supplementary Table 4.
#
# Non-O2 variables are paired within station using maximum-cardinality,
# minimum-total-depth-offset one-to-one matching without replacement.
# O2 is linearly interpolated to target horizons and set to zero at and
# below the first measured zero-O2 depth.
#
# All 15 P values are adjusted together using Benjamini-Hochberg.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(readr)
})

file_path <- "../metadata/PRT_geochem_long_format.xlsx"
out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Parameter names
param_cells_rep1 <- "Cells rep1"
param_cells_rep2 <- "Cells rep2"
param_atp        <- "ATP"
param_asv        <- "Rarefied ASV"
param_toc        <- "TOC"
param_tn         <- "TN"
param_d13c       <- "d13C"

# Maximum permitted depth offsets (cmbsf)
tol_atp_cells  <- 50
tol_asv_cells  <- 15
tol_atp_asv    <- 10

tol_toc_cells  <- 150
tol_d13c_cells <- 60
tol_tn_cells   <- 20

tol_toc_atp    <- 50
tol_d13c_atp   <- 60
tol_tn_atp     <- 50

tol_toc_asv    <- 15
tol_d13c_asv   <- 60
tol_tn_asv     <- 15


# -------------------------------------------------------------------------
# Read source data
# -------------------------------------------------------------------------

source_df <- read_excel(file_path, sheet = 1) %>%
  rename_with(~ gsub("\\s+", "_", trimws(.x))) %>%
  mutate(
    Station        = as.character(Station),
    water_depth    = suppressWarnings(as.numeric(water_depth)),
    sediment_depth = suppressWarnings(as.numeric(sediment_depth)),
    Parameter      = as.character(Parameter),
    Value          = suppressWarnings(as.numeric(Value))
  ) %>%
  fill(Station, water_depth)


# -------------------------------------------------------------------------
# Detect oxygen parameter
# -------------------------------------------------------------------------

detect_o2_param <- function(long_df) {

  p <- unique(long_df$Parameter)
  p_lower <- tolower(trimws(p))

  exact <- c(
    "o2",
    "oxygen",
    "dissolved oxygen"
  )

  hit <- which(p_lower %in% exact)

  if (length(hit) == 1)
    return(p[hit])

  hit <- which(
    str_detect(
      p_lower,
      "(^|[^a-z0-9])o2([^a-z0-9]|$)|oxygen"
    )
  )

  if (length(hit) == 1)
    return(p[hit])

  stop("Could not uniquely detect oxygen parameter.")
}

param_o2 <- detect_o2_param(source_df)


# -------------------------------------------------------------------------
# Build variable tables
# -------------------------------------------------------------------------

build_var_df <- function(long_df, param_name, positive_only = FALSE) {

  out <- long_df %>%
    filter(Parameter == param_name) %>%
    group_by(
      Station,
      water_depth,
      sediment_depth
    ) %>%
    summarise(
      value = mean(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(
      is.finite(value),
      is.finite(sediment_depth),
      is.finite(water_depth)
    ) %>%
    arrange(
      Station,
      sediment_depth
    )

  if (positive_only)
    out <- out %>% filter(value > 0)

  out
}


cells_df <- source_df %>%
  filter(
    Parameter %in%
      c(param_cells_rep1, param_cells_rep2)
  ) %>%
  group_by(
    Station,
    water_depth,
    sediment_depth
  ) %>%
  summarise(
    value = mean(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    is.finite(value),
    value > 0,
    is.finite(sediment_depth),
    is.finite(water_depth)
  ) %>%
  arrange(
    Station,
    sediment_depth
  )

atp_df  <- build_var_df(source_df, param_atp,  TRUE)
asv_df  <- build_var_df(source_df, param_asv)
toc_df  <- build_var_df(source_df, param_toc)
tn_df   <- build_var_df(source_df, param_tn,   TRUE)
d13c_df <- build_var_df(source_df, param_d13c)
o2_df   <- build_var_df(source_df, param_o2)


# -------------------------------------------------------------------------
# Optimal one-to-one depth matching
#
# Maximizes number of pairs first.
# Among maximum-cardinality solutions, minimizes total absolute
# depth offset.
# -------------------------------------------------------------------------

one_to_one_match_optimal <- function(
    anchor_df,
    candidate_df,
    tolerance
) {

  station_ids <- intersect(
    unique(anchor_df$Station),
    unique(candidate_df$Station)
  )

  matched_list <- lapply(
    station_ids,
    function(stn) {

      a <- anchor_df %>%
        filter(Station == stn) %>%
        arrange(sediment_depth)

      b <- candidate_df %>%
        filter(Station == stn) %>%
        arrange(sediment_depth)

      n_a <- nrow(a)
      n_b <- nrow(b)

      if (n_a == 0 || n_b == 0)
        return(tibble())

      best_n <- matrix(
        0L,
        nrow = n_a + 1L,
        ncol = n_b + 1L
      )

      best_cost <- matrix(
        0,
        nrow = n_a + 1L,
        ncol = n_b + 1L
      )

      action <- matrix(
        "",
        nrow = n_a + 1L,
        ncol = n_b + 1L
      )

      if (n_a > 0)
        action[2:(n_a + 1L), 1] <- "skip_anchor"

      if (n_b > 0)
        action[1, 2:(n_b + 1L)] <- "skip_candidate"


      for (i in seq_len(n_a)) {

        for (j in seq_len(n_b)) {

          options <- tibble(
            action = c(
              "skip_anchor",
              "skip_candidate"
            ),
            n_pair = c(
              best_n[i, j + 1L],
              best_n[i + 1L, j]
            ),
            cost = c(
              best_cost[i, j + 1L],
              best_cost[i + 1L, j]
            ),
            priority = c(2L, 3L)
          )

          dz <- abs(
            a$sediment_depth[i] -
              b$sediment_depth[j]
          )

          if (dz <= tolerance) {

            options <- bind_rows(
              options,
              tibble(
                action = "match",
                n_pair =
                  best_n[i, j] + 1L,
                cost =
                  best_cost[i, j] + dz,
                priority = 1L
              )
            )
          }

          winner <- options %>%
            arrange(
              desc(n_pair),
              cost,
              priority
            ) %>%
            slice(1)

          best_n[i + 1L, j + 1L] <-
            winner$n_pair

          best_cost[i + 1L, j + 1L] <-
            winner$cost

          action[i + 1L, j + 1L] <-
            winner$action
        }
      }


      i <- n_a
      j <- n_b
      pairs <- list()

      while (i > 0L || j > 0L) {

        move <- action[
          i + 1L,
          j + 1L
        ]

        if (move == "match") {

          pairs[[length(pairs) + 1L]] <-
            tibble(
              Station = stn,

              anchor_depth =
                a$sediment_depth[i],

              matched_depth =
                b$sediment_depth[j],

              depth_diff_cmbsf =
                abs(
                  a$sediment_depth[i] -
                    b$sediment_depth[j]
                ),

              x_value = a$value[i],
              y_value = b$value[j]
            )

          i <- i - 1L
          j <- j - 1L

        } else if (move == "skip_anchor") {

          i <- i - 1L

        } else if (move == "skip_candidate") {

          j <- j - 1L

        } else {

          stop(
            "Matching backtrack failed for ",
            stn
          )
        }
      }

      bind_rows(rev(pairs))
    }
  )

  bind_rows(matched_list)
}


# -------------------------------------------------------------------------
# Oxygen interpolation
#
# Linear interpolation to target horizons.
# At/below first measured zero-O2 depth, O2 = 0.
# -------------------------------------------------------------------------

interpolate_o2_to_target <- function(
    o2_tbl,
    target_tbl
) {

  station_ids <- intersect(
    unique(o2_tbl$Station),
    unique(target_tbl$Station)
  )

  map_dfr(
    station_ids,
    function(stn) {

      o <- o2_tbl %>%
        filter(Station == stn) %>%
        group_by(sediment_depth) %>%
        summarise(
          value = mean(value, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        filter(
          is.finite(value),
          is.finite(sediment_depth)
        ) %>%
        arrange(sediment_depth) %>%
        mutate(
          value = pmax(value, 0)
        )

      target <- target_tbl %>%
        filter(Station == stn) %>%
        filter(
          is.finite(value),
          is.finite(sediment_depth)
        ) %>%
        arrange(sediment_depth)

      if (nrow(o) < 2 ||
          nrow(target) == 0)
        return(tibble())


      interpolated <- approx(
        x = o$sediment_depth,
        y = o$value,
        xout = target$sediment_depth,
        method = "linear",
        rule = 2
      )$y


      if (any(o$value <= 0)) {

        first_zero_depth <-
          min(
            o$sediment_depth[
              o$value <= 0
            ]
          )

        interpolated[
          target$sediment_depth >=
            first_zero_depth
        ] <- 0
      }


      tibble(
        Station = stn,

        anchor_depth =
          target$sediment_depth,

        matched_depth =
          target$sediment_depth,

        depth_diff_cmbsf = 0,

        x_value =
          pmax(interpolated, 0),

        y_value =
          target$value
      )
    }
  )
}


# -------------------------------------------------------------------------
# Construct matched datasets
# -------------------------------------------------------------------------

matched_atp_cells <-
  one_to_one_match_optimal(
    atp_df,
    cells_df,
    tol_atp_cells
  )

matched_asv_cells <-
  one_to_one_match_optimal(
    asv_df,
    cells_df,
    tol_asv_cells
  )

matched_atp_asv <-
  one_to_one_match_optimal(
    atp_df,
    asv_df,
    tol_atp_asv
  )


matched_toc_cells <-
  one_to_one_match_optimal(
    toc_df,
    cells_df,
    tol_toc_cells
  )

matched_d13c_cells <-
  one_to_one_match_optimal(
    d13c_df,
    cells_df,
    tol_d13c_cells
  )

matched_tn_cells <-
  one_to_one_match_optimal(
    tn_df,
    cells_df,
    tol_tn_cells
  )


matched_toc_atp <-
  one_to_one_match_optimal(
    toc_df,
    atp_df,
    tol_toc_atp
  )

matched_d13c_atp <-
  one_to_one_match_optimal(
    d13c_df,
    atp_df,
    tol_d13c_atp
  )

matched_tn_atp <-
  one_to_one_match_optimal(
    tn_df,
    atp_df,
    tol_tn_atp
  )


matched_toc_asv <-
  one_to_one_match_optimal(
    toc_df,
    asv_df,
    tol_toc_asv
  )

matched_d13c_asv <-
  one_to_one_match_optimal(
    d13c_df,
    asv_df,
    tol_d13c_asv
  )

matched_tn_asv <-
  one_to_one_match_optimal(
    tn_df,
    asv_df,
    tol_tn_asv
  )


matched_o2_toc <-
  interpolate_o2_to_target(
    o2_df,
    toc_df
  )

matched_o2_d13c <-
  interpolate_o2_to_target(
    o2_df,
    d13c_df
  )

matched_o2_tn <-
  interpolate_o2_to_target(
    o2_df,
    tn_df
  )


# -------------------------------------------------------------------------
# Supplementary Table 4 registry
# -------------------------------------------------------------------------

registry <- tibble(
  Comparison = c(
    "Cells vs total adenylate",
    "Cells vs rarefied ASV richness",
    "Rarefied ASV richness vs total adenylate",
    "Cells vs TOC",
    "Cells vs δ13Corg",
    "Cells vs TN",
    "Total adenylate vs TOC",
    "Total adenylate vs δ13Corg",
    "Total adenylate vs TN",
    "Rarefied ASV richness vs TOC",
    "Rarefied ASV richness vs δ13Corg",
    "Rarefied ASV richness vs TN",
    "O2 vs TOC",
    "O2 vs δ13Corg",
    "O2 vs TN"
  ),

  matched_data = list(
    matched_atp_cells,
    matched_asv_cells,
    matched_atp_asv,
    matched_toc_cells,
    matched_d13c_cells,
    matched_tn_cells,
    matched_toc_atp,
    matched_d13c_atp,
    matched_tn_atp,
    matched_toc_asv,
    matched_d13c_asv,
    matched_tn_asv,
    matched_o2_toc,
    matched_o2_d13c,
    matched_o2_tn
  )
)


safe_spearman <- function(d) {

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
    rho = unname(test$estimate),
    raw_P = test$p.value,
    n = nrow(d)
  )
}


results <- registry %>%
  mutate(
    statistics =
      map(
        matched_data,
        safe_spearman
      )
  ) %>%
  select(
    -matched_data
  ) %>%
  unnest(statistics) %>%
  mutate(
    BH_adjusted_P =
      p.adjust(
        raw_P,
        method = "BH"
      )
  )


write_csv(
  results,
  file.path(
    out_dir,
    "Supplementary_Table_4_correlations.csv"
  )
)

print(results)

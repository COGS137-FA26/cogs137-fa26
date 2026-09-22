# cs01_functions.R ==========================================================
# Shared helpers for CS01: Biomarkers of Recent Use
# Used in the CS01 lecture notes and the Lab 02 answer key.
# Requires: tidyverse, janitor, patchwork (for plot_cutoffs)
#
# The `## ---- name` lines mark sections that the lecture notes display with
# knitr::read_chunk(); keep them when editing.

## ---- cs01_windows
# ==== Time windows ============================================================
# Windows (minutes from start of smoking) differ by matrix, following the paper.
# `cut()` intervals are right-closed: (0, 30] is "0-30 min"; anything <= 0 is
# pre-smoking; the last window is open-ended.
cs01_windows <- list(
  WB = list(breaks = c(-Inf, 0, 30, 70, 100, 180, 210, 240, 270, 300, Inf),
            labels = c("pre-smoking", "0-30 min", "31-70 min", "71-100 min",
                       "101-180 min", "181-210 min", "211-240 min",
                       "241-270 min", "271-300 min", "301+ min")),
  OF = list(breaks = c(-Inf, 0, 30, 90, 180, 210, 240, 270, Inf),
            labels = c("pre-smoking", "0-30 min", "31-90 min", "91-180 min",
                       "181-210 min", "211-240 min", "241-270 min", "271+ min")),
  BR = list(breaks = c(-Inf, 0, 40, 90, 180, 210, 240, 270, Inf),
            labels = c("pre-smoking", "0-40 min", "41-90 min", "91-180 min",
                       "181-210 min", "211-240 min", "241-270 min", "271+ min"))
)

# The first window after smoking, for each matrix
cs01_early_windows <- c("0-30 min", "0-40 min")

# Compounds measured (not every matrix has every compound)
cs01_compounds <- c("cbn", "cbd", "thc", "thcoh", "thccooh", "thccooh_gluc",
                    "cbg", "thcv", "thca_a")

## ---- cs01_clean
# ==== Cleaning ================================================================
# Clean one raw CS01 file (Blood.csv, OF.csv, or Breath.csv).
#   fluid: "WB", "OF", or "BR" (determines the time windows)
clean_cs01 <- function(df, fluid = c("WB", "OF", "BR")) {
  fluid <- match.arg(fluid)
  win <- cs01_windows[[fluid]]

  df |>
    janitor::clean_names() |>
    rename(any_of(c(fluid_type   = "fluid",
                    thcoh        = "x11_oh_thc",
                    thccooh      = "thc_cooh",
                    thccooh_gluc = "thc_cooh_gluc",
                    thcv         = "thc_v",
                    thc          = "thc_pg_pad"))) |>
    mutate(
      treatment = fct_recode(treatment,
                             "5.9% THC (low dose)"   = "5.90%",
                             "13.4% THC (high dose)" = "13.40%"),
      treatment = fct_relevel(treatment, "Placebo", "5.9% THC (low dose)"),
      group = case_match(group,
                         "Experienced user"     ~ "Frequent user",
                         "Not experienced user" ~ "Occasional user",
                         .default = group),
      timepoint = cut(time_from_start, breaks = win$breaks, labels = win$labels)
    )
}

## ---- cs01_drop_dups
# Keep one measurement per person per matrix per time window (the first one)
drop_dups <- function(df) {
  df |>
    filter(!is.na(timepoint)) |>
    distinct(fluid_type, timepoint, id, .keep_all = TRUE)
}

## ---- cs01_to_long
# Wide -> long: one row per measurement of one compound
cs01_to_long <- function(df) {
  df |>
    pivot_longer(any_of(cs01_compounds),
                 names_to = "compound", values_to = "value")
}

## ---- cs01_sens_spec
# ==== Sensitivity & specificity ===============================================
# For every matrix x compound x time window x cutoff:
#   "used"     = received THC AND sample taken after smoking
#   "detected" = value at or above the cutoff
# Pre-smoking samples count as true negatives/false positives only.
sens_spec <- function(df_long, cutoffs = c(0.5, 1, 2, 5, 10)) {
  df_long |>
    filter(!is.na(value), !is.na(timepoint)) |>
    cross_join(tibble(cutoff = cutoffs)) |>
    mutate(used = treatment != "Placebo" & timepoint != "pre-smoking",
           detected = value >= cutoff) |>
    group_by(fluid_type, compound, timepoint, cutoff) |>
    summarize(TP = sum(used & detected),
              FN = sum(used & !detected),
              FP = sum(!used & detected),
              TN = sum(!used & !detected),
              .groups = "drop") |>
    mutate(N = TP + FN + FP + TN,
           sensitivity = TP / (TP + FN),
           specificity = TN / (TN + FP),
           ppv         = TP / (TP + FP),
           npv         = TN / (TN + FN),
           efficiency  = (TP + TN) / N)
}

## ---- cs01_theme
# ==== Plotting ================================================================
cs01_colors <- c("Frequent user" = "#19831C", "Occasional user" = "#A27FC9")

theme_cs01 <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          panel.grid = element_blank(),
          plot.title.position = "plot",
          strip.background = element_blank())
}

## ---- cs01_plot_time
# Every measurement over time for one matrix, one panel per compound
plot_scatter_time <- function(df_long, matrix) {
  df_long |>
    filter(fluid_type == matrix, !is.na(time_from_start), !is.na(value)) |>
    ggplot(aes(x = time_from_start, y = value, color = group)) +
    geom_point() +
    facet_wrap(~compound, scales = "free") +
    scale_color_manual(values = cs01_colors) +
    labs(title = matrix, x = "Time From Start (min)", y = "Measurement") +
    theme_cs01()
}

## ---- cs01_plot_boxplots
# Boxplots at the first post-smoking window for one matrix.
#   x: the variable to compare (e.g., group or treatment), unquoted
plot_boxplots <- function(df_long, matrix, x) {
  df_long |>
    filter(fluid_type == matrix,
           timepoint %in% cs01_early_windows,
           !is.na(value)) |>
    ggplot(aes(x = {{ x }}, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0),
               size = 0.8, color = "gray40") +
    facet_wrap(~compound, scales = "free") +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = cs01_colors) +
    labs(title = matrix, x = NULL, y = "Measurement (ng/mL)") +
    theme_cs01()
}

## ---- cs01_plot_cutoffs
# Sensitivity (left) & specificity (right) across time windows, one line per cutoff
plot_cutoffs <- function(ss, fluid, cpd, title = fluid) {
  df <- ss |>
    filter(fluid_type == fluid, compound == cpd) |>
    mutate(timepoint = factor(timepoint, levels = cs01_windows[[fluid]]$labels),
           cutoff = factor(cutoff))

  base <- function(y, ylab) {
    ggplot(df |> filter(!is.na({{ y }})),
           aes(x = timepoint, y = {{ y }} * 100, color = cutoff, group = cutoff)) +
      geom_line(linewidth = 1.2) +
      geom_point() +
      scale_y_continuous(limits = c(0, 100)) +
      scale_x_discrete(drop = FALSE) +
      # colorblind-friendly sequential palette (darker = higher cutoff)
      scale_color_viridis_d(name = "Cutoff (ng/mL)", end = 0.85, direction = -1) +
      labs(x = "Time Window", y = ylab) +
      theme_cs01() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  (base(sensitivity, "Sensitivity (%)") + labs(title = paste0(title, ": ", toupper(cpd)))) +
    base(specificity, "Specificity (%)") +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.title = element_text())
}

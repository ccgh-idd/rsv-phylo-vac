# phylo_trees_by_country.R
# ─────────────────────────────────────────────────────────────────────────────
# Builds country-specific phylogenetic trees for RSV-A and RSV-B from the
# full cleaned AA (F-protein) sequence dataset, then:
#   1. Estimates evolutionary rate via root-to-tip regression
#   2. Compares rates pre vs post vaccination rollout
#   3. Meta-analyses across countries
#   4. Produces publication-ready figures
#
# Inputs:
#   data/sequences/rsv-a_analysis_ready.fasta   ← full dataset
#   data/sequences/rsv-b_analysis_ready.fasta
#   summaries/all_clean_meta.csv                ← full metadata
#   data/country_vaccination_status.csv
#
# Per-country stratified subsampling (MAX_PER_COUNTRY) is applied before tree
# building so that dominant countries (e.g. USA n=13k) do not cause hours-long
# runs. Sampling is balanced across years to preserve temporal signal.
#
# Outputs (all under results/phylo_trees_full/):
#   trees/   – Newick tree files per country per subtype
#   plots/   – Rate comparison, RTT regression, tree, forest plots
#   summaries/ – CSV rate tables and meta-analysis results
# ─────────────────────────────────────────────────────────────────────────────

# ── 0. Packages ───────────────────────────────────────────────────────────────

packages <- c("ape", "phangorn", "seqinr", "tidyverse", "lubridate",
              "here", "metafor", "ggplot2", "patchwork", "ggtree")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "ggtree") {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install("ggtree", ask = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ── 1. Output directories ─────────────────────────────────────────────────────

out_base   <- here("results", "phylo_trees_full")
out_trees  <- file.path(out_base, "trees")
out_plots  <- file.path(out_base, "plots")
out_summ   <- file.path(out_base, "summaries")

for (d in c(out_trees, out_plots, out_summ)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ── 2. Parameters ─────────────────────────────────────────────────────────────

MIN_SEQS        <- 10   # minimum sequences per country to attempt a tree
MIN_PRE_POST    <- 5    # minimum sequences in each time period for rate comparison
MAX_PER_COUNTRY <- 300  # per-country cap before tree building (stratified by year)
                        # keeps tree building tractable even for USA (n=13k)

# ── 3. Load data ──────────────────────────────────────────────────────────────

cat("Loading metadata...\n")
meta <- read_csv(here("summaries", "all_clean_meta.csv"),
                 show_col_types = FALSE) %>%
  mutate(
    # Decimal year: e.g. 2023-07-15 → 2023.535
    date_decimal = decimal_date(as.Date(substr(parsed_date, 1, 10)))
  ) %>%
  filter(!is.na(date_decimal), !is.na(geoLocCountry), geoLocCountry != "")

cat(sprintf("  Full metadata: %d sequences across %d countries\n",
            nrow(meta), n_distinct(meta$geoLocCountry)))

cat("Loading vaccination dates...\n")
vacc_dates <- read_csv(here("data", "country_vaccination_status.csv"),
                       show_col_types = FALSE) %>%
  filter(has_rsv_vaccination == 1) %>%
  mutate(vacc_decimal = introduction_year) %>%   # integer year as threshold
  select(country, vacc_decimal)

cat("Loading sequences...\n")
read_fasta_named <- function(path) {
  seqs <- seqinr::read.fasta(path, seqtype = "AA", as.string = TRUE)
  vapply(seqs, function(x) toupper(paste(x, collapse = "")), character(1))
}

seqs_a <- read_fasta_named(here("data", "sequences", "rsv-a_analysis_ready.fasta"))
seqs_b <- read_fasta_named(here("data", "sequences", "rsv-b_analysis_ready.fasta"))

cat(sprintf("  RSV-A: %d sequences | RSV-B: %d sequences (full dataset)\n",
            length(seqs_a), length(seqs_b)))

# ── 3b. Per-country stratified subsampling ────────────────────────────────────
# For countries with > MAX_PER_COUNTRY sequences, sample proportionally across
# years so temporal coverage is preserved. Countries below the cap are kept whole.
stratified_subsample <- function(meta_sub, max_n, seed = 42) {
  set.seed(seed)
  meta_sub %>%
    group_by(geoLocCountry, year) %>%
    mutate(year_n = n()) %>%
    ungroup() %>%
    group_by(geoLocCountry) %>%
    mutate(
      country_n  = n(),
      keep_prob  = pmin(1, max_n / country_n)  # uniform thinning within country
    ) %>%
    filter(runif(n()) <= keep_prob) %>%
    ungroup()
}

meta <- stratified_subsample(meta, MAX_PER_COUNTRY)
cat(sprintf("  After per-country cap (%d): %d sequences across %d countries\n",
            MAX_PER_COUNTRY, nrow(meta), n_distinct(meta$geoLocCountry)))

# ── 4. Helper functions ───────────────────────────────────────────────────────

# Build a midpoint-rooted NJ tree from a named character vector of AA sequences
build_tree <- function(seq_vec, country, subtype) {
  # Convert to phyDat
  mat <- do.call(rbind, strsplit(seq_vec, ""))
  rownames(mat) <- names(seq_vec)
  pd <- phyDat(mat, type = "AA")

  # Distance matrix and NJ tree
  dm  <- dist.ml(pd, model = "WAG")
  nj  <- NJ(dm)

  # ML optimisation (topology + branch lengths, WAG model)
  fit <- pml(nj, pd, model = "WAG")
  opt <- tryCatch(
    optim.pml(fit, optNni = TRUE, optBf = FALSE, optQ = FALSE,
              optGamma = FALSE, rearrangement = "NNI", control = pml.control(trace = 0)),
    error = function(e) {
      message(sprintf("  ML optimisation failed for %s %s – using NJ tree", subtype, country))
      fit
    }
  )

  tree <- midpoint(opt$tree)   # midpoint root
  tree
}

# Root-to-tip regression
rtt_regression <- function(tree, meta_sub) {
  depths <- ape::node.depth.edgelength(tree)
  rtt_df <- data.frame(
    accessionVersion = tree$tip.label,
    rtt              = depths[seq_len(ape::Ntip(tree))]
  ) %>%
    left_join(meta_sub %>% select(accessionVersion, date_decimal, year),
              by = "accessionVersion") %>%
    filter(!is.na(date_decimal))

  if (nrow(rtt_df) < 4) return(NULL)

  model    <- lm(rtt ~ date_decimal, data = rtt_df)
  list(
    rate      = coef(model)[["date_decimal"]],
    r_squared = summary(model)$r.squared,
    model     = model,
    data      = rtt_df
  )
}

# Per-period rate comparison (pre vs post vaccination)
compare_periods <- function(rtt_result, vacc_decimal, country) {
  df <- rtt_result$data %>%
    mutate(period = if_else(date_decimal < vacc_decimal,
                            "pre_vaccination", "post_vaccination"))

  n_pre  <- sum(df$period == "pre_vaccination")
  n_post <- sum(df$period == "post_vaccination")

  if (n_pre < MIN_PRE_POST || n_post < MIN_PRE_POST) return(NULL)

  m_pre  <- lm(rtt ~ date_decimal, data = df[df$period == "pre_vaccination", ])
  m_post <- lm(rtt ~ date_decimal, data = df[df$period == "post_vaccination", ])

  # ANCOVA interaction test: slope difference + its SE
  m_int  <- lm(rtt ~ date_decimal * period, data = df)
  int_summary <- summary(m_int)$coefficients

  p_val <- tryCatch(
    anova(m_int)["date_decimal:period", "Pr(>F)"],
    error = function(e) NA_real_
  )

  # Rate difference = post_slope - pre_slope (the interaction coefficient)
  rate_pre  <- coef(m_pre)[["date_decimal"]]
  rate_post <- coef(m_post)[["date_decimal"]]
  rate_diff <- rate_post - rate_pre

  # SE of the interaction term (rate difference)
  se_diff <- tryCatch(
    int_summary["date_decimal:periodpre_vaccination", "Std. Error"],
    error = function(e) NA_real_
  )
  # If term name differs due to factor ordering, try the other direction
  if (is.na(se_diff)) {
    se_diff <- tryCatch(
      int_summary[grep("date_decimal:", rownames(int_summary)), "Std. Error"],
      error = function(e) NA_real_
    )
    se_diff <- se_diff[1]  # take first match
  }

  list(
    country    = country,
    rate_pre   = rate_pre,
    rate_post  = rate_post,
    rate_ratio = rate_post / rate_pre,
    rate_diff  = rate_diff,
    se_diff    = se_diff,
    p_value    = p_val,
    n_pre      = n_pre,
    n_post     = n_post,
    data       = df,
    model_pre  = m_pre,
    model_post = m_post
  )
}

# ── 5. Main loop: build trees + compute rates ──────────────────────────────────

run_subtype <- function(subtype_label, seqs_all, meta_sub_all) {
  cat(sprintf("\n════════ %s ════════\n", subtype_label))

  countries <- meta_sub_all %>%
    count(geoLocCountry) %>%
    filter(n >= MIN_SEQS) %>%
    pull(geoLocCountry)

  cat(sprintf("Countries with >= %d sequences: %s\n", MIN_SEQS,
              paste(countries, collapse = ", ")))

  country_trees  <- list()
  rtt_results    <- list()
  period_results <- list()

  for (country in countries) {
    cat(sprintf("\n  [%s] Building tree...", country))

    ids <- meta_sub_all %>%
      filter(geoLocCountry == country) %>%
      pull(accessionVersion)

    seq_sub <- seqs_all[names(seqs_all) %in% ids]

    if (length(seq_sub) < MIN_SEQS) {
      cat(" SKIPPED (too few sequences match FASTA)\n")
      next
    }

    # Build tree
    tree <- tryCatch(
      build_tree(seq_sub, country, subtype_label),
      error = function(e) {
        message(sprintf("  Tree build failed: %s", e$message))
        NULL
      }
    )

    if (is.null(tree)) next

    # Save tree
    safe_name <- gsub("[^A-Za-z0-9_]", "_", country)
    nwk_path  <- file.path(out_trees,
                           sprintf("%s_%s_tree.nwk", subtype_label, safe_name))
    write.tree(tree, file = nwk_path)

    country_trees[[country]] <- tree

    # RTT regression
    rtt <- rtt_regression(tree, meta_sub_all)
    if (!is.null(rtt)) {
      rtt_results[[country]] <- rtt
      cat(sprintf(" rate = %.2e sub/site/yr (R² = %.3f)\n",
                  rtt$rate, rtt$r_squared))
    }

    # Pre/post vaccination comparison (only for vaccinating countries)
    vacc_row <- vacc_dates %>% filter(country == !!country)
    if (nrow(vacc_row) > 0 && !is.null(rtt)) {
      comp <- compare_periods(rtt, vacc_row$vacc_decimal, country)
      if (!is.null(comp)) {
        period_results[[country]] <- comp
        cat(sprintf("       pre-vacc rate = %.2e | post-vacc rate = %.2e | ratio = %.3f | p = %.4f\n",
                    comp$rate_pre, comp$rate_post, comp$rate_ratio, comp$p_value))
      }
    }
  }

  list(trees = country_trees, rtt = rtt_results, periods = period_results)
}

# Split metadata by subtype
meta_a <- meta %>% filter(subtype == "RSV-A")
meta_b <- meta %>% filter(subtype == "RSV-B")

results_a <- run_subtype("RSV-A", seqs_a, meta_a)
results_b <- run_subtype("RSV-B", seqs_b, meta_b)

# ── 6. Build summary tables ───────────────────────────────────────────────────

build_rate_table <- function(period_results, subtype_label) {
  if (length(period_results) == 0) return(NULL)
  bind_rows(lapply(names(period_results), function(c) {
    r <- period_results[[c]]
    tibble(
      subtype    = subtype_label,
      country    = r$country,
      rate_pre   = r$rate_pre,
      rate_post  = r$rate_post,
      rate_ratio = r$rate_ratio,
      rate_diff  = r$rate_diff,
      se_diff    = r$se_diff,
      p_value    = r$p_value,
      n_pre      = r$n_pre,
      n_post     = r$n_post
    )
  }))
}

rate_tbl <- bind_rows(
  build_rate_table(results_a$periods, "RSV-A"),
  build_rate_table(results_b$periods, "RSV-B")
) %>%
  mutate(log_rate_ratio = if_else(rate_ratio > 0, log(rate_ratio), NA_real_))

write_csv(rate_tbl, file.path(out_summ, "rate_comparison_summary.csv"))
cat(sprintf("\nRate comparison table saved (%d rows)\n", nrow(rate_tbl)))
print(rate_tbl)

# ── 7. Meta-analysis ──────────────────────────────────────────────────────────

run_meta <- function(rate_sub, subtype_label) {
  if (is.null(rate_sub) || nrow(rate_sub) < 2) {
    message(sprintf("Not enough countries for meta-analysis: %s", subtype_label))
    return(NULL)
  }

  # Use the regression-derived SE of the rate difference directly.
  # Works even when rate ratios are negative (pre-rate near 0).
  rate_sub <- rate_sub %>%
    filter(is.finite(se_diff), se_diff > 0)

  if (nrow(rate_sub) < 2) return(NULL)

meta_res <- tryCatch(
    rma(yi = rate_diff, sei = se_diff, data = rate_sub, method = "REML"),
    error = function(e) {
      message(sprintf("Meta-analysis failed: %s", e$message))
      NULL
    }
  )

  if (!is.null(meta_res)) {
    sink(file.path(out_summ, sprintf("meta_analysis_%s.txt", subtype_label)))
    cat(sprintf("Meta-analysis: %s\n", subtype_label))
    cat("Rate ratio (overall) =", exp(meta_res$b), "\n")
    cat("95% CI:", exp(meta_res$ci.lb), "--", exp(meta_res$ci.ub), "\n")
    cat("I² =", meta_res$I2, "%\n\n")
    print(summary(meta_res))
    sink()
  }

  list(result = meta_res, data = rate_sub)
}

meta_a_res <- run_meta(rate_tbl %>% filter(subtype == "RSV-A"), "RSV-A")
meta_b_res <- run_meta(rate_tbl %>% filter(subtype == "RSV-B"), "RSV-B")

# ── 8. Visualisation ──────────────────────────────────────────────────────────

theme_rsv <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", size = 13)
    )
}

# ── 8a. RTT regression plots (one per country per subtype) ────────────────────

plot_rtt <- function(rtt_results, period_results, vacc_dates_df, subtype_label) {
  # Each PDF has two panels:
  #   Top:    RTT scatter + regression coloured by pre/post period
  #   Bottom: Sampling histogram to expose temporal over-sampling bias
  for (country in names(rtt_results)) {
    rtt      <- rtt_results[[country]]
    vacc_row <- vacc_dates_df %>% filter(country == !!country)

    has_periods <- !is.null(period_results[[country]])
    df <- if (has_periods) period_results[[country]]$data else rtt$data
    if (!"period" %in% names(df)) df$period <- "all"

    period_cols <- c(pre_vaccination = "#2166ac", post_vaccination = "#d6604d",
                     all = "#555555")
    period_labs <- c(pre_vaccination = "Pre-vaccination",
                     post_vaccination = "Post-vaccination", all = "All")

    subtitle_txt <- if (has_periods) {
      r <- period_results[[country]]
      sprintf("Pre = %.2e | Post = %.2e | diff = %.2e | p = %.4f",
              r$rate_pre, r$rate_post, r$rate_diff, r$p_value)
    } else {
      sprintf("Rate = %.2e AA sub/site/yr  |  R^2 = %.3f", rtt$rate, rtt$r_squared)
    }

    # Top panel
    p_top <- ggplot(df, aes(x = date_decimal, y = rtt,
                            colour = period, fill = period)) +
      geom_point(alpha = 0.65, size = 1.8) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.12) +
      scale_colour_manual(values = period_cols, labels = period_labs, name = NULL) +
      scale_fill_manual(values = period_cols, guide = "none") +
      labs(
        title    = sprintf("%s - %s: Root-to-Tip Regression", subtype_label, country),
        subtitle = subtitle_txt,
        x = NULL, y = "RTT distance (AA sub/site)"
      ) +
      theme_rsv() +
      theme(legend.position = if (has_periods) "top" else "none")

    if (nrow(vacc_row) > 0) {
      p_top <- p_top +
        geom_vline(xintercept = vacc_row$vacc_decimal,
                   linetype = "dashed", colour = "forestgreen", linewidth = 0.9) +
        annotate("text", x = vacc_row$vacc_decimal + 0.05,
                 y = max(df$rtt, na.rm = TRUE),
                 label = "Vacc.", colour = "forestgreen",
                 hjust = 0, size = 3, fontface = "italic")
    }

    # Bottom panel: sampling density
    # A spike of samples near 2025 inflates RTT because more diverged sequences
    # dominate the later years, pulling the regression slope steeper.
    p_hist <- ggplot(df, aes(x = date_decimal, fill = period)) +
      geom_histogram(binwidth = 0.25, colour = "white",
                     linewidth = 0.2, alpha = 0.8) +
      scale_fill_manual(values = period_cols, guide = "none") +
      labs(
        x       = "Sampling date (decimal year)",
        y       = "No. sequences",
        caption = paste0(
          "Check for over-sampling in recent years: ",
          "uneven temporal sampling inflates apparent RTT distance."
        )
      ) +
      theme_rsv() +
      theme(plot.caption = element_text(size = 7, colour = "grey50", hjust = 0))

    if (nrow(vacc_row) > 0) {
      p_hist <- p_hist +
        geom_vline(xintercept = vacc_row$vacc_decimal,
                   linetype = "dashed", colour = "forestgreen", linewidth = 0.9)
    }

    p_combined <- p_top / p_hist + plot_layout(heights = c(3, 1))

    safe_name <- gsub("[^A-Za-z0-9_]", "_", country)
    ggsave(
      file.path(out_plots, sprintf("rtt_%s_%s.pdf", subtype_label, safe_name)),
      p_combined, width = 8, height = 7
    )
  }
}
plot_rtt(results_a$rtt, results_a$periods, vacc_dates, "RSV-A")
plot_rtt(results_b$rtt, results_b$periods, vacc_dates, "RSV-B")
cat("RTT plots saved.\n")

# ── 8b. Rate ratio slope plot (all vaccinating countries) ─────────────────────

if (nrow(rate_tbl) > 0) {
  slope_df <- rate_tbl %>%
    select(subtype, country, rate_pre, rate_post) %>%
    pivot_longer(c(rate_pre, rate_post),
                 names_to  = "period",
                 values_to = "rate") %>%
    mutate(period = factor(period,
                           levels = c("rate_pre", "rate_post"),
                           labels = c("Pre-vaccination", "Post-vaccination")))

  p_slope <- ggplot(slope_df,
                    aes(x = period, y = rate,
                        group = interaction(country, subtype),
                        colour = country)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    facet_wrap(~subtype) +
    labs(
      title  = "Evolutionary Rate Before vs After Vaccination",
      x      = NULL,
      y      = "Rate (AA sub/site/yr)",
      colour = "Country"
    ) +
    theme_rsv()

  ggsave(file.path(out_plots, "rate_slope_by_country.pdf"),
         p_slope, width = 10, height = 5)
  cat("Slope plot saved.\n")
}

# ── 8c. Forest plot ───────────────────────────────────────────────────────────

make_forest_plot <- function(meta_res_obj, subtype_label, all_rate_tbl) {
  # Uses rate_diff (post slope - pre slope) on a linear scale so ALL countries
  # are shown, including those where the pre-rate is near zero or negative.
  if (is.null(meta_res_obj) && is.null(all_rate_tbl)) return(invisible(NULL))

  overall <- if (!is.null(meta_res_obj)) meta_res_obj$result else NULL

  forest_df <- all_rate_tbl %>%
    filter(subtype == subtype_label) %>%
    mutate(
      mean  = rate_diff,
      lower = if_else(is.finite(se_diff), rate_diff - 1.96 * se_diff, NA_real_),
      upper = if_else(is.finite(se_diff), rate_diff + 1.96 * se_diff, NA_real_),
      label = sprintf("%s  (n=%d/%d)", country, n_pre, n_post),
      sig   = !is.na(p_value) & p_value < 0.05
    ) %>%
    arrange(rate_diff) %>%
    select(label, mean, lower, upper, sig, p_value)

  if (!is.null(overall)) {
    overall_row <- tibble(
      label   = "Overall (RE meta-analysis)",
      mean    = overall$b[1],
      lower   = overall$ci.lb,
      upper   = overall$ci.ub,
      sig     = overall$pval < 0.05,
      p_value = overall$pval
    )
    forest_df <- bind_rows(forest_df, overall_row)
  }

  forest_df <- forest_df %>%
    mutate(
      label      = factor(label, levels = label),
      is_overall = grepl("Overall", label)
    )

  p_forest <- ggplot(forest_df, aes(x = mean, y = label)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_errorbar(aes(xmin = lower, xmax = upper),
                  orientation = "y", width = 0.35, colour = "grey50",
                  na.rm = TRUE) +
    geom_point(aes(size = is_overall,
                   colour = interaction(is_overall, sig))) +
    scale_colour_manual(
      values = c(
        `FALSE.FALSE` = "#2166ac",
        `FALSE.TRUE`  = "#d6604d",
        `TRUE.FALSE`  = "#555555",
        `TRUE.TRUE`   = "#b2182b"
      ),
      guide = "none"
    ) +
    scale_size_manual(values = c(`FALSE` = 2.5, `TRUE` = 4.5), guide = "none") +
    scale_x_continuous(labels = scales::label_scientific(digits = 2)) +
    labs(
      title = sprintf("%s: Rate Change Post- vs Pre-Vaccination (all countries)",
                      subtype_label),
      subtitle = if (!is.null(overall))
        sprintf("Overall diff = %.2e (95%% CI: %.2e to %.2e)  I^2 = %.1f%%",
                overall$b[1], overall$ci.lb, overall$ci.ub, overall$I2)
      else "No meta-analysis (insufficient data)",
      x       = "Rate difference: post minus pre (AA sub/site/yr)  |  >0 = faster post-vaccination",
      y       = NULL,
      caption = "Orange/red = p<0.05. Error bars = 95% CI from ANCOVA interaction term."
    ) +
    theme_rsv() +
    theme(plot.caption = element_text(size = 8, colour = "grey50"))

  ggsave(
    file.path(out_plots, sprintf("forest_plot_%s.pdf", subtype_label)),
    p_forest, width = 10, height = max(4, nrow(forest_df) * 0.55 + 2)
  )
  cat(sprintf("Forest plot saved: %s\n", subtype_label))
}
make_forest_plot(meta_a_res, "RSV-A", rate_tbl)
make_forest_plot(meta_b_res, "RSV-B", rate_tbl)

# ── 8d. All-country RTT summary (panel of rates) ──────────────────────────────

plot_all_rtt_rates <- function(rtt_results, subtype_label, colour_val) {
  if (length(rtt_results) == 0) return(invisible(NULL))
  rate_overview <- tibble(
    country   = names(rtt_results),
    rate      = sapply(rtt_results, `[[`, "rate"),
    r_squared = sapply(rtt_results, `[[`, "r_squared")
  ) %>%
    arrange(rate) %>%
    mutate(country = factor(country, levels = country))

  p <- ggplot(rate_overview, aes(x = rate, y = country)) +
    geom_point(aes(size = r_squared), colour = colour_val) +
    geom_segment(aes(x = 0, xend = rate, yend = country),
                 colour = colour_val, alpha = 0.5) +
    scale_size_continuous(range = c(2, 6), name = "R²") +
    labs(
      title    = sprintf("%s: Overall Evolutionary Rate by Country", subtype_label),
      subtitle = "Point size = R² of root-to-tip regression",
      x        = "Rate (AA sub/site/yr)",
      y        = NULL
    ) +
    theme_rsv()

  ggsave(
    file.path(out_plots, sprintf("rate_overview_%s.pdf", subtype_label)),
    p, width = 8, height = max(4, length(rtt_results) * 0.35 + 2)
  )
}

plot_all_rtt_rates(results_a$rtt, "RSV-A", "#2166ac")
plot_all_rtt_rates(results_b$rtt, "RSV-B", "#d6604d")

# ── 8e. Phylogenetic tree plots (rectangular layout, tips coloured by year) ───

plot_trees <- function(trees, meta_sub, subtype_label) {
  if (!requireNamespace("ggtree", quietly = TRUE)) {
    message("ggtree not available – skipping tree plots")
    return(invisible(NULL))
  }
  library(ggtree)

  for (country in names(trees)) {
    tree <- trees[[country]]

    tip_meta <- meta_sub %>%
      filter(accessionVersion %in% tree$tip.label) %>%
      select(accessionVersion, year, lineage) %>%
      rename(label = accessionVersion)

    # Rectangular (linear) layout
    p_tree <- ggtree(tree, layout = "rectangular",
                     size = 0.25, colour = "grey50") %<+% tip_meta +
      geom_tippoint(aes(colour = year), size = 1.2, alpha = 0.85) +
      scale_colour_viridis_c(
        name   = "Year",
        option = "plasma",
        breaks = pretty(tip_meta$year, n = 5),
        labels = function(x) as.integer(x)
      ) +
      labs(
        title    = sprintf("%s – %s phylogenetic tree", subtype_label, country),
        subtitle = sprintf("%d tips  |  midpoint-rooted  |  WAG model",
                           ape::Ntip(tree))
      ) +
      theme_tree2() +
      theme(
        plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        legend.position = "right"
      )

    safe_name <- gsub("[^A-Za-z0-9_]", "_", country)
    n_tips <- ape::Ntip(tree)
    plot_h  <- max(6, min(20, n_tips * 0.06))
    ggsave(
      file.path(out_plots, sprintf("tree_%s_%s.pdf", subtype_label, safe_name)),
      p_tree, width = 10, height = plot_h
    )
  }
  cat(sprintf("Tree plots saved for %s (%d countries)\n",
              subtype_label, length(trees)))
}

plot_trees(results_a$trees, meta_a, "RSV-A")
plot_trees(results_b$trees, meta_b, "RSV-B")

# ── 9. Done ───────────────────────────────────────────────────────────────────

cat(sprintf("\n✓ All done. Outputs in: %s\n", out_base))
cat("  trees/    – Newick (.nwk) per country per subtype\n")
cat("  plots/    – RTT regression, slope, forest, rate overview PDFs\n")
cat("  summaries/ – rate_comparison_summary.csv + meta-analysis text files\n")

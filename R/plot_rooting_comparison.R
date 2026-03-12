# plot_rooting_comparison.R
# Visualises how much root placement (outgroup vs midpoint) affects RTT slope
# estimates. Three panels:
#   1. Scatter: outgroup rate vs midpoint rate (identity line = no difference)
#   2. Dot plot: difference in rate per country, sorted by magnitude
#   3. R² comparison: does rooting change how well the clock fits?

library(tidyverse)
library(here)
library(patchwork)

cmp <- read_csv(here("results", "phylo_trees_full", "summaries", "rooting_comparison.csv"),
                show_col_types = FALSE)

subtype_cols <- c("RSV-A" = "#2166ac", "RSV-B" = "#d6604d")

# ── 1. Scatter: outgroup rate vs midpoint rate ────────────────────────────────
lims <- range(c(cmp$rate_outgroup, cmp$rate_midpoint), na.rm = TRUE)

# Only label countries with the largest absolute deviation from the identity line
top_dev <- cmp %>%
  mutate(dev = abs(rate_outgroup - rate_midpoint)) %>%
  slice_max(dev, n = 15)

p1 <- ggplot(cmp, aes(x = rate_midpoint, y = rate_outgroup, colour = subtype)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 2.5, alpha = 0.8) +
  ggrepel::geom_text_repel(
    data = top_dev,
    aes(label = paste0(country, "\n(", subtype, ")")),
    size = 2.3, max.overlaps = 30, show.legend = FALSE,
    segment.colour = "grey60", segment.size = 0.3
  ) +
  scale_colour_manual(values = subtype_cols, name = NULL) +
  coord_equal(xlim = lims, ylim = lims) +
  labs(
    title = "Outgroup vs midpoint rooting: RTT slope",
    x = "Midpoint-rooted rate (AA sub/site/yr)",
    y = "Outgroup-rooted rate (AA sub/site/yr)"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# ── 2. Dot plot: rate difference (outgroup − midpoint) ───────────────────────
cmp2 <- cmp %>%
  arrange(subtype, rate_diff_root) %>%
  mutate(country_label = fct_inorder(paste0(country, " (", subtype, ")")))

p2 <- ggplot(cmp2, aes(x = rate_diff_root, y = country_label, colour = subtype)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_segment(aes(x = 0, xend = rate_diff_root,
                   y = country_label, yend = country_label),
               linewidth = 0.4, alpha = 0.5) +
  scale_colour_manual(values = subtype_cols, name = NULL) +
  labs(
    title = "Rate difference by rooting method",
    subtitle = "Outgroup-rooted minus midpoint-rooted (AA sub/site/yr)",
    x = "Rate difference", y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(face = "bold"),
        legend.position = "none")

# ── 3. R² comparison ──────────────────────────────────────────────────────────
r2_long <- cmp %>%
  select(country, subtype, r2_outgroup, r2_midpoint) %>%
  pivot_longer(c(r2_outgroup, r2_midpoint),
               names_to = "method", values_to = "r2") %>%
  mutate(method = recode(method,
                         r2_outgroup = "Outgroup",
                         r2_midpoint = "Midpoint"))

p3 <- ggplot(r2_long, aes(x = method, y = r2, colour = subtype,
                           group = paste(country, subtype))) +
  geom_line(alpha = 0.35, linewidth = 0.4) +
  geom_point(size = 1.8, alpha = 0.7) +
  scale_colour_manual(values = subtype_cols, name = NULL) +
  labs(
    title = "Clock fit (R²) by rooting method",
    x = NULL, y = "R² (RTT regression)"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# ── Combine and save ──────────────────────────────────────────────────────────
out_dir <- here("results", "phylo_trees_full", "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p_all <- (p1 | p3) / p2 +
  plot_annotation(
    title   = "Effect of root placement on RTT rate estimates",
    subtitle = "Points off the dashed line (panel 1) or far from zero (panel 2) indicate\nrooting materially changes rate estimates for that country.",
    theme = theme(plot.title    = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 9, colour = "grey40"))
  )

ggsave(file.path(out_dir, "rooting_comparison.pdf"),
       p_all, width = 13, height = 14)
ggsave(file.path(out_dir, "rooting_comparison.png"),
       p_all, width = 13, height = 14, dpi = 150)

cat("Saved rooting_comparison.pdf/.png\n")

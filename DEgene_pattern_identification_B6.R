# ===================================================================
# GENE PATTERN IDENTIFICATION - WITH MAX vs POST COMPARISONS
# Classify genes by their developmental trajectories
# ===================================================================

library(tidyverse)

cat("===== IDENTIFYING GENE EXPRESSION PATTERNS =====\n\n")

# ===================================================================
# PART 1: MAIN EFFECTS - DOMAIN DIFFERENCES (NO INTERACTION)
# These genes show different expression levels across domains
# but maintain parallel developmental trajectories
# ===================================================================

cat("PART 1: DOMAIN MAIN EFFECTS (Parallel Trajectories)\n")
cat("=" %>% rep(60) %>% paste(collapse=""), "\n\n")

# ----- MAX vs PM Comparisons -----

genes_MAX_higher_than_PM <- b6_domain_main_effects_clean %>%
  filter(Comparison == "MAX vs PM", log2FoldChange > 0) %>%
  pull(gene_id)

genes_MAX_lower_than_PM <- b6_domain_main_effects_clean %>%
  filter(Comparison == "MAX vs PM", log2FoldChange < 0) %>%
  pull(gene_id)

cat("MAX vs PM:\n")
cat("  MAX higher than PM:", length(genes_MAX_higher_than_PM), "genes\n")
cat("  MAX lower than PM:", length(genes_MAX_lower_than_PM), "genes\n\n")

# ----- POST vs PM Comparisons -----

genes_POST_higher_than_PM <- b6_domain_main_effects_clean %>%
  filter(Comparison == "POST vs PM", log2FoldChange > 0) %>%
  pull(gene_id)

genes_POST_lower_than_PM <- b6_domain_main_effects_clean %>%
  filter(Comparison == "POST vs PM", log2FoldChange < 0) %>%
  pull(gene_id)

cat("POST vs PM:\n")
cat("  POST higher than PM:", length(genes_POST_higher_than_PM), "genes\n")
cat("  POST lower than PM:", length(genes_POST_lower_than_PM), "genes\n\n")

# ----- POST vs MAX Comparisons -----

genes_POST_higher_than_MAX <- b6_domain_main_effects_clean %>%
  filter(Comparison == "POST vs MAX", log2FoldChange > 0) %>%
  pull(gene_id)

genes_POST_lower_than_MAX <- b6_domain_main_effects_clean %>%
  filter(Comparison == "POST vs MAX", log2FoldChange < 0) %>%
  pull(gene_id)

cat("POST vs MAX:\n")
cat("  POST higher than MAX:", length(genes_POST_higher_than_MAX), "genes\n")
cat("  POST lower than MAX:", length(genes_POST_lower_than_MAX), "genes\n\n")

# ----- Complex Domain Patterns (Combinations) -----

# Genes highest in MAX (compared to both PM and POST)
genes_MAX_highest <- intersect(genes_MAX_higher_than_PM, 
                               genes_POST_lower_than_MAX)

# Genes highest in POST (compared to both PM and MAX)
genes_POST_highest <- intersect(genes_POST_higher_than_PM, 
                                genes_POST_higher_than_MAX)

# Genes highest in PM (compared to both MAX and POST)
genes_PM_highest <- intersect(genes_MAX_lower_than_PM, 
                              genes_POST_lower_than_PM)

# Genes lowest in MAX
genes_MAX_lowest <- intersect(genes_MAX_lower_than_PM, 
                              genes_POST_higher_than_MAX)

# Genes lowest in POST
genes_POST_lowest <- intersect(genes_POST_lower_than_PM, 
                               genes_POST_lower_than_MAX)

# Genes lowest in PM
genes_PM_lowest <- intersect(genes_MAX_higher_than_PM, 
                             genes_POST_higher_than_PM)

cat("Complex Domain Patterns:\n")
cat("  Highest in MAX:", length(genes_MAX_highest), "genes\n")
cat("  Highest in POST:", length(genes_POST_highest), "genes\n")
cat("  Highest in PM:", length(genes_PM_highest), "genes\n")
cat("  Lowest in MAX:", length(genes_MAX_lowest), "genes\n")
cat("  Lowest in POST:", length(genes_POST_lowest), "genes\n")
cat("  Lowest in PM:", length(genes_PM_lowest), "genes\n\n")

# ----- Monotonic Gradients -----

# A→P decreasing (PM > MAX > POST)
genes_AP_decreasing <- intersect(
  genes_MAX_lower_than_PM,
  genes_POST_lower_than_MAX
)

# A→P increasing (POST > MAX > PM)
genes_AP_increasing <- intersect(
  genes_POST_higher_than_PM,
  genes_POST_higher_than_MAX
) %>% intersect(genes_MAX_higher_than_PM)

cat("Monotonic A-P Gradients:\n")
cat("  A→P decreasing (PM > MAX > POST):", length(genes_AP_decreasing), "genes\n")
cat("  A→P increasing (POST > MAX > PM):", length(genes_AP_increasing), "genes\n\n")

# ===================================================================
# PART 2: MAIN EFFECTS - DEVELOPMENTAL CHANGES (NO INTERACTION)
# ===================================================================

cat("\nPART 2: STAGE MAIN EFFECTS (Common Trajectories)\n")
cat("=" |> rep(60) |> paste(collapse=""), "\n\n")

genes_increasing_with_development <- b6_stage_main_effects_clean |>
  filter(Direction == "Increase") |>
  pull(gene_id)

genes_decreasing_with_development <- b6_stage_main_effects_clean |>
  filter(Direction == "Decrease") |>
  pull(gene_id)

cat("Developmental Changes (all domains together):\n")
cat("  Increasing:", length(genes_increasing_with_development), "genes\n")
cat("  Decreasing:", length(genes_decreasing_with_development), "genes\n\n")

# Strong changes
strong_increasers <- b6_stage_main_effects_clean |>
  filter(Direction == "Increase", abs(log2FoldChange) > 1) |>
  pull(gene_id)

strong_decreasers <- b6_stage_main_effects_clean |>
  filter(Direction == "Decrease", abs(log2FoldChange) > 1) |>
  pull(gene_id)

cat("Strong changes (|log2FC| > 1):\n")
cat("  Strong increasers:", length(strong_increasers), "genes\n")
cat("  Strong decreasers:", length(strong_decreasers), "genes\n\n")

# ===================================================================
# PART 3: INTERACTION EFFECTS - DIVERGING TRAJECTORIES
# NOTE: Using b6_domain_slopes (UNFILTERED) here because we need BOTH
# slopes to classify the pattern, even if one slope isn't individually
# significant. The slope DIFFERENCE is what matters for diverging genes.
# ===================================================================

cat("\nPART 3: DIVERGING TRAJECTORIES\n")
cat("=" |> rep(60) |> paste(collapse=""), "\n\n")

# ----- MAX vs PM Divergence -----

genes_MAX_diverging_from_PM <- b6_slope_differences |>
  filter(Comparison == "MAX vs PM slope difference",
         padj < 0.05, log2FoldChange > 0.5) |>
  pull(gene_id)

genes_MAX_diverging_details_PM <- b6_domain_slopes |>  # USE UNFILTERED
  filter(gene_id %in% genes_MAX_diverging_from_PM) |>
  select(gene_id, symbol, Domain, log2FoldChange, padj) |>
  pivot_wider(names_from = Domain, values_from = c(log2FoldChange, padj), 
              names_glue = "{.value}_{Domain}") |>
  mutate(
    pattern = case_when(
      !is.na(log2FoldChange_MAX) & !is.na(log2FoldChange_PM) & 
        log2FoldChange_MAX > 0 & log2FoldChange_PM >= 0 ~ "MAX increases faster",
      !is.na(log2FoldChange_MAX) & !is.na(log2FoldChange_PM) & 
        log2FoldChange_MAX > 0 & log2FoldChange_PM < 0 ~ "MAX increases, PM decreases",
      TRUE ~ "Other"
    )
  )

# Diagnostic output
if(nrow(genes_MAX_diverging_details_PM) > 0) {
  cat("  Diagnostic - Pattern breakdown:\n")
  cat("    ", table(genes_MAX_diverging_details_PM$pattern), "\n")
}

genes_MAX_increases_faster_than_PM <- genes_MAX_diverging_details_PM |>
  filter(pattern == "MAX increases faster") |> pull(gene_id)

genes_MAX_up_PM_down <- genes_MAX_diverging_details_PM |>
  filter(pattern == "MAX increases, PM decreases") |> pull(gene_id)

cat("MAX diverging from PM:\n")
cat("  Total:", length(genes_MAX_diverging_from_PM), "genes\n")
cat("  MAX increases faster:", length(genes_MAX_increases_faster_than_PM), "genes\n")
cat("  MAX up, PM down:", length(genes_MAX_up_PM_down), "genes\n\n")

# ----- POST vs PM Divergence -----

genes_POST_diverging_from_PM <- b6_slope_differences |>
  filter(Comparison == "POST vs PM slope difference",
         padj < 0.05, log2FoldChange > 0.5) |>
  pull(gene_id)

genes_POST_diverging_details_PM <- b6_domain_slopes |>  # USE UNFILTERED
  filter(gene_id %in% genes_POST_diverging_from_PM) |>
  select(gene_id, symbol, Domain, log2FoldChange, padj) |>
  pivot_wider(names_from = Domain, values_from = c(log2FoldChange, padj), 
              names_glue = "{.value}_{Domain}") |>
  mutate(
    pattern = case_when(
      !is.na(log2FoldChange_POST) & !is.na(log2FoldChange_PM) & 
        log2FoldChange_POST > 0 & log2FoldChange_PM >= 0 ~ "POST increases faster",
      !is.na(log2FoldChange_POST) & !is.na(log2FoldChange_PM) & 
        log2FoldChange_POST > 0 & log2FoldChange_PM < 0 ~ "POST increases, PM decreases",
      TRUE ~ "Other"
    )
  )

genes_POST_increases_faster_than_PM <- genes_POST_diverging_details_PM |>
  filter(pattern == "POST increases faster") |> pull(gene_id)

genes_POST_up_PM_down <- genes_POST_diverging_details_PM |>
  filter(pattern == "POST increases, PM decreases") |> pull(gene_id)

cat("POST diverging from PM:\n")
cat("  Total:", length(genes_POST_diverging_from_PM), "genes\n")
cat("  POST increases faster:", length(genes_POST_increases_faster_than_PM), "genes\n")
cat("  POST up, PM down:", length(genes_POST_up_PM_down), "genes\n\n")

# ----- POST vs MAX Divergence (NEW) -----

genes_POST_diverging_from_MAX <- b6_slope_differences |>
  filter(Comparison == "POST vs MAX slope difference",
         padj < 0.05, log2FoldChange > 0.5) |>
  pull(gene_id)

genes_POST_diverging_details_MAX <- b6_domain_slopes |>  # USE UNFILTERED
  filter(gene_id %in% genes_POST_diverging_from_MAX) |>
  select(gene_id, symbol, Domain, log2FoldChange, padj) |>
  pivot_wider(names_from = Domain, values_from = c(log2FoldChange, padj), 
              names_glue = "{.value}_{Domain}") |>
  mutate(
    pattern = case_when(
      !is.na(log2FoldChange_POST) & !is.na(log2FoldChange_MAX) & 
        log2FoldChange_POST > 0 & log2FoldChange_MAX >= 0 ~ "POST increases faster",
      !is.na(log2FoldChange_POST) & !is.na(log2FoldChange_MAX) & 
        log2FoldChange_POST > 0 & log2FoldChange_MAX < 0 ~ "POST increases, MAX decreases",
      TRUE ~ "Other"
    )
  )

genes_POST_increases_faster_than_MAX <- genes_POST_diverging_details_MAX |>
  filter(pattern == "POST increases faster") |> pull(gene_id)

genes_POST_up_MAX_down <- genes_POST_diverging_details_MAX |>
  filter(pattern == "POST increases, MAX decreases") |> pull(gene_id)

cat("POST diverging from MAX:\n")
cat("  Total:", length(genes_POST_diverging_from_MAX), "genes\n")
cat("  POST increases faster:", length(genes_POST_increases_faster_than_MAX), "genes\n")
cat("  POST up, MAX down:", length(genes_POST_up_MAX_down), "genes\n\n")

# ===================================================================
# PART 4: INTERACTION EFFECTS - CONVERGING TRAJECTORIES
# ===================================================================

cat("\nPART 4: CONVERGING TRAJECTORIES\n")
cat("=" %>% rep(60) %>% paste(collapse=""), "\n\n")

# ----- MAX vs PM Convergence -----

genes_MAX_converging_to_PM <- b6_slope_differences %>%
  filter(Comparison == "MAX vs PM slope difference",
         padj < 0.05, log2FoldChange < -0.5) %>%
  pull(gene_id)

genes_MAX_converging_details_PM <- b6_domain_slopes_clean %>%
  filter(gene_id %in% genes_MAX_converging_to_PM) %>%
  select(gene_id, symbol, Domain, log2FoldChange, padj) %>%
  pivot_wider(names_from = Domain, values_from = c(log2FoldChange, padj), 
              names_glue = "{.value}_{Domain}") %>%
  mutate(
    pattern = case_when(
      !is.na(log2FoldChange_MAX) & !is.na(log2FoldChange_PM) & 
        log2FoldChange_MAX < 0 & log2FoldChange_PM <= 0 ~ "MAX decreases faster",
      !is.na(log2FoldChange_MAX) & !is.na(log2FoldChange_PM) & 
        log2FoldChange_MAX < 0 & log2FoldChange_PM > 0 ~ "MAX decreases, PM increases",
      TRUE ~ "Other"
    )
  )

genes_MAX_decreases_faster_than_PM <- genes_MAX_converging_details_PM %>%
  filter(pattern == "MAX decreases faster") %>% pull(gene_id)

genes_MAX_down_PM_up <- genes_MAX_converging_details_PM %>%
  filter(pattern == "MAX decreases, PM increases") %>% pull(gene_id)

cat("MAX converging to PM:\n")
cat("  Total:", length(genes_MAX_converging_to_PM), "genes\n")
cat("  MAX decreases faster:", length(genes_MAX_decreases_faster_than_PM), "genes\n")
cat("  MAX down, PM up:", length(genes_MAX_down_PM_up), "genes\n\n")

# ----- POST vs PM Convergence -----

genes_POST_converging_to_PM <- b6_slope_differences %>%
  filter(Comparison == "POST vs PM slope difference",
         padj < 0.05, log2FoldChange < -0.5) %>%
  pull(gene_id)

genes_POST_converging_details_PM <- b6_domain_slopes_clean %>%
  filter(gene_id %in% genes_POST_converging_to_PM) %>%
  select(gene_id, symbol, Domain, log2FoldChange, padj) %>%
  pivot_wider(names_from = Domain, values_from = c(log2FoldChange, padj), 
              names_glue = "{.value}_{Domain}") %>%
  mutate(
    pattern = case_when(
      !is.na(log2FoldChange_POST) & !is.na(log2FoldChange_PM) & 
        log2FoldChange_POST < 0 & log2FoldChange_PM <= 0 ~ "POST decreases faster",
      !is.na(log2FoldChange_POST) & !is.na(log2FoldChange_PM) & 
        log2FoldChange_POST < 0 & log2FoldChange_PM > 0 ~ "POST decreases, PM increases",
      TRUE ~ "Other"
    )
  )

genes_POST_decreases_faster_than_PM <- genes_POST_converging_details_PM %>%
  filter(pattern == "POST decreases faster") %>% pull(gene_id)

genes_POST_down_PM_up <- genes_POST_converging_details_PM %>%
  filter(pattern == "POST decreases, PM increases") %>% pull(gene_id)

cat("POST converging to PM:\n")
cat("  Total:", length(genes_POST_converging_to_PM), "genes\n")
cat("  POST decreases faster:", length(genes_POST_decreases_faster_than_PM), "genes\n")
cat("  POST down, PM up:", length(genes_POST_down_PM_up), "genes\n\n")

# ----- POST vs MAX Convergence (NEW) -----

genes_POST_converging_to_MAX <- b6_slope_differences %>%
  filter(Comparison == "POST vs MAX slope difference",
         padj < 0.05, log2FoldChange < -0.5) %>%
  pull(gene_id)

genes_POST_converging_details_MAX <- b6_domain_slopes_clean %>%
  filter(gene_id %in% genes_POST_converging_to_MAX) %>%
  select(gene_id, symbol, Domain, log2FoldChange, padj) %>%
  pivot_wider(names_from = Domain, values_from = c(log2FoldChange, padj), 
              names_glue = "{.value}_{Domain}") %>%
  mutate(
    pattern = case_when(
      !is.na(log2FoldChange_POST) & !is.na(log2FoldChange_MAX) & 
        log2FoldChange_POST < 0 & log2FoldChange_MAX <= 0 ~ "POST decreases faster",
      !is.na(log2FoldChange_POST) & !is.na(log2FoldChange_MAX) & 
        log2FoldChange_POST < 0 & log2FoldChange_MAX > 0 ~ "POST decreases, MAX increases",
      TRUE ~ "Other"
    )
  )

genes_POST_decreases_faster_than_MAX <- genes_POST_converging_details_MAX %>%
  filter(pattern == "POST decreases faster") %>% pull(gene_id)

genes_POST_down_MAX_up <- genes_POST_converging_details_MAX %>%
  filter(pattern == "POST decreases, MAX increases") %>% pull(gene_id)

cat("POST converging to MAX:\n")
cat("  Total:", length(genes_POST_converging_to_MAX), "genes\n")
cat("  POST decreases faster:", length(genes_POST_decreases_faster_than_MAX), "genes\n")
cat("  POST down, MAX up:", length(genes_POST_down_MAX_up), "genes\n\n")

# ===================================================================
# PART 5: DOMAIN-SPECIFIC DYNAMICS (DETAILED)
# ===================================================================

cat("\nPART 5: DOMAIN-SPECIFIC DYNAMICS\n")
cat("=" |> rep(60) |> paste(collapse=""), "\n\n")

domain_specific_detailed <- b6_domain_slopes_clean |>
  select(gene_id, symbol, Domain, log2FoldChange, padj) |>
  # Ensure unique rows per gene-domain combination
  distinct(gene_id, symbol, Domain, .keep_all = TRUE) |>
  pivot_wider(
    id_cols = c(gene_id, symbol),
    names_from = Domain,
    values_from = c(log2FoldChange, padj),
    names_glue = "{.value}_{Domain}"
  ) |>
  mutate(
    # Check if dynamic (significant slope with meaningful effect)
    PM_dynamic = !is.na(padj_PM) & padj_PM < 0.05 & abs(log2FoldChange_PM) > 0.5,
    MAX_dynamic = !is.na(padj_MAX) & padj_MAX < 0.05 & abs(log2FoldChange_MAX) > 0.5,
    POST_dynamic = !is.na(padj_POST) & padj_POST < 0.05 & abs(log2FoldChange_POST) > 0.5,
    
    # Get slope values (will be NA if not present)
    PM_slope = log2FoldChange_PM,
    MAX_slope = log2FoldChange_MAX,
    POST_slope = log2FoldChange_POST,
    
    # Count dynamic domains
    n_dynamic = PM_dynamic + MAX_dynamic + POST_dynamic,
    
    # Classify pattern
    pattern = case_when(
      n_dynamic == 0 ~ "None dynamic",
      n_dynamic == 3 ~ "All dynamic",
      PM_dynamic & !MAX_dynamic & !POST_dynamic ~ "PM only",
      !PM_dynamic & MAX_dynamic & !POST_dynamic ~ "MAX only",
      !PM_dynamic & !MAX_dynamic & POST_dynamic ~ "POST only",
      PM_dynamic & MAX_dynamic & !POST_dynamic ~ "PM & MAX",
      PM_dynamic & !MAX_dynamic & POST_dynamic ~ "PM & POST",
      !PM_dynamic & MAX_dynamic & POST_dynamic ~ "MAX & POST",
      TRUE ~ "Other"
    )
  )

# Extract gene lists
genes_PM_only_dynamic <- domain_specific_detailed |>
  filter(pattern == "PM only") |> pull(gene_id)

genes_MAX_only_dynamic <- domain_specific_detailed |>
  filter(pattern == "MAX only") |> pull(gene_id)

genes_POST_only_dynamic <- domain_specific_detailed |>
  filter(pattern == "POST only") |> pull(gene_id)

genes_PM_MAX_dynamic <- domain_specific_detailed |>
  filter(pattern == "PM & MAX") |> pull(gene_id)

genes_PM_POST_dynamic <- domain_specific_detailed |>
  filter(pattern == "PM & POST") |> pull(gene_id)

genes_MAX_POST_dynamic <- domain_specific_detailed |>
  filter(pattern == "MAX & POST") |> pull(gene_id)

cat("Single domain dynamics:\n")
cat("  PM only:", length(genes_PM_only_dynamic), "genes\n")
cat("  MAX only:", length(genes_MAX_only_dynamic), "genes\n")
cat("  POST only:", length(genes_POST_only_dynamic), "genes\n\n")

cat("Two domain dynamics:\n")
cat("  PM & MAX:", length(genes_PM_MAX_dynamic), "genes\n")
cat("  PM & POST:", length(genes_PM_POST_dynamic), "genes\n")
cat("  MAX & POST:", length(genes_MAX_POST_dynamic), "genes\n\n")

# Directional subgroups
genes_PM_only_increasing <- domain_specific_detailed |>
  filter(pattern == "PM only", !is.na(PM_slope), PM_slope > 0.5) |> 
  pull(gene_id)

genes_PM_only_decreasing <- domain_specific_detailed |>
  filter(pattern == "PM only", !is.na(PM_slope), PM_slope < -0.5) |> 
  pull(gene_id)

genes_MAX_only_increasing <- domain_specific_detailed |>
  filter(pattern == "MAX only", !is.na(MAX_slope), MAX_slope > 0.5) |> 
  pull(gene_id)

genes_MAX_only_decreasing <- domain_specific_detailed |>
  filter(pattern == "MAX only", !is.na(MAX_slope), MAX_slope < -0.5) |> 
  pull(gene_id)

genes_POST_only_increasing <- domain_specific_detailed |>
  filter(pattern == "POST only", !is.na(POST_slope), POST_slope > 0.5) |> 
  pull(gene_id)

genes_POST_only_decreasing <- domain_specific_detailed |>
  filter(pattern == "POST only", !is.na(POST_slope), POST_slope < -0.5) |> 
  pull(gene_id)

cat("Directional single-domain dynamics:\n")
cat("  PM only ↑:", length(genes_PM_only_increasing), "genes\n")
cat("  PM only ↓:", length(genes_PM_only_decreasing), "genes\n")
cat("  MAX only ↑:", length(genes_MAX_only_increasing), "genes\n")
cat("  MAX only ↓:", length(genes_MAX_only_decreasing), "genes\n")
cat("  POST only ↑:", length(genes_POST_only_increasing), "genes\n")
cat("  POST only ↓:", length(genes_POST_only_decreasing), "genes\n\n")

# ===================================================================
# PART 6: COMPLEX PATTERNS
# ===================================================================

cat("\nPART 6: COMPLEX PATTERNS\n")
cat("=" |> rep(60) |> paste(collapse=""), "\n\n")

# Opposite slopes
genes_opposite_slopes <- b6_domain_slopes_clean |>
  group_by(gene_id) |>
  summarize(
    has_positive = any(padj < 0.05 & log2FoldChange > 0.5),
    has_negative = any(padj < 0.05 & log2FoldChange < -0.5),
    .groups = "drop"
  ) |>
  filter(has_positive & has_negative) |>
  pull(gene_id)

# Three-way divergence
genes_three_way_divergent <- b6_domain_slopes_clean |>
  filter(padj < 0.05, abs(log2FoldChange) > 0.5) |>
  group_by(gene_id) |>
  filter(n() == 3) |>
  pull(gene_id) |>
  unique()

cat("Opposite slopes:", length(genes_opposite_slopes), "genes\n")
cat("Three-way divergent:", length(genes_three_way_divergent), "genes\n\n")

# ===================================================================
# CREATE COMPREHENSIVE SUMMARY TABLE
# ===================================================================

cat("\n===== CREATING SUMMARY TABLES =====\n\n")

pattern_summary_detailed <- tibble(
  Category = c(
    rep("Domain Main Effects", 18),
    rep("Stage Main Effects", 4),
    rep("Diverging", 9),
    rep("Converging", 9),
    rep("Domain-Specific", 12),
    rep("Complex", 2)
  ),
  Pattern = c(
    # Domain main effects
    "MAX higher than PM", "MAX lower than PM",
    "POST higher than PM", "POST lower than PM",
    "POST higher than MAX", "POST lower than MAX",
    "Highest in MAX", "Highest in POST", "Highest in PM",
    "Lowest in MAX", "Lowest in POST", "Lowest in PM",
    "A→P decreasing", "A→P increasing",
    "Any MAX vs PM", "Any POST vs PM", "Any POST vs MAX", "Any domain difference",
    # Stage main effects
    "Increasing", "Decreasing", "Strong increasers", "Strong decreasers",
    # Diverging
    "MAX diverging from PM", "MAX faster than PM", "MAX up PM down",
    "POST diverging from PM", "POST faster than PM", "POST up PM down",
    "POST diverging from MAX", "POST faster than MAX", "POST up MAX down",
    # Converging
    "MAX converging to PM", "MAX faster decrease than PM", "MAX down PM up",
    "POST converging to PM", "POST faster decrease than PM", "POST down PM up",
    "POST converging to MAX", "POST faster decrease than MAX", "POST down MAX up",
    # Domain-specific
    "PM only", "MAX only", "POST only",
    "PM & MAX", "PM & POST", "MAX & POST",
    "PM only ↑", "PM only ↓", "MAX only ↑", "MAX only ↓", 
    "POST only ↑", "POST only ↓",
    # Complex
    "Opposite slopes", "Three-way divergent"
  ),
  N_Genes = c(
    # Domain main effects
    length(genes_MAX_higher_than_PM), length(genes_MAX_lower_than_PM),
    length(genes_POST_higher_than_PM), length(genes_POST_lower_than_PM),
    length(genes_POST_higher_than_MAX), length(genes_POST_lower_than_MAX),
    length(genes_MAX_highest), length(genes_POST_highest), length(genes_PM_highest),
    length(genes_MAX_lowest), length(genes_POST_lowest), length(genes_PM_lowest),
    length(genes_AP_decreasing), length(genes_AP_increasing),
    length(unique(c(genes_MAX_higher_than_PM, genes_MAX_lower_than_PM))),
    length(unique(c(genes_POST_higher_than_PM, genes_POST_lower_than_PM))),
    length(unique(c(genes_POST_higher_than_MAX, genes_POST_lower_than_MAX))),
    length(unique(b6_domain_main_effects_clean$gene_id)),
    # Stage main effects
    length(genes_increasing_with_development), length(genes_decreasing_with_development),
    length(strong_increasers), length(strong_decreasers),
    # Diverging
    length(genes_MAX_diverging_from_PM), length(genes_MAX_increases_faster_than_PM),
    length(genes_MAX_up_PM_down),
    length(genes_POST_diverging_from_PM), length(genes_POST_increases_faster_than_PM),
    length(genes_POST_up_PM_down),
    length(genes_POST_diverging_from_MAX), length(genes_POST_increases_faster_than_MAX),
    length(genes_POST_up_MAX_down),
    # Converging
    length(genes_MAX_converging_to_PM), length(genes_MAX_decreases_faster_than_PM),
    length(genes_MAX_down_PM_up),
    length(genes_POST_converging_to_PM), length(genes_POST_decreases_faster_than_PM),
    length(genes_POST_down_PM_up),
    length(genes_POST_converging_to_MAX), length(genes_POST_decreases_faster_than_MAX),
    length(genes_POST_down_MAX_up),
    # Domain-specific
    length(genes_PM_only_dynamic), length(genes_MAX_only_dynamic),
    length(genes_POST_only_dynamic),
    length(genes_PM_MAX_dynamic), length(genes_PM_POST_dynamic),
    length(genes_MAX_POST_dynamic),
    length(genes_PM_only_increasing), length(genes_PM_only_decreasing),
    length(genes_MAX_only_increasing), length(genes_MAX_only_decreasing),
    length(genes_POST_only_increasing), length(genes_POST_only_decreasing),
    # Complex
    length(genes_opposite_slopes), length(genes_three_way_divergent)
  )
)

print(pattern_summary_detailed)
write_csv(pattern_summary_detailed, "pattern_summary_detailed.csv")

# ===================================================================
# CREATE GENE-LEVEL PATTERN ASSIGNMENTS
# ===================================================================

all_pattern_genes <- unique(c(
  genes_MAX_higher_than_PM, genes_MAX_lower_than_PM,
  genes_POST_higher_than_PM, genes_POST_lower_than_PM,
  genes_POST_higher_than_MAX, genes_POST_lower_than_MAX,
  genes_increasing_with_development, genes_decreasing_with_development,
  genes_MAX_diverging_from_PM, genes_POST_diverging_from_PM, genes_POST_diverging_from_MAX,
  genes_MAX_converging_to_PM, genes_POST_converging_to_PM, genes_POST_converging_to_MAX,
  genes_PM_only_dynamic, genes_MAX_only_dynamic, genes_POST_only_dynamic,
  genes_opposite_slopes, genes_three_way_divergent
))

gene_pattern_assignments <- tibble(gene_id = all_pattern_genes) |>
  left_join(all_genes_palate |> select(gene_id, symbol), by = "gene_id") |>
  mutate(
    # Domain main effects
    MAX_higher_PM = gene_id %in% genes_MAX_higher_than_PM,
    MAX_lower_PM = gene_id %in% genes_MAX_lower_than_PM,
    POST_higher_PM = gene_id %in% genes_POST_higher_than_PM,
    POST_lower_PM = gene_id %in% genes_POST_lower_than_PM,
    POST_higher_MAX = gene_id %in% genes_POST_higher_than_MAX,
    POST_lower_MAX = gene_id %in% genes_POST_lower_than_MAX,
    Highest_MAX = gene_id %in% genes_MAX_highest,
    Highest_POST = gene_id %in% genes_POST_highest,
    Highest_PM = gene_id %in% genes_PM_highest,
    AP_decreasing = gene_id %in% genes_AP_decreasing,
    AP_increasing = gene_id %in% genes_AP_increasing,
    # Stage main effects
    Dev_increasing = gene_id %in% genes_increasing_with_development,
    Dev_decreasing = gene_id %in% genes_decreasing_with_development,
    # Interaction - Diverging
    MAX_div_PM = gene_id %in% genes_MAX_diverging_from_PM,
    POST_div_PM = gene_id %in% genes_POST_diverging_from_PM,
    POST_div_MAX = gene_id %in% genes_POST_diverging_from_MAX,
    # Interaction - Converging
    MAX_conv_PM = gene_id %in% genes_MAX_converging_to_PM,
    POST_conv_PM = gene_id %in% genes_POST_converging_to_PM,
    POST_conv_MAX = gene_id %in% genes_POST_converging_to_MAX,
    # Domain-specific
    PM_only = gene_id %in% genes_PM_only_dynamic,
    MAX_only = gene_id %in% genes_MAX_only_dynamic,
    POST_only = gene_id %in% genes_POST_only_dynamic,
    PM_only_up = gene_id %in% genes_PM_only_increasing,
    MAX_only_up = gene_id %in% genes_MAX_only_increasing,
    POST_only_up = gene_id %in% genes_POST_only_increasing,
    # Complex
    Opposite_slopes = gene_id %in% genes_opposite_slopes,
    Three_way = gene_id %in% genes_three_way_divergent
  ) |>
  mutate(n_patterns = rowSums(across(MAX_higher_PM:Three_way)))

write_csv(gene_pattern_assignments, "gene_pattern_assignments_detailed.csv")
write_csv(domain_specific_detailed, "domain_specific_dynamics_detailed.csv")

cat("\n===== RESULTS SAVED =====\n")
cat("Files:\n")
cat("  1. pattern_summary_detailed.csv\n")
cat("  2. gene_pattern_assignments_detailed.csv\n")
cat("  3. domain_specific_dynamics_detailed.csv\n\n")
cat("Genes classified:", nrow(gene_pattern_assignments), "\n")
cat("Multi-pattern genes:", sum(gene_pattern_assignments$n_patterns > 1), "\n")
cat("\n===== COMPLETE =====\n")
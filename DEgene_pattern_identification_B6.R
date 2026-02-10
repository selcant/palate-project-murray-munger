# ===================================================================
# GENE PATTERN IDENTIFICATION - ENHANCED WITH DIRECTIONALITY
# Classify genes by their developmental trajectories
# Updated version with comprehensive directional patterns
# ===================================================================

library(tidyverse)

cat("===== IDENTIFYING GENE EXPRESSION PATTERNS (ENHANCED) =====\n\n")

# ===================================================================
# PARAMETERS
# ===================================================================
slope_threshold <- 0.5
padj_threshold <- 0.05
rate_diff_threshold <- 0.25  # For distinguishing "faster" vs "similar" rates

# ===================================================================
# PART 1: MAIN EFFECTS - DOMAIN DIFFERENCES (NO INTERACTION)
# These genes show different expression levels across domains
# but maintain parallel developmental trajectories
# ===================================================================

cat("PART 1: DOMAIN MAIN EFFECTS (Parallel Trajectories)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# ----- MAX vs PM Comparisons -----
genes_MAX_higher_than_PM <- b6_domain_main_effects_clean |>
  
  filter(Comparison == "MAX vs PM", log2FoldChange > 0) |>
  pull(gene_id)

genes_MAX_lower_than_PM <- b6_domain_main_effects_clean |>
  filter(Comparison == "MAX vs PM", log2FoldChange < 0) |>
  pull(gene_id)

cat("MAX vs PM:\n")
cat("  MAX higher than PM:", length(genes_MAX_higher_than_PM), "genes\n")
cat("  MAX lower than PM:", length(genes_MAX_lower_than_PM), "genes\n\n")

# ----- POST vs PM Comparisons -----
genes_POST_higher_than_PM <- b6_domain_main_effects_clean |>
  filter(Comparison == "POST vs PM", log2FoldChange > 0) |>
  pull(gene_id)

genes_POST_lower_than_PM <- b6_domain_main_effects_clean |>
  filter(Comparison == "POST vs PM", log2FoldChange < 0) |>
  pull(gene_id)

cat("POST vs PM:\n")
cat("  POST higher than PM:", length(genes_POST_higher_than_PM), "genes\n")
cat("  POST lower than PM:", length(genes_POST_lower_than_PM), "genes\n\n")

# ----- POST vs MAX Comparisons -----
genes_POST_higher_than_MAX <- b6_domain_main_effects_clean |>
  filter(Comparison == "POST vs MAX", log2FoldChange > 0) |>
  pull(gene_id)

genes_POST_lower_than_MAX <- b6_domain_main_effects_clean |>
  filter(Comparison == "POST vs MAX", log2FoldChange < 0) |>
  pull(gene_id)

cat("POST vs MAX:\n")
cat("  POST higher than MAX:", length(genes_POST_higher_than_MAX), "genes\n")
cat("  POST lower than MAX:", length(genes_POST_lower_than_MAX), "genes\n\n")

# ----- Complex Domain Patterns -----
genes_MAX_highest <- intersect(genes_MAX_higher_than_PM, genes_POST_lower_than_MAX)
genes_POST_highest <- intersect(genes_POST_higher_than_PM, genes_POST_higher_than_MAX)
genes_PM_highest <- intersect(genes_MAX_lower_than_PM, genes_POST_lower_than_PM)
genes_MAX_lowest <- intersect(genes_MAX_lower_than_PM, genes_POST_higher_than_MAX)
genes_POST_lowest <- intersect(genes_POST_lower_than_PM, genes_POST_lower_than_MAX)
genes_PM_lowest <- intersect(genes_MAX_higher_than_PM, genes_POST_higher_than_PM)

cat("Complex Domain Patterns:\n")
cat("  Highest in MAX:", length(genes_MAX_highest), "genes\n")
cat("  Highest in POST:", length(genes_POST_highest), "genes\n")
cat("  Highest in PM:", length(genes_PM_highest), "genes\n")
cat("  Lowest in MAX:", length(genes_MAX_lowest), "genes\n")
cat("  Lowest in POST:", length(genes_POST_lowest), "genes\n")
cat("  Lowest in PM:", length(genes_PM_lowest), "genes\n\n")

# ----- Monotonic Gradients -----
genes_AP_decreasing <- intersect(genes_MAX_lower_than_PM, genes_POST_lower_than_MAX)
genes_AP_increasing <- intersect(
  intersect(genes_POST_higher_than_PM, genes_POST_higher_than_MAX),
  genes_MAX_higher_than_PM
)

cat("Monotonic A-P Gradients:\n")
cat("  A→P decreasing (PM > MAX > POST):", length(genes_AP_decreasing), "genes\n")
cat("  A→P increasing (POST > MAX > PM):", length(genes_AP_increasing), "genes\n\n")

# ===================================================================
# PART 2: MAIN EFFECTS - DEVELOPMENTAL CHANGES (NO INTERACTION)
# ===================================================================

cat("\nPART 2: STAGE MAIN EFFECTS (Common Trajectories)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

genes_increasing_with_development <- b6_stage_main_effects_clean |>
  filter(Direction == "Increase") |>
  pull(gene_id)

genes_decreasing_with_development <- b6_stage_main_effects_clean |>
  filter(Direction == "Decrease") |>
  pull(gene_id)

strong_increasers <- b6_stage_main_effects_clean |>
  filter(Direction == "Increase", abs(log2FoldChange) > 1) |>
  pull(gene_id)

strong_decreasers <- b6_stage_main_effects_clean |>
  filter(Direction == "Decrease", abs(log2FoldChange) > 1) |>
  pull(gene_id)

cat("Developmental Changes (all domains together):\n")
cat("  Increasing:", length(genes_increasing_with_development), "genes\n")
cat("  Decreasing:", length(genes_decreasing_with_development), "genes\n")
cat("  Strong increasers (|log2FC| > 1):", length(strong_increasers), "genes\n")
cat("  Strong decreasers (|log2FC| > 1):", length(strong_decreasers), "genes\n\n")

# ===================================================================
# PART 3: COMPREHENSIVE TRAJECTORY PATTERN CLASSIFICATION
# For genes WITH interaction - classify by what slopes are doing
# ===================================================================

cat("\nPART 3: TRAJECTORY PATTERNS (Interaction Genes)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Get slopes in wide format (use UNFILTERED slopes for pattern detection)
slopes_wide <- b6_domain_slopes |>
  filter(gene_id %in% genes_with_interaction) |>
  select(gene_id, symbol, Domain, log2FoldChange, padj) |>
  pivot_wider(
    names_from = Domain,
    values_from = c(log2FoldChange, padj),
    names_glue = "{Domain}_{.value}"
  )

# Classify trajectory patterns comprehensively
trajectory_patterns <- slopes_wide |>
  mutate(
    # === Direction flags ===
    PM_up = PM_log2FoldChange > slope_threshold,
    PM_down = PM_log2FoldChange < -slope_threshold,
    PM_stable = !PM_up & !PM_down,
    PM_sig = PM_padj < padj_threshold & !PM_stable,
    
    MAX_up = MAX_log2FoldChange > slope_threshold,
    MAX_down = MAX_log2FoldChange < -slope_threshold,
    MAX_stable = !MAX_up & !MAX_down,
    MAX_sig = MAX_padj < padj_threshold & !MAX_stable,
    
    POST_up = POST_log2FoldChange > slope_threshold,
    POST_down = POST_log2FoldChange < -slope_threshold,
    POST_stable = !POST_up & !POST_down,
    POST_sig = POST_padj < padj_threshold & !POST_stable,
    
    # === MAX vs PM Trajectory Pattern ===
    MAX_PM_pattern = case_when(
      # Opposite trajectories
      MAX_up & PM_down ~ "Opposite: MAX↑ PM↓",
      MAX_down & PM_up ~ "Opposite: MAX↓ PM↑",
      
      # Co-increasing with rate differences
      MAX_up & PM_up & (MAX_log2FoldChange - PM_log2FoldChange) > rate_diff_threshold ~ 
        "Co-increasing: MAX faster",
      MAX_up & PM_up & (PM_log2FoldChange - MAX_log2FoldChange) > rate_diff_threshold ~ 
        "Co-increasing: PM faster",
      MAX_up & PM_up ~ "Co-increasing: similar",
      
      # Co-decreasing with rate differences
      MAX_down & PM_down & (PM_log2FoldChange - MAX_log2FoldChange) > rate_diff_threshold ~ 
        "Co-decreasing: MAX faster",
      MAX_down & PM_down & (MAX_log2FoldChange - PM_log2FoldChange) > rate_diff_threshold ~ 
        "Co-decreasing: PM faster",
      MAX_down & PM_down ~ "Co-decreasing: similar",
      
      # One dynamic, one stable
      MAX_up & PM_stable ~ "MAX↑ only",
      MAX_down & PM_stable ~ "MAX↓ only",
      PM_up & MAX_stable ~ "PM↑ only",
      PM_down & MAX_stable ~ "PM↓ only",
      
      MAX_stable & PM_stable ~ "Both stable",
      TRUE ~ "Unclassified"
    ),
    
    # === POST vs PM Trajectory Pattern ===
    POST_PM_pattern = case_when(
      POST_up & PM_down ~ "Opposite: POST↑ PM↓",
      POST_down & PM_up ~ "Opposite: POST↓ PM↑",
      
      POST_up & PM_up & (POST_log2FoldChange - PM_log2FoldChange) > rate_diff_threshold ~ 
        "Co-increasing: POST faster",
      POST_up & PM_up & (PM_log2FoldChange - POST_log2FoldChange) > rate_diff_threshold ~ 
        "Co-increasing: PM faster",
      POST_up & PM_up ~ "Co-increasing: similar",
      
      POST_down & PM_down & (PM_log2FoldChange - POST_log2FoldChange) > rate_diff_threshold ~ 
        "Co-decreasing: POST faster",
      POST_down & PM_down & (POST_log2FoldChange - PM_log2FoldChange) > rate_diff_threshold ~ 
        "Co-decreasing: PM faster",
      POST_down & PM_down ~ "Co-decreasing: similar",
      
      POST_up & PM_stable ~ "POST↑ only",
      POST_down & PM_stable ~ "POST↓ only",
      PM_up & POST_stable ~ "PM↑ only",
      PM_down & POST_stable ~ "PM↓ only",
      
      POST_stable & PM_stable ~ "Both stable",
      TRUE ~ "Unclassified"
    ),
    
    # === POST vs MAX Trajectory Pattern ===
    POST_MAX_pattern = case_when(
      POST_up & MAX_down ~ "Opposite: POST↑ MAX↓",
      POST_down & MAX_up ~ "Opposite: POST↓ MAX↑",
      
      POST_up & MAX_up & (POST_log2FoldChange - MAX_log2FoldChange) > rate_diff_threshold ~ 
        "Co-increasing: POST faster",
      POST_up & MAX_up & (MAX_log2FoldChange - POST_log2FoldChange) > rate_diff_threshold ~ 
        "Co-increasing: MAX faster",
      POST_up & MAX_up ~ "Co-increasing: similar",
      
      POST_down & MAX_down & (MAX_log2FoldChange - POST_log2FoldChange) > rate_diff_threshold ~ 
        "Co-decreasing: POST faster",
      POST_down & MAX_down & (POST_log2FoldChange - MAX_log2FoldChange) > rate_diff_threshold ~ 
        "Co-decreasing: MAX faster",
      POST_down & MAX_down ~ "Co-decreasing: similar",
      
      POST_up & MAX_stable ~ "POST↑ only",
      POST_down & MAX_stable ~ "POST↓ only",
      MAX_up & POST_stable ~ "MAX↑ only",
      MAX_down & POST_stable ~ "MAX↓ only",
      
      POST_stable & MAX_stable ~ "Both stable",
      TRUE ~ "Unclassified"
    ),
    
    # === Overall 3-Domain Pattern ===
    overall_pattern = case_when(
      # All same direction
      PM_up & MAX_up & POST_up ~ "All↑",
      PM_down & MAX_down & POST_down ~ "All↓",
      
      # Two up, one down
      PM_up & MAX_up & POST_down ~ "PM+MAX↑, POST↓",
      PM_up & POST_up & MAX_down ~ "PM+POST↑, MAX↓",
      MAX_up & POST_up & PM_down ~ "MAX+POST↑, PM↓",
      
      # Two down, one up
      PM_down & MAX_down & POST_up ~ "PM+MAX↓, POST↑",
      PM_down & POST_down & MAX_up ~ "PM+POST↓, MAX↑",
      MAX_down & POST_down & PM_up ~ "MAX+POST↓, PM↑",
      
      # One dynamic, two stable
      PM_up & MAX_stable & POST_stable ~ "PM↑ only",
      PM_down & MAX_stable & POST_stable ~ "PM↓ only",
      MAX_up & PM_stable & POST_stable ~ "MAX↑ only",
      MAX_down & PM_stable & POST_stable ~ "MAX↓ only",
      POST_up & PM_stable & MAX_stable ~ "POST↑ only",
      POST_down & PM_stable & MAX_stable ~ "POST↓ only",
      
      # Two dynamic same direction, one stable
      PM_up & MAX_up & POST_stable ~ "PM+MAX↑, POST stable",
      PM_up & POST_up & MAX_stable ~ "PM+POST↑, MAX stable",
      MAX_up & POST_up & PM_stable ~ "MAX+POST↑, PM stable",
      PM_down & MAX_down & POST_stable ~ "PM+MAX↓, POST stable",
      PM_down & POST_down & MAX_stable ~ "PM+POST↓, MAX stable",
      MAX_down & POST_down & PM_stable ~ "MAX+POST↓, PM stable",
      
      # Two dynamic opposite, one stable
      PM_up & MAX_down & POST_stable ~ "PM↑ MAX↓, POST stable",
      PM_down & MAX_up & POST_stable ~ "PM↓ MAX↑, POST stable",
      PM_up & POST_down & MAX_stable ~ "PM↑ POST↓, MAX stable",
      PM_down & POST_up & MAX_stable ~ "PM↓ POST↑, MAX stable",
      MAX_up & POST_down & PM_stable ~ "MAX↑ POST↓, PM stable",
      MAX_down & POST_up & PM_stable ~ "MAX↓ POST↑, PM stable",
      
      TRUE ~ "Complex/Other"
    )
  )

# Print summary of trajectory patterns
cat("MAX vs PM Trajectory Patterns:\n")
print(table(trajectory_patterns$MAX_PM_pattern))
cat("\n")

cat("POST vs PM Trajectory Patterns:\n
")
print(table(trajectory_patterns$POST_PM_pattern))
cat("\n")

cat("POST vs MAX Trajectory Patterns:\n")
print(table(trajectory_patterns$POST_MAX_pattern))
cat("\n")

cat("Overall 3-Domain Patterns:\n")
print(table(trajectory_patterns$overall_pattern))
cat("\n")

# ===================================================================
# PART 4: EXTRACT GENE LISTS BY TRAJECTORY PATTERN
# ===================================================================

cat("\nPART 4: EXTRACTING GENE LISTS BY PATTERN\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# --- MAX vs PM patterns ---
genes_MAX_up_PM_down <- trajectory_patterns |>
  filter(MAX_PM_pattern == "Opposite: MAX↑ PM↓") |> pull(gene_id)
genes_MAX_down_PM_up <- trajectory_patterns |>
  filter(MAX_PM_pattern == "Opposite: MAX↓ PM↑") |> pull(gene_id)
genes_coinc_MAX_faster_PM <- trajectory_patterns |>
  filter(MAX_PM_pattern == "Co-increasing: MAX faster") |> pull(gene_id)
genes_coinc_PM_faster_MAX <- trajectory_patterns |>
  filter(MAX_PM_pattern == "Co-increasing: PM faster") |> pull(gene_id)
genes_codec_MAX_faster_PM <- trajectory_patterns |>
  filter(MAX_PM_pattern == "Co-decreasing: MAX faster") |> pull(gene_id)
genes_codec_PM_faster_MAX <- trajectory_patterns |>
  filter(MAX_PM_pattern == "Co-decreasing: PM faster") |> pull(gene_id)

cat("MAX vs PM:\n")
cat("  Opposite (MAX↑ PM↓):", length(genes_MAX_up_PM_down), "\n")
cat("  Opposite (MAX↓ PM↑):", length(genes_MAX_down_PM_up), "\n")
cat("  Co-increasing (MAX faster):", length(genes_coinc_MAX_faster_PM), "\n")
cat("  Co-increasing (PM faster):", length(genes_coinc_PM_faster_MAX), "\n")
cat("  Co-decreasing (MAX faster):", length(genes_codec_MAX_faster_PM), "\n")
cat("  Co-decreasing (PM faster):", length(genes_codec_PM_faster_MAX), "\n\n")

# --- POST vs PM patterns ---
genes_POST_up_PM_down <- trajectory_patterns |>
  filter(POST_PM_pattern == "Opposite: POST↑ PM↓") |> pull(gene_id)
genes_POST_down_PM_up <- trajectory_patterns |>
  filter(POST_PM_pattern == "Opposite: POST↓ PM↑") |> pull(gene_id)
genes_coinc_POST_faster_PM <- trajectory_patterns |>
  filter(POST_PM_pattern == "Co-increasing: POST faster") |> pull(gene_id)
genes_coinc_PM_faster_POST <- trajectory_patterns |>
  filter(POST_PM_pattern == "Co-increasing: PM faster") |> pull(gene_id)
genes_codec_POST_faster_PM <- trajectory_patterns |>
  filter(POST_PM_pattern == "Co-decreasing: POST faster") |> pull(gene_id)
genes_codec_PM_faster_POST <- trajectory_patterns |>
  filter(POST_PM_pattern == "Co-decreasing: PM faster") |> pull(gene_id)

cat("POST vs PM:\n")
cat("  Opposite (POST↑ PM↓):", length(genes_POST_up_PM_down), "\n")
cat("  Opposite (POST↓ PM↑):", length(genes_POST_down_PM_up), "\n")
cat("  Co-increasing (POST faster):", length(genes_coinc_POST_faster_PM), "\n")
cat("  Co-increasing (PM faster):", length(genes_coinc_PM_faster_POST), "\n")
cat("  Co-decreasing (POST faster):", length(genes_codec_POST_faster_PM), "\n")
cat("  Co-decreasing (PM faster):", length(genes_codec_PM_faster_POST), "\n\n")

# --- POST vs MAX patterns ---
genes_POST_up_MAX_down <- trajectory_patterns |>
  filter(POST_MAX_pattern == "Opposite: POST↑ MAX↓") |> pull(gene_id)
genes_POST_down_MAX_up <- trajectory_patterns |>
  filter(POST_MAX_pattern == "Opposite: POST↓ MAX↑") |> pull(gene_id)
genes_coinc_POST_faster_MAX <- trajectory_patterns |>
  filter(POST_MAX_pattern == "Co-increasing: POST faster") |> pull(gene_id)
genes_coinc_MAX_faster_POST <- trajectory_patterns |>
  filter(POST_MAX_pattern == "Co-increasing: MAX faster") |> pull(gene_id)
genes_codec_POST_faster_MAX <- trajectory_patterns |>
  filter(POST_MAX_pattern == "Co-decreasing: POST faster") |> pull(gene_id)
genes_codec_MAX_faster_POST <- trajectory_patterns |>
  filter(POST_MAX_pattern == "Co-decreasing: MAX faster") |> pull(gene_id)

cat("POST vs MAX:\n")
cat("  Opposite (POST↑ MAX↓):", length(genes_POST_up_MAX_down), "\n")
cat("  Opposite (POST↓ MAX↑):", length(genes_POST_down_MAX_up), "\n")
cat("  Co-increasing (POST faster):", length(genes_coinc_POST_faster_MAX), "\n")
cat("  Co-increasing (MAX faster):", length(genes_coinc_MAX_faster_POST), "\n")
cat("  Co-decreasing (POST faster):", length(genes_codec_POST_faster_MAX), "\n")
cat("  Co-decreasing (MAX faster):", length(genes_codec_MAX_faster_POST), "\n\n")

# --- Overall pattern gene lists ---
genes_all_increasing <- trajectory_patterns |>
  filter(overall_pattern == "All↑") |> pull(gene_id)
genes_all_decreasing <- trajectory_patterns |>
  filter(overall_pattern == "All↓") |> pull(gene_id)

cat("Overall patterns:\n")
cat("  All domains increasing:", length(genes_all_increasing), "\n")
cat("  All domains decreasing:", length(genes_all_decreasing), "\n\n")

# ===================================================================
# PART 5: DOMAIN-SPECIFIC DYNAMICS (FROM TRAJECTORY PATTERNS)
# ===================================================================

cat("\nPART 5: DOMAIN-SPECIFIC DYNAMICS\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Single domain dynamics from overall pattern
genes_PM_only_up <- trajectory_patterns |>
  filter(overall_pattern == "PM↑ only") |> pull(gene_id)
genes_PM_only_down <- trajectory_patterns |>
  filter(overall_pattern == "PM↓ only") |> pull(gene_id)
genes_MAX_only_up <- trajectory_patterns |>
  filter(overall_pattern == "MAX↑ only") |> pull(gene_id)
genes_MAX_only_down <- trajectory_patterns |>
  filter(overall_pattern == "MAX↓ only") |> pull(gene_id)
genes_POST_only_up <- trajectory_patterns |>
  filter(overall_pattern == "POST↑ only") |> pull(gene_id)
genes_POST_only_down <- trajectory_patterns |>
  filter(overall_pattern == "POST↓ only") |> pull(gene_id)

cat("Single domain dynamics:\n")
cat("  PM↑ only:", length(genes_PM_only_up), "\n")
cat("  PM↓ only:", length(genes_PM_only_down), "\n")
cat("  MAX↑ only:", length(genes_MAX_only_up), "\n")
cat("  MAX↓ only:", length(genes_MAX_only_down), "\n")
cat("  POST↑ only:", length(genes_POST_only_up), "\n")
cat("  POST↓ only:", length(genes_POST_only_down), "\n\n")

# Two domain dynamics
genes_PM_MAX_up <- trajectory_patterns |>
  filter(overall_pattern == "PM+MAX↑, POST stable") |> pull(gene_id)
genes_PM_MAX_down <- trajectory_patterns |>
  filter(overall_pattern == "PM+MAX↓, POST stable") |> pull(gene_id)
genes_PM_POST_up <- trajectory_patterns |>
  filter(overall_pattern == "PM+POST↑, MAX stable") |> pull(gene_id)
genes_PM_POST_down <- trajectory_patterns |>
  filter(overall_pattern == "PM+POST↓, MAX stable") |> pull(gene_id)
genes_MAX_POST_up <- trajectory_patterns |>
  filter(overall_pattern == "MAX+POST↑, PM stable") |> pull(gene_id)
genes_MAX_POST_down <- trajectory_patterns |>
  filter(overall_pattern == "MAX+POST↓, PM stable") |> pull(gene_id)

cat("Two domain dynamics (same direction):\n")
cat("  PM+MAX↑:", length(genes_PM_MAX_up), "\n")
cat("  PM+MAX↓:", length(genes_PM_MAX_down), "\n")
cat("  PM+POST↑:", length(genes_PM_POST_up), "\n")
cat("  PM+POST↓:", length(genes_PM_POST_down), "\n")
cat("  MAX+POST↑:", length(genes_MAX_POST_up), "\n")
cat("  MAX+POST↓:", length(genes_MAX_POST_down), "\n\n")

# ===================================================================
# PART 6: COMPLEX PATTERNS
# ===================================================================

cat("\nPART 6: COMPLEX PATTERNS\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Genes with opposite slopes (any pair going opposite directions)
genes_opposite_slopes <- trajectory_patterns |>
  filter(
    grepl("Opposite", MAX_PM_pattern) |
      grepl("Opposite", POST_PM_pattern) |
      grepl("Opposite", POST_MAX_pattern)
  ) |>
  pull(gene_id) |>
  unique()

# Three-way: all domains significant and different
genes_three_way_divergent <- trajectory_patterns |>
  filter(PM_sig & MAX_sig & POST_sig) |>
  filter(
    !(PM_up & MAX_up & POST_up) &
      !(PM_down & MAX_down & POST_down)
  ) |>
  pull(gene_id)

cat("Complex patterns:\n")
cat("  Any opposite slopes:", length(genes_opposite_slopes), "\n")
cat("  Three-way divergent:", length(genes_three_way_divergent), "\n\n")

# ===================================================================
# CREATE COMPREHENSIVE SUMMARY TABLES
# ===================================================================

cat("\n===== CREATING SUMMARY TABLES =====\n\n")

# Pattern summary with new categories
pattern_summary_detailed <- tibble(
  Category = c(
    # Domain Main Effects
    rep("Domain Main Effects", 14),
    # Stage Main Effects
    rep("Stage Main Effects", 4),
    # Trajectory Patterns - MAX vs PM
    rep("Trajectories: MAX vs PM", 6),
    # Trajectory Patterns - POST vs PM
    rep("Trajectories: POST vs PM", 6),
    # Trajectory Patterns - POST vs MAX
    rep("Trajectories: POST vs MAX", 6),
    # Domain-Specific Dynamics
    rep("Domain-Specific", 12),
    # Complex
    rep("Complex", 2)
  ),
  Pattern = c(
    # Domain main effects
    "MAX higher than PM", "MAX lower than PM",
    "POST higher than PM", "POST lower than PM",
    "POST higher than MAX", "POST lower than MAX",
    "Highest in MAX", "Highest in POST", "Highest in PM",
    "Lowest in MAX", "Lowest in POST",
    "A→P decreasing", "A→P increasing",
    "Any domain difference",
    # Stage main effects
    "All↑ (common)", "All↓ (common)", 
    "Strong increasers", "Strong decreasers",
    # MAX vs PM trajectories
    "Opposite: MAX↑ PM↓", "Opposite: MAX↓ PM↑",
    "Co-increasing: MAX faster", "Co-increasing: PM faster",
    "Co-decreasing: MAX faster", "Co-decreasing: PM faster",
    # POST vs PM trajectories
    "Opposite: POST↑ PM↓", "Opposite: POST↓ PM↑",
    "Co-increasing: POST faster", "Co-increasing: PM faster",
    "Co-decreasing: POST faster", "Co-decreasing: PM faster",
    # POST vs MAX trajectories
    "Opposite: POST↑ MAX↓", "Opposite: POST↓ MAX↑",
    "Co-increasing: POST faster", "Co-increasing: MAX faster",
    "Co-decreasing: POST faster", "Co-decreasing: MAX faster",
    # Domain-specific
    "PM↑ only", "PM↓ only", "MAX↑ only", "MAX↓ only", "POST↑ only", "POST↓ only",
    "PM+MAX↑", "PM+MAX↓", "PM+POST↑", "PM+POST↓", "MAX+POST↑", "MAX+POST↓",
    # Complex
    "Any opposite slopes", "Three-way divergent"
  ),
  N_Genes = c(
    # Domain main effects
    length(genes_MAX_higher_than_PM), length(genes_MAX_lower_than_PM),
    length(genes_POST_higher_than_PM), length(genes_POST_lower_than_PM),
    length(genes_POST_higher_than_MAX), length(genes_POST_lower_than_MAX),
    length(genes_MAX_highest), length(genes_POST_highest), length(genes_PM_highest),
    length(genes_MAX_lowest), length(genes_POST_lowest),
    length(genes_AP_decreasing), length(genes_AP_increasing),
    length(unique(b6_domain_main_effects_clean$gene_id)),
    # Stage main effects
    length(genes_increasing_with_development), length(genes_decreasing_with_development),
    length(strong_increasers), length(strong_decreasers),
    # MAX vs PM trajectories
    length(genes_MAX_up_PM_down), length(genes_MAX_down_PM_up),
    length(genes_coinc_MAX_faster_PM), length(genes_coinc_PM_faster_MAX),
    length(genes_codec_MAX_faster_PM), length(genes_codec_PM_faster_MAX),
    # POST vs PM trajectories
    length(genes_POST_up_PM_down), length(genes_POST_down_PM_up),
    length(genes_coinc_POST_faster_PM), length(genes_coinc_PM_faster_POST),
    length(genes_codec_POST_faster_PM), length(genes_codec_PM_faster_POST),
    # POST vs MAX trajectories
    length(genes_POST_up_MAX_down), length(genes_POST_down_MAX_up),
    length(genes_coinc_POST_faster_MAX), length(genes_coinc_MAX_faster_POST),
    length(genes_codec_POST_faster_MAX), length(genes_codec_MAX_faster_POST),
    # Domain-specific
    length(genes_PM_only_up), length(genes_PM_only_down),
    length(genes_MAX_only_up), length(genes_MAX_only_down),
    length(genes_POST_only_up), length(genes_POST_only_down),
    length(genes_PM_MAX_up), length(genes_PM_MAX_down),
    length(genes_PM_POST_up), length(genes_PM_POST_down),
    length(genes_MAX_POST_up), length(genes_MAX_POST_down),
    # Complex
    length(genes_opposite_slopes), length(genes_three_way_divergent)
  )
)

print(pattern_summary_detailed)

# ===================================================================
# CREATE GENE-LEVEL PATTERN ASSIGNMENTS
# ===================================================================

# Get all genes
all_pattern_genes <- unique(c(
  genes_MAX_higher_than_PM, genes_MAX_lower_than_PM,
  genes_POST_higher_than_PM, genes_POST_lower_than_PM,
  genes_POST_higher_than_MAX, genes_POST_lower_than_MAX,
  genes_increasing_with_development, genes_decreasing_with_development,
  trajectory_patterns$gene_id
))

gene_pattern_assignments <- tibble(gene_id = all_pattern_genes) |>
  left_join(all_genes_palate |> select(gene_id, symbol), by = "gene_id") |>
  # Domain main effects
  mutate(
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
    # MAX vs PM trajectories
    MAX_up_PM_down = gene_id %in% genes_MAX_up_PM_down,
    MAX_down_PM_up = gene_id %in% genes_MAX_down_PM_up,
    Coinc_MAX_faster_PM = gene_id %in% genes_coinc_MAX_faster_PM,
    Coinc_PM_faster_MAX = gene_id %in% genes_coinc_PM_faster_MAX,
    Codec_MAX_faster_PM = gene_id %in% genes_codec_MAX_faster_PM,
    Codec_PM_faster_MAX = gene_id %in% genes_codec_PM_faster_MAX,
    # POST vs PM trajectories
    POST_up_PM_down = gene_id %in% genes_POST_up_PM_down,
    POST_down_PM_up = gene_id %in% genes_POST_down_PM_up,
    Coinc_POST_faster_PM = gene_id %in% genes_coinc_POST_faster_PM,
    Coinc_PM_faster_POST = gene_id %in% genes_coinc_PM_faster_POST,
    Codec_POST_faster_PM = gene_id %in% genes_codec_POST_faster_PM,
    Codec_PM_faster_POST = gene_id %in% genes_codec_PM_faster_POST,
    # POST vs MAX trajectories
    POST_up_MAX_down = gene_id %in% genes_POST_up_MAX_down,
    POST_down_MAX_up = gene_id %in% genes_POST_down_MAX_up,
    Coinc_POST_faster_MAX = gene_id %in% genes_coinc_POST_faster_MAX,
    Coinc_MAX_faster_POST = gene_id %in% genes_coinc_MAX_faster_POST,
    Codec_POST_faster_MAX = gene_id %in% genes_codec_POST_faster_MAX,
    Codec_MAX_faster_POST = gene_id %in% genes_codec_MAX_faster_POST,
    # Domain-specific
    PM_only_up = gene_id %in% genes_PM_only_up,
    PM_only_down = gene_id %in% genes_PM_only_down,
    MAX_only_up = gene_id %in% genes_MAX_only_up,
    MAX_only_down = gene_id %in% genes_MAX_only_down,
    POST_only_up = gene_id %in% genes_POST_only_up,
    POST_only_down = gene_id %in% genes_POST_only_down,
    # Complex
    Opposite_slopes = gene_id %in% genes_opposite_slopes,
    Three_way = gene_id %in% genes_three_way_divergent
  ) |>
  # Add trajectory pattern labels
  left_join(
    trajectory_patterns |> 
      select(gene_id, MAX_PM_pattern, POST_PM_pattern, POST_MAX_pattern, overall_pattern),
    by = "gene_id"
  ) |>
  mutate(n_patterns = rowSums(across(MAX_higher_PM:Three_way)))

# ===================================================================
# SAVE OUTPUT FILES
# ===================================================================

write_csv(pattern_summary_detailed, "pattern_summary_detailed.csv")
write_csv(gene_pattern_assignments, "gene_pattern_assignments_detailed.csv")
write_csv(trajectory_patterns, "trajectory_patterns_detailed.csv")

cat("\n===== RESULTS SAVED =====\n")
cat("Files:\n")
cat("  1. pattern_summary_detailed.csv\n")
cat("  2. gene_pattern_assignments_detailed.csv\n")
cat("  3. trajectory_patterns_detailed.csv\n\n")
cat("Genes classified:", nrow(gene_pattern_assignments), "\n")
cat("Interaction genes with trajectories:", nrow(trajectory_patterns), "\n")
cat("\n===== COMPLETE =====\n")
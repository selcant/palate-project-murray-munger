# ===================================================================
# STRAIN GENE PATTERN IDENTIFICATION
# Classify genes by their strain-specific expression patterns
# ===================================================================

library(tidyverse)

cat("===== IDENTIFYING STRAIN EXPRESSION PATTERNS =====\n\n")

# ===================================================================
# HELPER FUNCTIONS
# ===================================================================

# Function to get normalized counts for plotting
get_normalized_counts <- function(dds, genes) {
  vsd <- vst(dds, blind = FALSE)
  counts_norm <- assay(vsd)[genes, , drop = FALSE]
  
  counts_df <- counts_norm |>
    as.data.frame() |>
    rownames_to_column("gene_id") |>
    pivot_longer(-gene_id, names_to = "sample", values_to = "norm_count") |>
    left_join(colData(dds) |> 
                as.data.frame() |> 
                rownames_to_column("sample"))
  
  return(counts_df)
}

# Function to calculate average expression by condition
calc_avg_by_condition <- function(counts_df, group_vars) {
  counts_df |>
    group_by(gene_id, !!!syms(group_vars)) |>
    summarize(
      mean_expr = mean(norm_count),
      se_expr = sd(norm_count) / sqrt(n()),
      .groups = "drop"
    )
}

# ===================================================================
# PART 1: MAIN STRAIN EFFECT PATTERNS
# These genes show consistent CAST vs B6 differences
# ===================================================================

cat("\nPART 1: MAIN STRAIN EFFECT PATTERNS\n")
cat("=" |> rep(60) |> paste(collapse=""), "\n\n")

genes_main_strain <- genes_main_only$gene_id

if(length(genes_main_strain) > 0) {
  # Get average expression by strain and domain
  main_counts <- get_normalized_counts(strain_dds_lrt, genes_main_strain)
  
  main_avg <- main_counts |>
    group_by(gene_id, strain, AP_domain) |>
    summarize(mean_expr = mean(norm_count), .groups = "drop") |>
    pivot_wider(names_from = c(strain, AP_domain), 
                values_from = mean_expr,
                names_sep = "_")
  
  # Calculate overall strain differences
  main_patterns <- genes_main_only |>
    select(gene_id, symbol, main_effect_lfc) |>
    mutate(
      strain_pattern = case_when(
        main_effect_lfc > 1 ~ "CAST >> B6 (strong)",
        main_effect_lfc > 0.5 ~ "CAST > B6 (moderate)",
        main_effect_lfc > 0 ~ "CAST > B6 (weak)",
        main_effect_lfc < -1 ~ "B6 >> CAST (strong)",
        main_effect_lfc < -0.5 ~ "B6 > CAST (moderate)",
        main_effect_lfc < 0 ~ "B6 > CAST (weak)",
        TRUE ~ "No difference"
      )
    )
  
  # Summary
  pattern_counts <- main_patterns |> count(strain_pattern)
  print(pattern_counts)
  
  # Extract gene lists
  genes_CAST_high_strong <- main_patterns |>
    filter(strain_pattern == "CAST >> B6 (strong)") |> pull(gene_id)
  
  genes_CAST_high_mod <- main_patterns |>
    filter(strain_pattern == "CAST > B6 (moderate)") |> pull(gene_id)
  
  genes_B6_high_strong <- main_patterns |>
    filter(strain_pattern == "B6 >> CAST (strong)") |> pull(gene_id)
  
  genes_B6_high_mod <- main_patterns |>
    filter(strain_pattern == "B6 > CAST (moderate)") |> pull(gene_id)
  
  cat("\nMain effect genes:\n")
  cat("  CAST >> B6 (strong):", length(genes_CAST_high_strong), "\n")
  cat("  CAST > B6 (moderate):", length(genes_CAST_high_mod), "\n")
  cat("  B6 >> CAST (strong):", length(genes_B6_high_strong), "\n")
  cat("  B6 > CAST (moderate):", length(genes_B6_high_mod), "\n\n")
}

# ===================================================================
# PART 2: DOMAIN × STRAIN INTERACTION PATTERNS
# Which domain shows the biggest strain difference?
# ===================================================================

cat("\nPART 2: DOMAIN × STRAIN INTERACTION PATTERNS\n")
cat("=" |> rep(60) |> paste(collapse=""), "\n\n")

genes_domain_strain <- genes_domain_int$gene_id

if(length(genes_domain_strain) > 0) {
  # For each gene, calculate strain differences in each domain
  domain_counts <- get_normalized_counts(strain_dds_lrt, genes_domain_strain)
  
  domain_strain_diffs <- domain_counts |>
    group_by(gene_id, strain, AP_domain) |>
    summarize(mean_expr = mean(norm_count), .groups = "drop") |>
    pivot_wider(names_from = strain, values_from = mean_expr) |>
    mutate(strain_diff = CAST - B6) |>
    select(gene_id, AP_domain, strain_diff) |>
    pivot_wider(names_from = AP_domain, 
                values_from = strain_diff,
                names_prefix = "diff_")
  
  # Classify patterns
  domain_patterns <- domain_strain_diffs |>
    left_join(genes_domain_int |> select(gene_id, symbol)) |>
    mutate(
      # Which domain has biggest absolute difference?
      max_diff_domain = case_when(
        abs(diff_PM) > abs(diff_MAX) & abs(diff_PM) > abs(diff_POST) ~ "PM",
        abs(diff_MAX) > abs(diff_PM) & abs(diff_MAX) > abs(diff_POST) ~ "MAX",
        abs(diff_POST) > abs(diff_PM) & abs(diff_POST) > abs(diff_MAX) ~ "POST",
        TRUE ~ "Similar"
      ),
      
      # Direction of biggest difference
      max_diff_direction = case_when(
        max_diff_domain == "PM" & diff_PM > 0 ~ "PM: CAST > B6",
        max_diff_domain == "PM" & diff_PM < 0 ~ "PM: B6 > CAST",
        max_diff_domain == "MAX" & diff_MAX > 0 ~ "MAX: CAST > B6",
        max_diff_domain == "MAX" & diff_MAX < 0 ~ "MAX: B6 > CAST",
        max_diff_domain == "POST" & diff_POST > 0 ~ "POST: CAST > B6",
        max_diff_domain == "POST" & diff_POST < 0 ~ "POST: B6 > CAST",
        TRUE ~ "Similar across domains"
      ),
      
      # Pattern: are differences in same direction?
      consistency = case_when(
        sign(diff_PM) == sign(diff_MAX) & 
          sign(diff_MAX) == sign(diff_POST) ~ "Consistent direction",
        TRUE ~ "Opposite directions"
      ),
      
      # Detailed pattern
      detailed_pattern = paste(max_diff_direction, consistency, sep = " | ")
    )
  
  # Summary
  pattern_counts <- domain_patterns |> count(detailed_pattern) |>
    arrange(desc(n))
  print(pattern_counts)
  
  # Extract key gene lists
  genes_PM_biggest_diff <- domain_patterns |>
    filter(max_diff_domain == "PM") |> pull(gene_id)
  
  genes_MAX_biggest_diff <- domain_patterns |>
    filter(max_diff_domain == "MAX") |> pull(gene_id)
  
  genes_POST_biggest_diff <- domain_patterns |>
    filter(max_diff_domain == "POST") |> pull(gene_id)
  
  genes_PM_CAST_high <- domain_patterns |>
    filter(max_diff_domain == "PM", diff_PM > 0) |> pull(gene_id)
  
  genes_MAX_CAST_high <- domain_patterns |>
    filter(max_diff_domain == "MAX", diff_MAX > 0) |> pull(gene_id)
  
  genes_POST_CAST_high <- domain_patterns |>
    filter(max_diff_domain == "POST", diff_POST > 0) |> pull(gene_id)
  
  genes_opposite_directions <- domain_patterns |>
    filter(consistency == "Opposite directions") |> pull(gene_id)
  
  cat("\nDomain × Strain patterns:\n")
  cat("  Biggest difference in PM:", length(genes_PM_biggest_diff), "\n")
  cat("  Biggest difference in MAX:", length(genes_MAX_biggest_diff), "\n")
  cat("  Biggest difference in POST:", length(genes_POST_biggest_diff), "\n")
  cat("  PM shows CAST > B6:", length(genes_PM_CAST_high), "\n")
  cat("  MAX shows CAST > B6:", length(genes_MAX_CAST_high), "\n")
  cat("  POST shows CAST > B6:", length(genes_POST_CAST_high), "\n")
  cat("  Opposite directions across domains:", length(genes_opposite_directions), "\n\n")
}

# ===================================================================
# PART 3: TIME × STRAIN INTERACTION PATTERNS
# Do strains diverge or converge over development?
# ===================================================================

cat("\nPART 3: TIME × STRAIN INTERACTION PATTERNS\n")
cat("=" |> rep(60) |> paste(collapse=""), "\n\n")

genes_time_strain <- genes_time_int$gene_id

if(length(genes_time_strain) > 0) {
  # Calculate slopes for each strain
  time_counts <- get_normalized_counts(strain_dds_lrt, genes_time_strain)
  
  # Fit linear models to get slopes for each gene and strain
  time_slopes <- time_counts |>
    group_by(gene_id, strain) |>
    do({
      model <- lm(norm_count ~ LB_stage_from_start, data = .)
      tibble(
        slope = coef(model)[2],
        intercept = coef(model)[1],
        r_squared = summary(model)$r.squared
      )
    }) |>
    ungroup()
  
  # Compare slopes between strains
  time_patterns <- time_slopes |>
    pivot_wider(names_from = strain, 
                values_from = c(slope, intercept, r_squared),
                names_sep = "_") |>
    mutate(
      slope_diff = slope_CAST - slope_B6,
      
      # Classify trajectory pattern
      trajectory_pattern = case_when(
        # Diverging (slopes have same sign, CAST steeper)
        slope_B6 > 0 & slope_CAST > 0 & slope_diff > 0.1 ~ "Both increase, CAST faster",
        slope_B6 < 0 & slope_CAST < 0 & slope_diff < -0.1 ~ "Both decrease, CAST faster",
        slope_B6 > 0 & slope_CAST > 0 & slope_diff < -0.1 ~ "Both increase, B6 faster",
        slope_B6 < 0 & slope_CAST < 0 & slope_diff > 0.1 ~ "Both decrease, B6 faster",
        
        # Opposite directions
        slope_B6 > 0 & slope_CAST < 0 ~ "B6 increases, CAST decreases",
        slope_B6 < 0 & slope_CAST > 0 ~ "CAST increases, B6 decreases",
        
        # Flat vs dynamic
        abs(slope_B6) < 0.1 & abs(slope_CAST) > 0.1 ~ "CAST dynamic, B6 flat",
        abs(slope_CAST) < 0.1 & abs(slope_B6) > 0.1 ~ "B6 dynamic, CAST flat",
        
        TRUE ~ "Parallel"
      ),
      
      # Diverging vs converging
      divergence_pattern = case_when(
        abs(slope_diff) < 0.1 ~ "Parallel",
        slope_diff > 0.1 ~ "Strains diverging (CAST steeper)",
        slope_diff < -0.1 ~ "Strains diverging (B6 steeper)",
        TRUE ~ "Similar"
      )
    ) |>
    left_join(genes_time_int |> select(gene_id, symbol))
  
  # Summary
  pattern_counts <- time_patterns |> count(trajectory_pattern) |>
    arrange(desc(n))
  print(pattern_counts)
  
  # Extract gene lists
  genes_both_increase_CAST_faster <- time_patterns |>
    filter(trajectory_pattern == "Both increase, CAST faster") |> pull(gene_id)
  
  genes_both_decrease_CAST_faster <- time_patterns |>
    filter(trajectory_pattern == "Both decrease, CAST faster") |> pull(gene_id)
  
  genes_both_increase_B6_faster <- time_patterns |>
    filter(trajectory_pattern == "Both increase, B6 faster") |> pull(gene_id)
  
  genes_both_decrease_B6_faster <- time_patterns |>
    filter(trajectory_pattern == "Both decrease, B6 faster") |> pull(gene_id)
  
  genes_B6_up_CAST_down <- time_patterns |>
    filter(trajectory_pattern == "B6 increases, CAST decreases") |> pull(gene_id)
  
  genes_CAST_up_B6_down <- time_patterns |>
    filter(trajectory_pattern == "CAST increases, B6 decreases") |> pull(gene_id)
  
  genes_CAST_dynamic_B6_flat <- time_patterns |>
    filter(trajectory_pattern == "CAST dynamic, B6 flat") |> pull(gene_id)
  
  genes_B6_dynamic_CAST_flat <- time_patterns |>
    filter(trajectory_pattern == "B6 dynamic, CAST flat") |> pull(gene_id)
  
  cat("\nTime × Strain patterns:\n")
  cat("  Both ↑, CAST faster:", length(genes_both_increase_CAST_faster), "\n")
  cat("  Both ↓, CAST faster:", length(genes_both_decrease_CAST_faster), "\n")
  cat("  Both ↑, B6 faster:", length(genes_both_increase_B6_faster), "\n")
  cat("  Both ↓, B6 faster:", length(genes_both_decrease_B6_faster), "\n")
  cat("  B6 ↑, CAST ↓:", length(genes_B6_up_CAST_down), "\n")
  cat("  CAST ↑, B6 ↓:", length(genes_CAST_up_B6_down), "\n")
  cat("  CAST dynamic, B6 flat:", length(genes_CAST_dynamic_B6_flat), "\n")
  cat("  B6 dynamic, CAST flat:", length(genes_B6_dynamic_CAST_flat), "\n\n")
}

# ===================================================================
# PART 4: THREE-WAY INTERACTION PATTERNS
# Domain-specific strain differences that change over time
# ===================================================================

cat("\nPART 4: THREE-WAY INTERACTION PATTERNS\n")
cat("=" |> rep(60) |> paste(collapse=""), "\n\n")

genes_three_way <- genes_three_ways$gene_id

if(length(genes_three_way) > 0) {
  # For each gene, fit slopes in each domain-strain combination
  threeway_counts <- get_normalized_counts(strain_dds_lrt, genes_three_way)
  
  threeway_slopes <- threeway_counts |>
    group_by(gene_id, strain, AP_domain) |>
    do({
      if(nrow(.) > 3) {  # Changed from n() to nrow(.)
        model <- lm(norm_count ~ LB_stage_from_start, data = .)
        tibble(
          slope = coef(model)[2],
          r_squared = summary(model)$r.squared
        )
      } else {
        tibble(slope = NA, r_squared = NA)
      }
    }) |>
    ungroup()
  
  # Reshape to compare slopes
  threeway_wide <- threeway_slopes |>
    unite("condition", strain, AP_domain, sep = "_") |>
    select(gene_id, condition, slope) |>
    pivot_wider(names_from = condition, values_from = slope, names_prefix = "slope_")
  
  # Calculate strain differences in slopes for each domain
  threeway_patterns <- threeway_wide |>
    mutate(
      slope_diff_PM = slope_CAST_PM - slope_B6_PM,
      slope_diff_MAX = slope_CAST_MAX - slope_B6_MAX,
      slope_diff_POST = slope_CAST_POST - slope_B6_POST,
      
      # Which domain shows biggest slope difference?
      max_slope_diff_domain = case_when(
        abs(slope_diff_PM) > abs(slope_diff_MAX) & 
          abs(slope_diff_PM) > abs(slope_diff_POST) ~ "PM",
        abs(slope_diff_MAX) > abs(slope_diff_PM) & 
          abs(slope_diff_MAX) > abs(slope_diff_POST) ~ "MAX",
        abs(slope_diff_POST) > abs(slope_diff_PM) & 
          abs(slope_diff_POST) > abs(slope_diff_MAX) ~ "POST",
        TRUE ~ "Similar"
      ),
      
      # Are slope differences consistent across domains?
      consistency = case_when(
        !is.na(slope_diff_PM) & !is.na(slope_diff_MAX) & !is.na(slope_diff_POST) &
          sign(slope_diff_PM) == sign(slope_diff_MAX) & 
          sign(slope_diff_MAX) == sign(slope_diff_POST) ~ "Consistent",
        TRUE ~ "Variable"
      ),
      
      # Detailed pattern
      pattern = case_when(
        max_slope_diff_domain == "PM" & slope_diff_PM > 0.1 ~ 
          "PM: Strains diverge most (CAST steeper)",
        max_slope_diff_domain == "PM" & slope_diff_PM < -0.1 ~ 
          "PM: Strains diverge most (B6 steeper)",
        max_slope_diff_domain == "MAX" & slope_diff_MAX > 0.1 ~ 
          "MAX: Strains diverge most (CAST steeper)",
        max_slope_diff_domain == "MAX" & slope_diff_MAX < -0.1 ~ 
          "MAX: Strains diverge most (B6 steeper)",
        max_slope_diff_domain == "POST" & slope_diff_POST > 0.1 ~ 
          "POST: Strains diverge most (CAST steeper)",
        max_slope_diff_domain == "POST" & slope_diff_POST < -0.1 ~ 
          "POST: Strains diverge most (B6 steeper)",
        TRUE ~ "Complex"
      )
    ) |>
    left_join(genes_three_ways |> select(gene_id, symbol))
  
  # Summary
  pattern_counts <- threeway_patterns |> count(pattern) |>
    arrange(desc(n))
  print(pattern_counts)
  
  # Extract gene lists
  genes_PM_biggest_slope_diff <- threeway_patterns |>
    filter(max_slope_diff_domain == "PM") |> pull(gene_id)
  
  genes_MAX_biggest_slope_diff <- threeway_patterns |>
    filter(max_slope_diff_domain == "MAX") |> pull(gene_id)
  
  genes_POST_biggest_slope_diff <- threeway_patterns |>
    filter(max_slope_diff_domain == "POST") |> pull(gene_id)
  
  genes_consistent_slope_diffs <- threeway_patterns |>
    filter(consistency == "Consistent") |> pull(gene_id)
  
  cat("\n3-way interaction patterns:\n")
  cat("  Biggest slope diff in PM:", length(genes_PM_biggest_slope_diff), "\n")
  cat("  Biggest slope diff in MAX:", length(genes_MAX_biggest_slope_diff), "\n")
  cat("  Biggest slope diff in POST:", length(genes_POST_biggest_slope_diff), "\n")
  cat("  Consistent slope differences:", length(genes_consistent_slope_diffs), "\n\n")
}

# ===================================================================
# PART 5: DOMAIN-SPECIFIC STRAIN EFFECTS
# Are strain differences specific to certain domains?
# ===================================================================

cat("\nPART 5: DOMAIN-SPECIFIC STRAIN EFFECTS\n")
cat("=" |> rep(60) |> paste(collapse=""), "\n\n")

# For genes with domain×strain interaction, identify which domains show effects
if(length(genes_domain_strain) > 0) {
  
  # Calculate strain difference in each domain for each gene
  domain_specific <- domain_counts |>
    group_by(gene_id, strain, AP_domain) |>
    summarize(mean_expr = mean(norm_count), .groups = "drop") |>
    pivot_wider(names_from = strain, values_from = mean_expr) |>
    mutate(
      strain_diff = CAST - B6,
      abs_diff = abs(strain_diff),
      has_effect = abs_diff > 0.5  # threshold for meaningful difference
    ) |>
    group_by(gene_id) |>
    mutate(
      max_abs_diff = max(abs_diff),
      is_max = abs_diff == max_abs_diff
    ) |>
    ungroup()
  
  # Classify by which domains show effects
  domain_specificity <- domain_specific |>
    select(gene_id, AP_domain, has_effect) |>
    pivot_wider(names_from = AP_domain, 
                values_from = has_effect,
                names_prefix = "effect_") |>
    mutate(
      specificity_pattern = case_when(
        effect_PM & !effect_MAX & !effect_POST ~ "PM only",
        !effect_PM & effect_MAX & !effect_POST ~ "MAX only",
        !effect_PM & !effect_MAX & effect_POST ~ "POST only",
        effect_PM & effect_MAX & !effect_POST ~ "PM & MAX",
        effect_PM & !effect_MAX & effect_POST ~ "PM & POST",
        !effect_PM & effect_MAX & effect_POST ~ "MAX & POST",
        effect_PM & effect_MAX & effect_POST ~ "All domains",
        TRUE ~ "None"
      )
    ) |>
    left_join(genes_domain_int |> select(gene_id, symbol))
  
  # Summary
  pattern_counts <- domain_specificity |> count(specificity_pattern) |>
    arrange(desc(n))
  print(pattern_counts)
  
  # Extract gene lists
  genes_strain_PM_only <- domain_specificity |>
    filter(specificity_pattern == "PM only") |> pull(gene_id)
  
  genes_strain_MAX_only <- domain_specificity |>
    filter(specificity_pattern == "MAX only") |> pull(gene_id)
  
  genes_strain_POST_only <- domain_specificity |>
    filter(specificity_pattern == "POST only") |> pull(gene_id)
  
  genes_strain_all_domains <- domain_specificity |>
    filter(specificity_pattern == "All domains") |> pull(gene_id)
  
  cat("\nDomain-specific strain effects:\n")
  cat("  PM only:", length(genes_strain_PM_only), "\n")
  cat("  MAX only:", length(genes_strain_MAX_only), "\n")
  cat("  POST only:", length(genes_strain_POST_only), "\n")
  cat("  All domains:", length(genes_strain_all_domains), "\n\n")
}

# ===================================================================
# PART 6: CREATE COMPREHENSIVE SUMMARY
# ===================================================================

cat("\n===== CREATING SUMMARY TABLES =====\n\n")

# Create pattern summary
strain_pattern_summary <- tibble(
  Category = c(
    rep("Main Strain Effect", 4),
    rep("Domain × Strain", 7),
    rep("Time × Strain", 8),
    rep("3-Way Interaction", 4),
    rep("Domain Specificity", 4)
  ),
  Pattern = c(
    # Main strain
    "CAST >> B6 (strong)", "CAST > B6 (moderate)",
    "B6 >> CAST (strong)", "B6 > CAST (moderate)",
    # Domain × Strain
    "Biggest diff in PM", "Biggest diff in MAX", "Biggest diff in POST",
    "PM: CAST > B6", "MAX: CAST > B6", "POST: CAST > B6",
    "Opposite directions across domains",
    # Time × Strain
    "Both ↑, CAST faster", "Both ↓, CAST faster",
    "Both ↑, B6 faster", "Both ↓, B6 faster",
    "B6 ↑, CAST ↓", "CAST ↑, B6 ↓",
    "CAST dynamic, B6 flat", "B6 dynamic, CAST flat",
    # 3-way
    "Biggest slope diff in PM", "Biggest slope diff in MAX",
    "Biggest slope diff in POST", "Consistent slope differences",
    # Domain specificity
    "PM only", "MAX only", "POST only", "All domains"
  ),
  N_Genes = c(
    # Main strain
    if(exists("genes_CAST_high_strong")) length(genes_CAST_high_strong) else 0,
    if(exists("genes_CAST_high_mod")) length(genes_CAST_high_mod) else 0,
    if(exists("genes_B6_high_strong")) length(genes_B6_high_strong) else 0,
    if(exists("genes_B6_high_mod")) length(genes_B6_high_mod) else 0,
    # Domain × Strain
    if(exists("genes_PM_biggest_diff")) length(genes_PM_biggest_diff) else 0,
    if(exists("genes_MAX_biggest_diff")) length(genes_MAX_biggest_diff) else 0,
    if(exists("genes_POST_biggest_diff")) length(genes_POST_biggest_diff) else 0,
    if(exists("genes_PM_CAST_high")) length(genes_PM_CAST_high) else 0,
    if(exists("genes_MAX_CAST_high")) length(genes_MAX_CAST_high) else 0,
    if(exists("genes_POST_CAST_high")) length(genes_POST_CAST_high) else 0,
    if(exists("genes_opposite_directions")) length(genes_opposite_directions) else 0,
    # Time × Strain
    if(exists("genes_both_increase_CAST_faster")) length(genes_both_increase_CAST_faster) else 0,
    if(exists("genes_both_decrease_CAST_faster")) length(genes_both_decrease_CAST_faster) else 0,
    if(exists("genes_both_increase_B6_faster")) length(genes_both_increase_B6_faster) else 0,
    if(exists("genes_both_decrease_B6_faster")) length(genes_both_decrease_B6_faster) else 0,
    if(exists("genes_B6_up_CAST_down")) length(genes_B6_up_CAST_down) else 0,
    if(exists("genes_CAST_up_B6_down")) length(genes_CAST_up_B6_down) else 0,
    if(exists("genes_CAST_dynamic_B6_flat")) length(genes_CAST_dynamic_B6_flat) else 0,
    if(exists("genes_B6_dynamic_CAST_flat")) length(genes_B6_dynamic_CAST_flat) else 0,
    # 3-way
    if(exists("genes_PM_biggest_slope_diff")) length(genes_PM_biggest_slope_diff) else 0,
    if(exists("genes_MAX_biggest_slope_diff")) length(genes_MAX_biggest_slope_diff) else 0,
    if(exists("genes_POST_biggest_slope_diff")) length(genes_POST_biggest_slope_diff) else 0,
    if(exists("genes_consistent_slope_diffs")) length(genes_consistent_slope_diffs) else 0,
    # Domain specificity
    if(exists("genes_strain_PM_only")) length(genes_strain_PM_only) else 0,
    if(exists("genes_strain_MAX_only")) length(genes_strain_MAX_only) else 0,
    if(exists("genes_strain_POST_only")) length(genes_strain_POST_only) else 0,
    if(exists("genes_strain_all_domains")) length(genes_strain_all_domains) else 0
  )
)

print(strain_pattern_summary)
write_csv(strain_pattern_summary, "strain_pattern_summary.csv")

# ===================================================================
# PART 7: GENE-LEVEL PATTERN ASSIGNMENTS
# ===================================================================

# Create comprehensive gene annotation table
all_strain_pattern_genes <- unique(c(
  if(exists("genes_main_strain")) genes_main_strain else NULL,
  if(exists("genes_domain_strain")) genes_domain_strain else NULL,
  if(exists("genes_time_strain")) genes_time_strain else NULL,
  if(exists("genes_three_way")) genes_three_way else NULL
))

strain_gene_patterns <- tibble(gene_id = all_strain_pattern_genes) |>
  left_join(all_genes_palate |> select(gene_id, symbol), by = "gene_id") |>
  mutate(
    # Main effect patterns
    CAST_high_strong = gene_id %in% if(exists("genes_CAST_high_strong")) genes_CAST_high_strong else character(0),
    CAST_high_mod = gene_id %in% if(exists("genes_CAST_high_mod")) genes_CAST_high_mod else character(0),
    B6_high_strong = gene_id %in% if(exists("genes_B6_high_strong")) genes_B6_high_strong else character(0),
    B6_high_mod = gene_id %in% if(exists("genes_B6_high_mod")) genes_B6_high_mod else character(0),
    
    # Domain × Strain patterns
    PM_biggest_diff = gene_id %in% if(exists("genes_PM_biggest_diff")) genes_PM_biggest_diff else character(0),
    MAX_biggest_diff = gene_id %in% if(exists("genes_MAX_biggest_diff")) genes_MAX_biggest_diff else character(0),
    POST_biggest_diff = gene_id %in% if(exists("genes_POST_biggest_diff")) genes_POST_biggest_diff else character(0),
    PM_CAST_high = gene_id %in% if(exists("genes_PM_CAST_high")) genes_PM_CAST_high else character(0),
    MAX_CAST_high = gene_id %in% if(exists("genes_MAX_CAST_high")) genes_MAX_CAST_high else character(0),
    POST_CAST_high = gene_id %in% if(exists("genes_POST_CAST_high")) genes_POST_CAST_high else character(0),
    Opposite_directions = gene_id %in% if(exists("genes_opposite_directions")) genes_opposite_directions else character(0),
    
    # Time × Strain patterns
    Both_inc_CAST_faster = gene_id %in% if(exists("genes_both_increase_CAST_faster")) genes_both_increase_CAST_faster else character(0),
    Both_dec_CAST_faster = gene_id %in% if(exists("genes_both_decrease_CAST_faster")) genes_both_decrease_CAST_faster else character(0),
    Both_inc_B6_faster = gene_id %in% if(exists("genes_both_increase_B6_faster")) genes_both_increase_B6_faster else character(0),
    Both_dec_B6_faster = gene_id %in% if(exists("genes_both_decrease_B6_faster")) genes_both_decrease_B6_faster else character(0),
    B6_up_CAST_down = gene_id %in% if(exists("genes_B6_up_CAST_down")) genes_B6_up_CAST_down else character(0),
    CAST_up_B6_down = gene_id %in% if(exists("genes_CAST_up_B6_down")) genes_CAST_up_B6_down else character(0),
    CAST_dynamic_B6_flat = gene_id %in% if(exists("genes_CAST_dynamic_B6_flat")) genes_CAST_dynamic_B6_flat else character(0),
    B6_dynamic_CAST_flat = gene_id %in% if(exists("genes_B6_dynamic_CAST_flat")) genes_B6_dynamic_CAST_flat else character(0),
    
    # 3-way patterns
    ThreeWay_PM_biggest = gene_id %in% if(exists("genes_PM_biggest_slope_diff")) genes_PM_biggest_slope_diff else character(0),
    ThreeWay_MAX_biggest = gene_id %in% if(exists("genes_MAX_biggest_slope_diff")) genes_MAX_biggest_slope_diff else character(0),
    ThreeWay_POST_biggest = gene_id %in% if(exists("genes_POST_biggest_slope_diff")) genes_POST_biggest_slope_diff else character(0),
    ThreeWay_consistent = gene_id %in% if(exists("genes_consistent_slope_diffs")) genes_consistent_slope_diffs else character(0),
    
    # Domain specificity
    Strain_PM_only = gene_id %in% if(exists("genes_strain_PM_only")) genes_strain_PM_only else character(0),
    Strain_MAX_only = gene_id %in% if(exists("genes_strain_MAX_only")) genes_strain_MAX_only else character(0),
    Strain_POST_only = gene_id %in% if(exists("genes_strain_POST_only")) genes_strain_POST_only else character(0),
    Strain_all_domains = gene_id %in% if(exists("genes_strain_all_domains")) genes_strain_all_domains else character(0)
  ) |>
  mutate(n_patterns = rowSums(across(CAST_high_strong:Strain_all_domains)))

write_csv(strain_gene_patterns, "strain_gene_pattern_assignments.csv")

# # Save detailed pattern tables
if(exists("main_patterns")) write_csv(main_patterns, "main_strain_effect_patterns.csv")
if(exists("domain_patterns")) write_csv(domain_patterns, "domain_strain_patterns.csv")
if(exists("time_patterns")) write_csv(time_patterns, "time_strain_patterns.csv")
if(exists("threeway_patterns")) write_csv(threeway_patterns, "threeway_patterns.csv")
if(exists("domain_specificity")) write_csv(domain_specificity, "domain_specificity_patterns.csv")

cat("\n===== RESULTS SAVED =====\n")
cat("Files created:\n")
cat("  1. strain_pattern_summary.csv\n")
cat("  2. strain_gene_pattern_assignments.csv\n")
cat("  3. main_strain_effect_patterns.csv\n")
cat("  4. domain_strain_patterns.csv\n")
cat("  5. time_strain_patterns.csv\n")
cat("  6. threeway_patterns.csv\n")
cat("  7. domain_specificity_patterns.csv\n\n")

cat("Total genes classified:", nrow(strain_gene_patterns), "\n")
cat("Genes with multiple patterns:", sum(strain_gene_patterns$n_patterns > 1), "\n")
cat("\n===== COMPLETE =====\n")


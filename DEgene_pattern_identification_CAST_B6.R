# ===================================================================
# STRAIN GENE PATTERN CLASSIFICATION
# Uses objects created in Quarto document: strain_effects_classified
# Run this AFTER the DESeq2 analysis chunks in the Quarto document
# ===================================================================

library(tidyverse)
library(here)

# ===================================================================
# PART 1: MAIN STRAIN EFFECT PATTERNS
# All genes with a significant strain effect, classified by the 
# direction and magnitude of the marginal (additive model) LFC.
# This includes genes that also have interactions — the additive
# model LFC represents the average effect across domains and time.
# Uses strain_effects_filtered (any_strain_sig & main_strain_pass).
# ===================================================================

main_patterns <- strain_effects_classified |>
  select(gene_id, symbol, main_strain_lfc, main_strain_pval, main_strain_pattern,
         has_any_interaction) |>
  filter(main_strain_pattern != "No effect / filtered") |>
  arrange(main_strain_pval)


genes_CAST_high_strong <- main_patterns |>
  filter(main_strain_pattern == "CAST >> B6 (strong)") |> pull(gene_id)
genes_CAST_high_mod <- main_patterns |>
  filter(main_strain_pattern == "CAST > B6 (moderate)") |> pull(gene_id)
genes_B6_high_strong <- main_patterns |>
  filter(main_strain_pattern == "B6 >> CAST (strong)") |> pull(gene_id)
genes_B6_high_mod <- main_patterns |>
  filter(main_strain_pattern == "B6 > CAST (moderate)") |> pull(gene_id)

# Subset: genes with consistent effects only (no interactions)
genes_main_only <- main_patterns |>
  filter(!has_any_interaction) |> pull(gene_id)

# ===================================================================
# PART 2: DOMAIN × STRAIN INTERACTION PATTERNS
# For genes where strain effects differ across domains,
# characterize the pattern of domain differences.
# Note: Genes with threeway interactions are flagged but analyzed 
# in detail in the threeway_patterns section.
# ===================================================================

domain_strain_patterns <- strain_effects_classified |>
  filter(any_strain_sig & domain_int_sig) |>
  filter(domain_int_MAX_pass | domain_int_POST_pass) |>
  mutate(
    # Get absolute magnitudes of strain effects in each domain at t=0
    abs_PM = abs(strain_effect_PM),
    abs_MAX = abs(strain_effect_MAX),
    abs_POST = abs(strain_effect_POST),
    
    # Calculate range across domains (how different are they?)
    domain_strain_range = pmax(abs_PM, abs_MAX, abs_POST) - 
      pmin(abs_PM, abs_MAX, abs_POST),
    
    # Which domain has the strongest strain effect?
    strongest_domain = case_when(
      abs_PM >= abs_MAX & abs_PM >= abs_POST ~ "PM",
      abs_MAX >= abs_POST ~ "MAX",
      TRUE ~ "POST"
    ),
    
    # Which domain has the weakest strain effect?
    weakest_domain = case_when(
      abs_PM <= abs_MAX & abs_PM <= abs_POST ~ "PM",
      abs_MAX <= abs_POST ~ "MAX",
      TRUE ~ "POST"
    ),
    
    # Are strain effects in the same direction across domains?
    same_direction = (sign(strain_effect_PM) == sign(strain_effect_MAX)) & 
      (sign(strain_effect_MAX) == sign(strain_effect_POST)),
    
    direction_pattern = case_when(
      same_direction ~ "Consistent direction",
      !same_direction ~ "Opposite directions"
    ),
    
    # Classify specificity strength
    specificity_strength = case_when(
      domain_strain_range > 1.0 ~ "High",      # One domain >> others
      domain_strain_range > 0.5 ~ "Moderate",  # Clear differences
      TRUE ~ "Low"                              # Similar magnitudes
    ),
    
    # Pattern summary
    domain_pattern = paste0(
      strongest_domain, " strongest | ",
      direction_pattern, " | ",
      specificity_strength, " specificity"
    ),
    
    # Flag time-dependent genes (detailed in other sections)
    time_status = case_when(
      threeway_int_sig ~ "Threeway (see threeway_patterns)",
      time_int_sig ~ "Main time effect present",
      TRUE ~ "Static"
    )
  ) |>
  select(gene_id, symbol,
         strain_effect_PM, strain_effect_MAX, strain_effect_POST,
         abs_PM, abs_MAX, abs_POST,
         domain_strain_range, specificity_strength,
         strongest_domain, weakest_domain, 
         direction_pattern, domain_pattern,
         time_status) |>
  arrange(desc(domain_strain_range))

# Gene lists focused on domain patterns
genes_PM_strongest_strain <- domain_strain_patterns |>
  filter(strongest_domain == "PM") |> pull(gene_id)

genes_MAX_strongest_strain <- domain_strain_patterns |>
  filter(strongest_domain == "MAX") |> pull(gene_id)

genes_POST_strongest_strain <- domain_strain_patterns |>
  filter(strongest_domain == "POST") |> pull(gene_id)

genes_opposite_directions <- domain_strain_patterns |>
  filter(direction_pattern == "Opposite directions") |> pull(gene_id)

genes_high_domain_specificity <- domain_strain_patterns |>
  filter(specificity_strength == "High") |> pull(gene_id)

# Static vs time-dependent (pointer to other analyses)
genes_domain_strain_static <- domain_strain_patterns |>
  filter(time_status == "Static") |> pull(gene_id)

genes_domain_strain_time_affected <- domain_strain_patterns |>
  filter(time_status == "Main time effect present") |> pull(gene_id)

# Summary
cat("Domain-specific strain pattern summary:\n")
table(domain_strain_patterns$strongest_domain, 
      domain_strain_patterns$specificity_strength)
cat("\nDirection consistency:\n")
table(domain_strain_patterns$direction_pattern)
cat("\nTime-dependence:\n")
table(domain_strain_patterns$time_status)

# ===================================================================
# PART 3: TIME × STRAIN INTERACTION PATTERNS (DOMAIN-AWARE)
# These genes have strain effects that change over developmental time.
# Now accounting for whether these temporal dynamics are uniform across
# domains or domain-specific.
# Source: strain_effects_classified with significant time LRT.
# ===================================================================

time_patterns <- strain_effects_classified |>
  filter(any_strain_sig & time_int_sig & time_int_pass) |>
  select(gene_id, symbol, 
         time_int_slope_avg, time_int_slope_PM, time_int_slope_lfcSE, time_int_slope_pval, time_int_pval,
         main_strain_lfc, time_pattern,
         domain_int_sig, threeway_int_sig,
         threeway_MAX_slope, threeway_POST_slope, threeway_MAX_pass, threeway_POST_pass) |>
  mutate(
    # Base trajectory (this is the PM slope when domain is in the model)
    base_slope = time_int_slope_avg,
    
    # For all cases, calculate the range of slopes across domains
    slope_PM  = if_else(threeway_int_sig, time_int_slope_PM, time_int_slope_avg),
    slope_MAX = if_else(threeway_int_sig, time_int_slope_PM + threeway_MAX_slope, time_int_slope_avg),
    slope_POST = if_else(threeway_int_sig, time_int_slope_PM + threeway_POST_slope, time_int_slope_avg),
    
    # Calculate range of slopes across domains
    slope_range = pmax(slope_PM, slope_MAX, slope_POST) - 
      pmin(slope_PM, slope_MAX, slope_POST),
    
    # Classify spatial uniformity using BOTH significance AND magnitude
    spatial_pattern = case_when(
      # Domain-specific: threeway significant AND slopes actually differ
      threeway_int_sig & (threeway_MAX_pass | threeway_POST_pass) & slope_range > 0.15 ~ 
        "Domain-specific: Different temporal dynamics across domains",
      
      # Domain-affected: domains differ at t=0 but slopes are parallel (small range)
      domain_int_sig & slope_range <= 0.15 ~ 
        "Domain-affected: Domains differ at t=0, parallel trajectories",
      
      # Uniform: no domain interaction OR very small slope differences
      !domain_int_sig & slope_range <= 0.1 ~ 
        "Uniform: Same temporal dynamics across all domains",
      
      # Edge cases: has domain/threeway sig but doesn't meet magnitude thresholds
      TRUE ~ "Mixed: Borderline pattern"
    ),
    
    # Trajectory pattern (based on base/PM slope)
    trajectory_pattern = case_when(
      base_slope > 0.2 ~ "Strong divergence (CAST increases faster)",
      base_slope > 0.1 ~ "Moderate divergence (CAST increases faster)",
      base_slope < -0.2 ~ "Strong divergence (B6 increases faster)",
      base_slope < -0.1 ~ "Moderate divergence (B6 increases faster)",
      abs(base_slope) < 0.1 ~ "Weak/parallel"
    ),
    
    # Combined detailed pattern
    detailed_pattern = paste0(trajectory_pattern, " | ", spatial_pattern)
  ) |>
  arrange(desc(abs(base_slope)))


genes_time_uniform <- time_patterns |>
  filter(!domain_int_sig & !threeway_int_sig) |> pull(gene_id)
genes_time_domain_affected <- time_patterns |>
  filter(domain_int_sig & !threeway_int_sig) |> pull(gene_id)
genes_diverging_CAST_faster <- time_patterns |>
  filter(base_slope > 0.1) |> pull(gene_id)
genes_diverging_B6_faster <- time_patterns |>
  filter(base_slope < -0.1) |> pull(gene_id)

# genes_time_domain_specific is removed — those ARE the threeway genes
# and are handled entirely in the Three-Way section


# ===================================================================
# PART 4: THREE-WAY INTERACTION PATTERNS
# These genes have domain-specific temporal strain dynamics.
# Consider BOTH initial effects (t=0) AND temporal trajectories (slopes).
# Source: strain_effects_classified with significant three-way LRT.
# ===================================================================

threeway_patterns <- strain_effects_classified |>
  filter(any_strain_sig & threeway_int_sig) |>
  filter(threeway_MAX_pass | threeway_POST_pass) |>
  select(gene_id, symbol, 
         time_int_slope_PM, threeway_MAX_slope, threeway_POST_slope,
         threeway_MAX_pass, threeway_POST_pass,
         strain_effect_PM, strain_effect_MAX, strain_effect_POST) |>
  mutate(
    # Calculate temporal slopes for each domain
    slope_PM = time_int_slope_PM,
    slope_MAX = time_int_slope_PM + threeway_MAX_slope,
    slope_POST = time_int_slope_PM + threeway_POST_slope,
    
    # Initial effects at t=0
    init_PM = strain_effect_PM,
    init_MAX = strain_effect_MAX,
    init_POST = strain_effect_POST,
    
    # Magnitude of initial differences across domains (row-wise)
    init_range = pmax(init_PM, init_MAX, init_POST) - 
      pmin(init_PM, init_MAX, init_POST),
    
    # Magnitude of slope differences across domains (row-wise)
    slope_range = pmax(slope_PM, slope_MAX, slope_POST) - 
      pmin(slope_PM, slope_MAX, slope_POST),
    
    # Classify based on relative importance of initial vs dynamic differences
    pattern_type = case_when(
      init_range > 1 & slope_range > 0.8 ~ "Complex: Differ in both initial state & dynamics",
      init_range > 1 & slope_range <= 0.8 ~ "Primarily initial: Domains differ at t=0, similar trajectories",
      init_range <= 1 & slope_range > 0.8 ~ "Primarily dynamic: Similar at t=0, divergent trajectories",
      TRUE ~ "Subtle: Small differences in both"
    ),
    
    # Identify which domain is most different (combining initial + slope)
    # Use combined score: absolute deviation from mean in both dimensions
    combined_PM = abs(init_PM - mean(c(init_PM, init_MAX, init_POST))) + 
      abs(slope_PM - mean(c(slope_PM, slope_MAX, slope_POST))),
    combined_MAX = abs(init_MAX - mean(c(init_PM, init_MAX, init_POST))) + 
      abs(slope_MAX - mean(c(slope_PM, slope_MAX, slope_POST))),
    combined_POST = abs(init_POST - mean(c(init_PM, init_MAX, init_POST))) + 
      abs(slope_POST - mean(c(slope_PM, slope_MAX, slope_POST))),
    
    most_distinct_domain = case_when(
      combined_PM > combined_MAX & combined_PM > combined_POST ~ "PM",
      combined_MAX > combined_POST ~ "MAX",
      TRUE ~ "POST"
    ),
    
    # Check if slopes are consistent in direction
    slope_direction_consistent = (sign(slope_PM) == sign(slope_MAX) & 
                                    sign(slope_MAX) == sign(slope_POST)),
    
    detailed_pattern = paste0(pattern_type, " | ", most_distinct_domain, " most distinct")
  ) |>
  arrange(desc(slope_range + init_range))

genes_threeway_complex <- threeway_patterns |>
  filter(str_detect(pattern_type, "Complex")) |> pull(gene_id)
genes_threeway_initial <- threeway_patterns |>
  filter(str_detect(pattern_type, "Primarily initial")) |> pull(gene_id)
genes_threeway_dynamic <- threeway_patterns |>
  filter(str_detect(pattern_type, "Primarily dynamic")) |> pull(gene_id)


# ===================================================================
# PART 5: SUMMARY TABLE
# ===================================================================

pattern_summary <- tibble(
  Category = c(
    rep("Main Strain Effect", 5),
    rep("Domain × Strain", 8),
    rep("Time × Strain", 4),
    rep("Three-Way Interaction", 3)
  ),
  Pattern = c(
    "CAST >> B6 (strong)", "CAST > B6 (moderate)",
    "B6 >> CAST (strong)", "B6 > CAST (moderate)",
    "Consistent only (no interactions)",
    "PM strongest strain effect", "MAX strongest strain effect", 
    "POST strongest strain effect", "Opposite directions across domains",
    "High domain specificity (range > 1.0)",
    "Static (no time interactions)", "Time-affected (main time effect)", 
    "Time-varying (threeway interaction)",
    "Diverging (CAST faster)", "Diverging (B6 faster)",
    "Temporal dynamics: Uniform across domains",
    "Temporal dynamics: Domain-affected (parallel slopes)",
    "Complex (differ in initial & dynamics)", 
    "Primarily initial differences", 
    "Primarily dynamic differences"
  ),
  N_Genes = c(
    length(genes_CAST_high_strong), length(genes_CAST_high_mod),
    length(genes_B6_high_strong), length(genes_B6_high_mod),
    length(genes_main_only),
    length(genes_PM_strongest_strain), length(genes_MAX_strongest_strain),
    length(genes_POST_strongest_strain), length(genes_opposite_directions),
    length(genes_high_domain_specificity),
    length(genes_domain_strain_static), length(genes_domain_strain_time_affected),
    nrow(domain_strain_patterns |> filter(time_status == "Threeway (see threeway_patterns)")),
    length(genes_diverging_CAST_faster), length(genes_diverging_B6_faster),
    length(genes_time_uniform), length(genes_time_domain_affected),
    length(genes_threeway_complex), length(genes_threeway_initial),
    length(genes_threeway_dynamic)
  ),
  Description = c(
    "All strain-affected genes, marginal LFC > 1", 
    "All strain-affected genes, 0.5 < marginal LFC ≤ 1",
    "All strain-affected genes, marginal LFC < -1", 
    "All strain-affected genes, -1 ≤ marginal LFC < -0.5",
    "Subset with no domain, time, or three-way interactions",
    "Domain int. sig, PM has strongest |strain effect|", 
    "Domain int. sig, MAX has strongest |strain effect|",
    "Domain int. sig, POST has strongest |strain effect|", 
    "Domain int. sig, strain effect direction varies across domains",
    "Domain int. sig, strain effect range > 1.0 log2FC across domains",
    "Domain int. sig, no time/3-way interaction",
    "Domain int. sig, has time interaction (may shift pattern)",
    "Domain int. sig, has 3-way interaction (pattern changes over time)",
    "Time int. sig, positive slope (CAST increases faster)", 
    "Time int. sig, negative slope (B6 increases faster)",
    "Time int. sig, no domain/3-way interaction",
    "Time int. sig, has domain interaction (parallel slopes)",
    "3-way: Domains differ in both t=0 state and slope",
    "3-way: Domains differ mainly at t=0, similar slopes",
    "3-way: Similar at t=0, divergent slopes"
  )
)

print(pattern_summary)

# ===================================================================
# PART 6: GENE ASSIGNMENT TABLE
# ===================================================================

all_strain_genes <- strain_effects_classified |>
  filter(any_strain_sig) |>
  select(gene_id, symbol, main_strain_lfc, main_strain_pval, 
         main_strain_pattern, has_any_interaction,
         domain_int_sig, time_int_sig, threeway_int_sig)

gene_assignments <- all_strain_genes |>
  mutate(
    # Main effects (all strain-affected genes, from additive model)
    CAST_high_strong = gene_id %in% genes_CAST_high_strong,
    CAST_high_mod = gene_id %in% genes_CAST_high_mod,
    B6_high_strong = gene_id %in% genes_B6_high_strong,
    B6_high_mod = gene_id %in% genes_B6_high_mod,
    Main_effect_only = gene_id %in% genes_main_only,
    
    # Domain × strain interaction patterns
    PM_strongest_strain = gene_id %in% genes_PM_strongest_strain,
    MAX_strongest_strain = gene_id %in% genes_MAX_strongest_strain,
    POST_strongest_strain = gene_id %in% genes_POST_strongest_strain,
    Opposite_directions = gene_id %in% genes_opposite_directions,
    High_domain_specificity = gene_id %in% genes_high_domain_specificity,
    Domain_strain_static = gene_id %in% genes_domain_strain_static,
    Domain_strain_time_affected = gene_id %in% genes_domain_strain_time_affected,
    
    # Time × strain interaction patterns
    Diverging_CAST_faster = gene_id %in% genes_diverging_CAST_faster,
    Diverging_B6_faster = gene_id %in% genes_diverging_B6_faster,
    Time_uniform = gene_id %in% genes_time_uniform,
    Time_domain_affected = gene_id %in% genes_time_domain_affected,
    
    # Three-way patterns
    ThreeWay_complex = gene_id %in% genes_threeway_complex,
    ThreeWay_initial = gene_id %in% genes_threeway_initial,
    ThreeWay_dynamic = gene_id %in% genes_threeway_dynamic
  )


# ===================================================================
# PART 7: SAVE OUTPUTS
# ===================================================================

write_csv(pattern_summary, here("strain_pattern_summary.csv"))
write_csv(main_patterns, here("main_strain_effect_patterns.csv"))
write_csv(domain_strain_patterns, here("domain_strain_patterns.csv"))
write_csv(time_patterns, here("time_strain_patterns.csv"))
write_csv(threeway_patterns, here("threeway_patterns.csv"))
write_csv(gene_assignments, here("strain_gene_pattern_assignments.csv"))
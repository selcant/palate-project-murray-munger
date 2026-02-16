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
# PART 2: DOMAIN × STRAIN INTERACTION PATTERNS (TIME-AWARE)
# These genes have strain effects that differ across palate domains.
# Now accounting for whether these domain differences are static or
# change over developmental time.
# Source: strain_effects_classified with significant domain LRT.
# ===================================================================

domain_patterns <- strain_effects_classified |>
  filter(any_strain_sig & domain_int_sig) |>
  filter(domain_int_MAX_pass | domain_int_POST_pass) |>
  select(gene_id, symbol, 
         strain_effect_PM, strain_effect_MAX, strain_effect_POST,
         strain_effect_PM_pass, strain_effect_MAX_pass, strain_effect_POST_pass,
         strain_effect_PM_pval, strain_effect_MAX_pval, strain_effect_POST_pval,
         domain_int_MAX_lfc, domain_int_POST_lfc,
         domain_int_MAX_pass, domain_int_POST_pass,
         max_domain_effect, max_domain_direction, domain_consistency,
         time_int_sig, threeway_int_sig) |>
  mutate(
    # Classify temporal stability of domain pattern
    temporal_stability = case_when(
      threeway_int_sig ~ "Dynamic: Domain pattern changes over time (3-way)",
      time_int_sig ~ "Dynamic: Main time effect present, domain pattern may shift",
      TRUE ~ "Static: Domain pattern stable over time"
    ),
    
    # Create detailed pattern combining direction, consistency, and stability
    detailed_pattern = paste0(
      max_domain_direction, " | ", 
      domain_consistency, " | ",
      case_when(
        threeway_int_sig ~ "Time-varying",
        time_int_sig ~ "Time-affected", 
        TRUE ~ "Static"
      )
    ),
    
    # Magnitude of domain differences at t=0
    domain_diff_magnitude = abs(domain_int_MAX_lfc) + abs(domain_int_POST_lfc)
  ) |>
  arrange(desc(domain_diff_magnitude))

# Split by temporal stability
genes_domain_static <- domain_patterns |>
  filter(!time_int_sig & !threeway_int_sig) |> pull(gene_id)
genes_domain_time_affected <- domain_patterns |>
  filter(time_int_sig & !threeway_int_sig) |> pull(gene_id)
genes_domain_time_varying <- domain_patterns |>
  filter(threeway_int_sig) |> pull(gene_id)

# Original categories (at t=0)
genes_PM_biggest_effect <- domain_patterns |>
  filter(max_domain_effect == "PM") |> pull(gene_id)
genes_MAX_biggest_effect <- domain_patterns |>
  filter(max_domain_effect == "MAX") |> pull(gene_id)
genes_POST_biggest_effect <- domain_patterns |>
  filter(max_domain_effect == "POST") |> pull(gene_id)
genes_opposite_directions <- domain_patterns |>
  filter(domain_consistency == "Variable") |> pull(gene_id)

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
         time_int_slope, time_int_slope_lfcSE, time_int_slope_pval, time_int_pval,
         main_strain_lfc, time_pattern,
         domain_int_sig, threeway_int_sig,
         threeway_MAX_slope, threeway_POST_slope, threeway_MAX_pass, threeway_POST_pass) |>
  mutate(
    # Base trajectory (this is the PM slope when domain is in the model)
    base_slope = time_int_slope,
    
    # Classify spatial uniformity of temporal dynamics
    spatial_pattern = case_when(
      threeway_int_sig & (threeway_MAX_pass | threeway_POST_pass) ~ 
        "Domain-specific: Different temporal dynamics across domains",
      domain_int_sig ~ 
        "Domain-affected: Domains differ at t=0, parallel trajectories",
      TRUE ~ 
        "Uniform: Same temporal dynamics across all domains"
    ),
    
    # For domain-specific cases, calculate the range of slopes
    slope_PM = time_int_slope,
    slope_MAX = if_else(threeway_int_sig, time_int_slope + threeway_MAX_slope, time_int_slope),
    slope_POST = if_else(threeway_int_sig, time_int_slope + threeway_POST_slope, time_int_slope),
    slope_range = max(abs(c(slope_PM, slope_MAX, slope_POST))) - 
      min(abs(c(slope_PM, slope_MAX, slope_POST))),
    
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

# Split by spatial pattern
genes_time_uniform <- time_patterns |>
  filter(!domain_int_sig & !threeway_int_sig) |> pull(gene_id)
genes_time_domain_affected <- time_patterns |>
  filter(domain_int_sig & !threeway_int_sig) |> pull(gene_id)
genes_time_domain_specific <- time_patterns |>
  filter(threeway_int_sig) |> pull(gene_id)

# Original categories (based on base slope)
genes_diverging_CAST_faster <- time_patterns |>
  filter(base_slope > 0.1) |> pull(gene_id)
genes_diverging_B6_faster <- time_patterns |>
  filter(base_slope < -0.1) |> pull(gene_id)

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
         time_int_slope, threeway_MAX_slope, threeway_POST_slope,
         threeway_MAX_pass, threeway_POST_pass,
         strain_effect_PM, strain_effect_MAX, strain_effect_POST) |>
  mutate(
    # Calculate temporal slopes for each domain
    slope_PM = time_int_slope,
    slope_MAX = time_int_slope + threeway_MAX_slope,
    slope_POST = time_int_slope + threeway_POST_slope,
    
    # Initial effects at t=0
    init_PM = strain_effect_PM,
    init_MAX = strain_effect_MAX,
    init_POST = strain_effect_POST,
    
    # Magnitude of initial differences across domains
    init_range = max(abs(c(init_PM, init_MAX, init_POST))) - 
      min(abs(c(init_PM, init_MAX, init_POST))),
    
    # Magnitude of slope differences across domains
    slope_range = max(abs(c(slope_PM, slope_MAX, slope_POST))) - 
      min(abs(c(slope_PM, slope_MAX, slope_POST))),
    
    # Classify based on relative importance of initial vs dynamic differences
    pattern_type = case_when(
      init_range > 0.5 & slope_range > 0.15 ~ "Complex: Differ in both initial state & dynamics",
      init_range > 0.5 & slope_range <= 0.15 ~ "Primarily initial: Domains differ at t=0, similar trajectories",
      init_range <= 0.5 & slope_range > 0.15 ~ "Primarily dynamic: Similar at t=0, divergent trajectories",
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
# PART 5: DOMAIN SPECIFICITY PATTERNS (TIME-AWARE)
# For genes with domain × strain interactions, classify domain specificity
# while accounting for temporal dynamics when present.
# ===================================================================

domain_specificity <- strain_effects_classified |>
  filter(any_strain_sig & domain_int_sig) |>
  mutate(
    # Check if each domain has a significant AND well-estimated effect at t=0
    effect_PM_t0 = strain_effect_PM_pass,
    effect_MAX_t0 = strain_effect_MAX_pass,
    effect_POST_t0 = strain_effect_POST_pass,
    
    # For genes with time interactions, also consider if effect persists/grows over time
    # If time interaction exists, check if domain-specific slopes are significant
    has_time_int = time_int_sig | threeway_int_sig,
    
    # Simplified: Does the domain show a clear strain effect?
    # Consider both: (1) effect at t=0, or (2) strong temporal dynamics
    # For now, we'll use the t=0 effects but flag time-dynamic genes
    
    specificity_pattern = case_when(
      effect_PM_t0 & !effect_MAX_t0 & !effect_POST_t0 ~ "PM only",
      !effect_PM_t0 & effect_MAX_t0 & !effect_POST_t0 ~ "MAX only",
      !effect_PM_t0 & !effect_MAX_t0 & effect_POST_t0 ~ "POST only",
      effect_PM_t0 & effect_MAX_t0 & !effect_POST_t0 ~ "PM & MAX",
      effect_PM_t0 & !effect_MAX_t0 & effect_POST_t0 ~ "PM & POST",
      !effect_PM_t0 & effect_MAX_t0 & effect_POST_t0 ~ "MAX & POST",
      effect_PM_t0 & effect_MAX_t0 & effect_POST_t0 ~ "All domains",
      TRUE ~ "None (filtered)"
    ),
    
    # Add time-awareness qualifier
    specificity_with_time = if_else(
      has_time_int,
      paste0(specificity_pattern, " (time-varying)"),
      specificity_pattern
    )
  ) |>
  select(gene_id, symbol, 
         strain_effect_PM, strain_effect_MAX, strain_effect_POST,
         effect_PM_t0, effect_MAX_t0, effect_POST_t0, 
         has_time_int, specificity_pattern, specificity_with_time)

# Gene lists based on static domain specificity (at t=0)
genes_strain_PM_only <- domain_specificity |>
  filter(specificity_pattern == "PM only") |> pull(gene_id)
genes_strain_MAX_only <- domain_specificity |>
  filter(specificity_pattern == "MAX only") |> pull(gene_id)
genes_strain_POST_only <- domain_specificity |>
  filter(specificity_pattern == "POST only") |> pull(gene_id)
genes_strain_all_domains <- domain_specificity |>
  filter(specificity_pattern == "All domains") |> pull(gene_id)

# Additional: genes where specificity is time-dependent
genes_strain_domain_specific_static <- domain_specificity |>
  filter(!has_time_int & specificity_pattern != "All domains" & 
           specificity_pattern != "None (filtered)") |> pull(gene_id)
genes_strain_domain_specific_dynamic <- domain_specificity |>
  filter(has_time_int & specificity_pattern != "All domains" & 
           specificity_pattern != "None (filtered)") |> pull(gene_id)

# ===================================================================
# PART 6: SUMMARY TABLE (revised categories)
# ===================================================================

pattern_summary <- tibble(
  Category = c(
    rep("Main Strain Effect", 5),
    rep("Domain × Strain", 7),
    rep("Time × Strain", 5),
    rep("Three-Way Interaction", 3),
    rep("Domain Specificity", 6)
  ),
  Pattern = c(
    "CAST >> B6 (strong)", "CAST > B6 (moderate)",
    "B6 >> CAST (strong)", "B6 > CAST (moderate)",
    "Consistent only (no interactions)",
    "Biggest effect in PM", "Biggest effect in MAX", 
    "Biggest effect in POST", "Opposite directions",
    "Domain pattern: Static", "Domain pattern: Time-affected", 
    "Domain pattern: Time-varying (3-way)",
    "Diverging (CAST faster)", "Diverging (B6 faster)",
    "Temporal dynamics: Uniform across domains",
    "Temporal dynamics: Domain-affected (parallel slopes)",
    "Temporal dynamics: Domain-specific (different slopes)",
    "Complex (differ in initial & dynamics)", 
    "Primarily initial differences", 
    "Primarily dynamic differences",
    "PM only", "MAX only", "POST only", "All domains",
    "Domain-specific (static over time)", "Domain-specific (time-varying)"
  ),
  N_Genes = c(
    length(genes_CAST_high_strong), length(genes_CAST_high_mod),
    length(genes_B6_high_strong), length(genes_B6_high_mod),
    length(genes_main_only),
    length(genes_PM_biggest_effect), length(genes_MAX_biggest_effect),
    length(genes_POST_biggest_effect), length(genes_opposite_directions),
    length(genes_domain_static), length(genes_domain_time_affected),
    length(genes_domain_time_varying),
    length(genes_diverging_CAST_faster), length(genes_diverging_B6_faster),
    length(genes_time_uniform), length(genes_time_domain_affected),
    length(genes_time_domain_specific),
    length(genes_threeway_complex), length(genes_threeway_initial),
    length(genes_threeway_dynamic),
    length(genes_strain_PM_only), length(genes_strain_MAX_only),
    length(genes_strain_POST_only), length(genes_strain_all_domains),
    length(genes_strain_domain_specific_static), 
    length(genes_strain_domain_specific_dynamic)
  ),
  Description = c(
    "All strain-affected genes, marginal LFC > 1", 
    "All strain-affected genes, 0.5 < marginal LFC ≤ 1",
    "All strain-affected genes, marginal LFC < -1", 
    "All strain-affected genes, -1 ≤ marginal LFC < -0.5",
    "Subset with no domain, time, or three-way interactions",
    "Domain int. sig, largest |effect| in PM", 
    "Domain int. sig, largest |effect| in MAX",
    "Domain int. sig, largest |effect| in POST", 
    "Domain int. sig, effect direction varies",
    "Domain int. sig, no time/3-way interaction",
    "Domain int. sig, has time interaction (may shift)",
    "Domain int. sig, has 3-way interaction (changes over time)",
    "Time int. sig, positive slope", 
    "Time int. sig, negative slope",
    "Time int. sig, no domain/3-way interaction",
    "Time int. sig, has domain interaction (parallel slopes)",
    "Time int. sig, has 3-way interaction (domain-specific slopes)",
    "3-way: Domains differ in both t=0 state and slope",
    "3-way: Domains differ mainly at t=0, similar slopes",
    "3-way: Similar at t=0, divergent slopes",
    "Domain int. sig, strain effect only in PM at t=0", 
    "Domain int. sig, strain effect only in MAX at t=0",
    "Domain int. sig, strain effect only in POST at t=0", 
    "Domain int. sig, strain effect in all domains at t=0",
    "Domain-specific effect, no time interaction",
    "Domain-specific effect, changes over time"
  )
)

print(pattern_summary)

# ===================================================================
# PART 7: GENE ASSIGNMENT TABLE (add new categories)
# ===================================================================

gene_assignments <- all_strain_genes |>
  mutate(
    # Main effects (all strain-affected genes, from additive model)
    CAST_high_strong = gene_id %in% genes_CAST_high_strong,
    CAST_high_mod = gene_id %in% genes_CAST_high_mod,
    B6_high_strong = gene_id %in% genes_B6_high_strong,
    B6_high_mod = gene_id %in% genes_B6_high_mod,
    Main_effect_only = gene_id %in% genes_main_only,
    
    # Domain interaction patterns
    PM_biggest_effect = gene_id %in% genes_PM_biggest_effect,
    MAX_biggest_effect = gene_id %in% genes_MAX_biggest_effect,
    POST_biggest_effect = gene_id %in% genes_POST_biggest_effect,
    Opposite_directions = gene_id %in% genes_opposite_directions,
    Domain_static = gene_id %in% genes_domain_static,
    Domain_time_affected = gene_id %in% genes_domain_time_affected,
    Domain_time_varying = gene_id %in% genes_domain_time_varying,
    
    # Time interaction patterns
    Diverging_CAST_faster = gene_id %in% genes_diverging_CAST_faster,
    Diverging_B6_faster = gene_id %in% genes_diverging_B6_faster,
    Time_uniform = gene_id %in% genes_time_uniform,
    Time_domain_affected = gene_id %in% genes_time_domain_affected,
    Time_domain_specific = gene_id %in% genes_time_domain_specific,
    
    # Three-way patterns (updated)
    ThreeWay_complex = gene_id %in% genes_threeway_complex,
    ThreeWay_initial = gene_id %in% genes_threeway_initial,
    ThreeWay_dynamic = gene_id %in% genes_threeway_dynamic,
    
    # Domain specificity (updated)
    Strain_PM_only = gene_id %in% genes_strain_PM_only,
    Strain_MAX_only = gene_id %in% genes_strain_MAX_only,
    Strain_POST_only = gene_id %in% genes_strain_POST_only,
    Strain_all_domains = gene_id %in% genes_strain_all_domains,
    Domain_specific_static = gene_id %in% genes_strain_domain_specific_static,
    Domain_specific_dynamic = gene_id %in% genes_strain_domain_specific_dynamic
  )


# ===================================================================
# PART 8: SAVE OUTPUTS
# ===================================================================

write_csv(pattern_summary, here("strain_pattern_summary.csv"))
write_csv(main_patterns, here("main_strain_effect_patterns.csv"))
write_csv(domain_patterns, here("domain_strain_patterns.csv"))
write_csv(time_patterns, here("time_strain_patterns.csv"))
write_csv(threeway_patterns, here("threeway_patterns.csv"))
write_csv(domain_specificity, here("domain_specificity_patterns.csv"))
write_csv(gene_assignments, here("strain_gene_pattern_assignments.csv"))
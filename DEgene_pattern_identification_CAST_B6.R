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
         has_any_strain_interact) |>
  filter(main_strain_pattern != "No effect / filtered") |>
  arrange(main_strain_pval)

genes_CAST_high_strong <- main_patterns |> filter(main_strain_pattern == "CAST >> B6 (strong)")   |> pull(gene_id)
genes_CAST_high_mod    <- main_patterns |> filter(main_strain_pattern == "CAST > B6 (moderate)")  |> pull(gene_id)
genes_B6_high_strong   <- main_patterns |> filter(main_strain_pattern == "B6 >> CAST (strong)")   |> pull(gene_id)
genes_B6_high_mod      <- main_patterns |> filter(main_strain_pattern == "B6 > CAST (moderate)")  |> pull(gene_id)

# Subset: genes with consistent effects only (no strain interactions)
genes_main_only <- main_patterns |> filter(!has_any_strain_interact) |> pull(gene_id)


# ===================================================================
# PART 2: DOMAIN MAIN EFFECT PATTERNS
# Genes where expression differs across AP domains, with no
# domain×strain or domain×time interactions (i.e., the domain
# effect is consistent across strains and developmental time).
# Uses passes_domain_filter.
# ===================================================================

domain_main_patterns <- strain_effects_classified |>
  filter(passes_domain_filter) |>
  select(gene_id, symbol,
         domain_main_MAX_vs_PM_lfc, domain_main_POST_vs_PM_lfc, domain_main_MAX_vs_POST_lfc,
         domain_main_MAX_pass, domain_main_POST_pass, domain_main_MAX_vs_POST_pass,
         domain_main_pattern, domain_pval) |>
  mutate(
    # Overall domain range (how different is expression across domains?)
    domain_lfc_range = pmax(abs(domain_main_MAX_vs_PM_lfc),
                            abs(domain_main_POST_vs_PM_lfc),
                            abs(domain_main_MAX_vs_POST_lfc),
                            na.rm = TRUE),
    
    # Direction of differences
    anterior_posterior_gradient = case_when(
      domain_main_MAX_vs_PM_lfc  > 0 & domain_main_POST_vs_PM_lfc  > 0 ~ "Increasing posterior (PM lowest)",
      domain_main_MAX_vs_PM_lfc  < 0 & domain_main_POST_vs_PM_lfc  < 0 ~ "Decreasing posterior (PM highest)",
      domain_main_MAX_vs_PM_lfc  > 0 & domain_main_POST_vs_PM_lfc  < 0 ~ "MAX peak (non-monotone)",
      domain_main_MAX_vs_PM_lfc  < 0 & domain_main_POST_vs_PM_lfc  > 0 ~ "PM/POST peak (non-monotone)",
      TRUE                                                               ~ "Complex"
    ),
    
    domain_specificity = case_when(
      domain_lfc_range > 1.0 ~ "High",
      domain_lfc_range > 0.5 ~ "Moderate",
      TRUE                   ~ "Low"
    )
  ) |>
  arrange(domain_pval)

genes_domain_increasing_posterior <- domain_main_patterns |> filter(anterior_posterior_gradient == "Increasing posterior (PM lowest)")  |> pull(gene_id)
genes_domain_decreasing_posterior <- domain_main_patterns |> filter(anterior_posterior_gradient == "Decreasing posterior (PM highest)") |> pull(gene_id)
genes_domain_MAX_peak             <- domain_main_patterns |> filter(anterior_posterior_gradient == "MAX peak (non-monotone)")           |> pull(gene_id)
genes_domain_high_specificity     <- domain_main_patterns |> filter(domain_specificity == "High")                                      |> pull(gene_id)

cat("Domain main effect gradient summary:\n")
print(table(domain_main_patterns$anterior_posterior_gradient,
            domain_main_patterns$domain_specificity))


# ===================================================================
# PART 3: TIME MAIN EFFECT PATTERNS
# Genes where expression changes across developmental time, with no
# time×strain or domain×time interactions (i.e., the trajectory is
# shared across strains and AP domains).
# Uses passes_time_filter.
# ===================================================================

time_main_patterns <- strain_effects_classified |>
  filter(passes_time_filter) |>
  select(gene_id, symbol,
         time_main_slope, time_main_slope_lfcSE, time_main_slope_pval,
         time_main_pattern, time_pval) |>
  mutate(
    slope_magnitude = case_when(
      abs(time_main_slope) > 0.5 ~ "Strong",
      abs(time_main_slope) > 0.2 ~ "Moderate",
      TRUE                       ~ "Weak"
    )
  ) |>
  arrange(time_pval)

genes_time_increasing <- time_main_patterns |> filter(time_main_pattern == "Expression increases over time") |> pull(gene_id)
genes_time_decreasing <- time_main_patterns |> filter(time_main_pattern == "Expression decreases over time") |> pull(gene_id)

cat("Time main effect pattern summary:\n")
print(table(time_main_patterns$time_main_pattern,
            time_main_patterns$slope_magnitude))


# ===================================================================
# PART 4: DOMAIN × TIME INTERACTION PATTERNS
# Genes where the temporal trajectory differs across AP domains,
# but this is not further modified by strain. Characterises which
# domains have different slopes and how those slopes compare.
# Uses passes_domain_time_filter.
# ===================================================================

domain_time_patterns <- strain_effects_classified |>
  filter(passes_domain_time_filter) |>
  select(gene_id, symbol,
         time_slope_PM, time_slope_MAX, time_slope_POST,
         domain_time_int_MAX_slope, domain_time_int_POST_slope,
         domain_time_int_MAX_pass, domain_time_int_POST_pass,
         domain_time_pattern, domain_time_int_pval) |>
  mutate(
    # Which domains show significant slope differences vs PM?
    diverging_domains = case_when(
      domain_time_int_MAX_pass & domain_time_int_POST_pass ~ "MAX & POST",
      domain_time_int_MAX_pass                             ~ "MAX only",
      domain_time_int_POST_pass                            ~ "POST only",
      TRUE                                                 ~ "None (filtered)"
    ),
    
    # Are domains converging or diverging over time?
    # (i.e., does the interaction slope point away from or toward zero?)
    MAX_dynamics = case_when(
      !domain_time_int_MAX_pass          ~ "No MAX difference",
      domain_time_int_MAX_slope >  0.1   ~ "MAX increases faster than PM",
      domain_time_int_MAX_slope < -0.1   ~ "MAX increases slower than PM",
      TRUE                               ~ "Parallel"
    ),
    POST_dynamics = case_when(
      !domain_time_int_POST_pass         ~ "No POST difference",
      domain_time_int_POST_slope >  0.1  ~ "POST increases faster than PM",
      domain_time_int_POST_slope < -0.1  ~ "POST increases slower than PM",
      TRUE                               ~ "Parallel"
    ),
    
    # Overall slope range across domains
    slope_range = pmax(time_slope_PM, time_slope_MAX, time_slope_POST, na.rm = TRUE) -
      pmin(time_slope_PM, time_slope_MAX, time_slope_POST, na.rm = TRUE),
    
    slope_range_magnitude = case_when(
      slope_range > 0.5 ~ "High",
      slope_range > 0.2 ~ "Moderate",
      TRUE              ~ "Low"
    )
  ) |>
  arrange(domain_time_int_pval)

genes_domain_time_MAX_diverges  <- domain_time_patterns |> filter(domain_time_int_MAX_pass)                             |> pull(gene_id)
genes_domain_time_POST_diverges <- domain_time_patterns |> filter(domain_time_int_POST_pass)                            |> pull(gene_id)
genes_domain_time_both_diverge  <- domain_time_patterns |> filter(domain_time_int_MAX_pass & domain_time_int_POST_pass) |> pull(gene_id)

cat("Domain × Time pattern summary:\n")
print(table(domain_time_patterns$diverging_domains,
            domain_time_patterns$slope_range_magnitude))


# ===================================================================
# PART 5: DOMAIN × STRAIN INTERACTION PATTERNS
# For genes where strain effects differ across domains,
# characterize the pattern of domain differences.
# Note: Genes with three-way interactions are flagged but analysed
# in detail in Part 7.
# ===================================================================

domain_strain_patterns <- strain_effects_classified |>
  filter(any_strain_sig & strain_domain_int_sig) |>
  filter(domain_int_MAX_pass | domain_int_POST_pass) |>
  mutate(
    abs_PM   = abs(strain_effect_PM),
    abs_MAX  = abs(strain_effect_MAX),
    abs_POST = abs(strain_effect_POST),
    
    domain_strain_range = pmax(abs_PM, abs_MAX, abs_POST) -
      pmin(abs_PM, abs_MAX, abs_POST),
    
    strongest_domain = case_when(
      abs_PM >= abs_MAX & abs_PM >= abs_POST ~ "PM",
      abs_MAX >= abs_POST                    ~ "MAX",
      TRUE                                   ~ "POST"
    ),
    
    weakest_domain = case_when(
      abs_PM <= abs_MAX & abs_PM <= abs_POST ~ "PM",
      abs_MAX <= abs_POST                    ~ "MAX",
      TRUE                                   ~ "POST"
    ),
    
    same_direction = (sign(strain_effect_PM) == sign(strain_effect_MAX)) &
      (sign(strain_effect_MAX) == sign(strain_effect_POST)),
    
    direction_pattern = if_else(same_direction, "Consistent direction", "Opposite directions"),
    
    specificity_strength = case_when(
      domain_strain_range > 1.0 ~ "High",
      domain_strain_range > 0.5 ~ "Moderate",
      TRUE                      ~ "Low"
    ),
    
    domain_pattern = paste0(strongest_domain, " strongest | ",
                            direction_pattern, " | ",
                            specificity_strength, " specificity"),
    
    time_status = case_when(
      threeway_int_sig    ~ "Three-way (see threeway_patterns)",
      strain_time_int_sig ~ "Main time effect present",
      TRUE                ~ "Static"
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

genes_PM_strongest_strain    <- domain_strain_patterns |> filter(strongest_domain    == "PM")                   |> pull(gene_id)
genes_MAX_strongest_strain   <- domain_strain_patterns |> filter(strongest_domain    == "MAX")                  |> pull(gene_id)
genes_POST_strongest_strain  <- domain_strain_patterns |> filter(strongest_domain    == "POST")                 |> pull(gene_id)
genes_opposite_directions    <- domain_strain_patterns |> filter(direction_pattern   == "Opposite directions")  |> pull(gene_id)
genes_high_domain_specificity <- domain_strain_patterns |> filter(specificity_strength == "High")               |> pull(gene_id)
genes_domain_strain_static        <- domain_strain_patterns |> filter(time_status == "Static")                              |> pull(gene_id)
genes_domain_strain_time_affected <- domain_strain_patterns |> filter(time_status == "Main time effect present")            |> pull(gene_id)

cat("Domain-specific strain pattern summary:\n")
print(table(domain_strain_patterns$strongest_domain,
            domain_strain_patterns$specificity_strength))
cat("\nDirection consistency:\n")
print(table(domain_strain_patterns$direction_pattern))
cat("\nTime-dependence:\n")
print(table(domain_strain_patterns$time_status))


# ===================================================================
# PART 6: TIME × STRAIN INTERACTION PATTERNS (DOMAIN-AWARE)
# Genes with strain effects that change over developmental time.
# Accounts for whether temporal dynamics are uniform across domains
# or domain-specific (flagged for Part 7 if three-way sig).
# ===================================================================

time_strain_patterns <- strain_effects_classified |>
  filter(any_strain_sig & strain_time_int_sig & time_int_pass) |>
  select(gene_id, symbol,
         time_int_slope_avg, time_strain_slope_PM, time_int_slope_lfcSE, time_int_slope_pval,
         main_strain_lfc, time_strain_pattern,
         strain_domain_int_sig, threeway_int_sig,strain_time_int_pval,
         threeway_MAX_slope, threeway_POST_slope, threeway_MAX_pass, threeway_POST_pass) |>
  mutate(
    base_slope = time_int_slope_avg,
    
    # Domain-specific slopes (use full-model PM slope when 3-way is present)
    slope_PM   = if_else(threeway_int_sig, time_strain_slope_PM,                             time_int_slope_avg),
    slope_MAX  = if_else(threeway_int_sig, time_strain_slope_PM + threeway_MAX_slope,         time_int_slope_avg),
    slope_POST = if_else(threeway_int_sig, time_strain_slope_PM + threeway_POST_slope,        time_int_slope_avg),
    
    slope_range = pmax(slope_PM, slope_MAX, slope_POST) -
      pmin(slope_PM, slope_MAX, slope_POST),
    
    spatial_pattern = case_when(
      threeway_int_sig & (threeway_MAX_pass | threeway_POST_pass) & slope_range > 0.15 ~
        "Domain-specific: Different temporal dynamics across domains",
      strain_domain_int_sig & slope_range <= 0.15 ~
        "Domain-affected: Domains differ at t=0, parallel trajectories",
      !strain_domain_int_sig & slope_range <= 0.1 ~
        "Uniform: Same temporal dynamics across all domains",
      TRUE ~ "Mixed: Borderline pattern"
    ),
    
    trajectory_pattern = case_when(
      base_slope >  0.2  ~ "Strong divergence (CAST increases faster)",
      base_slope >  0.1  ~ "Moderate divergence (CAST increases faster)",
      base_slope < -0.2  ~ "Strong divergence (B6 increases faster)",
      base_slope < -0.1  ~ "Moderate divergence (B6 increases faster)",
      TRUE               ~ "Weak/parallel"
    ),
    
    detailed_pattern = paste0(trajectory_pattern, " | ", spatial_pattern)
  ) |>
  arrange(desc(abs(base_slope)))

genes_time_uniform          <- time_strain_patterns |> filter(!strain_domain_int_sig & !threeway_int_sig) |> pull(gene_id)
genes_time_domain_affected  <- time_strain_patterns |> filter(strain_domain_int_sig  & !threeway_int_sig) |> pull(gene_id)
genes_diverging_CAST_faster <- time_strain_patterns |> filter(base_slope >  0.1)                          |> pull(gene_id)
genes_diverging_B6_faster   <- time_strain_patterns |> filter(base_slope < -0.1)                          |> pull(gene_id)


# ===================================================================
# PART 7: THREE-WAY INTERACTION PATTERNS
# Genes with domain-specific temporal strain dynamics.
# Considers both initial effects (t=0) and temporal trajectories.
# ===================================================================

threeway_patterns <- strain_effects_classified |>
  filter(any_strain_sig & threeway_int_sig) |>
  filter(threeway_MAX_pass | threeway_POST_pass) |>
  select(gene_id, symbol,
         time_strain_slope_PM, threeway_MAX_slope, threeway_POST_slope,
         threeway_MAX_pass, threeway_POST_pass,
         strain_effect_PM, strain_effect_MAX, strain_effect_POST) |>
  mutate(
    slope_PM   = time_strain_slope_PM,
    slope_MAX  = time_strain_slope_PM + threeway_MAX_slope,
    slope_POST = time_strain_slope_PM + threeway_POST_slope,
    
    init_PM   = strain_effect_PM,
    init_MAX  = strain_effect_MAX,
    init_POST = strain_effect_POST,
    
    init_range  = pmax(init_PM,  init_MAX,  init_POST)  - pmin(init_PM,  init_MAX,  init_POST),
    slope_range = pmax(slope_PM, slope_MAX, slope_POST) - pmin(slope_PM, slope_MAX, slope_POST),
    
    pattern_type = case_when(
      init_range > 1   & slope_range > 0.8 ~ "Complex: Differ in both initial state & dynamics",
      init_range > 1   & slope_range <= 0.8 ~ "Primarily initial: Domains differ at t=0, similar trajectories",
      init_range <= 1  & slope_range > 0.8  ~ "Primarily dynamic: Similar at t=0, divergent trajectories",
      TRUE                                  ~ "Subtle: Small differences in both"
    ),
    
    combined_PM   = abs(init_PM   - mean(c(init_PM,   init_MAX,  init_POST))) +
      abs(slope_PM  - mean(c(slope_PM,  slope_MAX, slope_POST))),
    combined_MAX  = abs(init_MAX  - mean(c(init_PM,   init_MAX,  init_POST))) +
      abs(slope_MAX - mean(c(slope_PM,  slope_MAX, slope_POST))),
    combined_POST = abs(init_POST - mean(c(init_PM,   init_MAX,  init_POST))) +
      abs(slope_POST- mean(c(slope_PM,  slope_MAX, slope_POST))),
    
    most_distinct_domain = case_when(
      combined_PM > combined_MAX & combined_PM > combined_POST ~ "PM",
      combined_MAX > combined_POST                             ~ "MAX",
      TRUE                                                     ~ "POST"
    ),
    
    slope_direction_consistent = sign(slope_PM) == sign(slope_MAX) &
      sign(slope_MAX) == sign(slope_POST),
    
    detailed_pattern = paste0(pattern_type, " | ", most_distinct_domain, " most distinct")
  ) |>
  arrange(desc(slope_range + init_range))

genes_threeway_complex <- threeway_patterns |> filter(str_detect(pattern_type, "Complex"))           |> pull(gene_id)
genes_threeway_initial <- threeway_patterns |> filter(str_detect(pattern_type, "Primarily initial")) |> pull(gene_id)
genes_threeway_dynamic <- threeway_patterns |> filter(str_detect(pattern_type, "Primarily dynamic")) |> pull(gene_id)


# ===================================================================
# PART 8: SUMMARY TABLE
# ===================================================================

pattern_summary <- tibble(
  Category = c(
    rep("Strain Main Effect",        5),
    rep("Domain Main Effect",        4),
    rep("Time Main Effect",          2),
    rep("Domain × Time Interaction", 3),
    rep("Domain × Strain",           8),
    rep("Time × Strain",             4),
    rep("Three-Way Interaction",     3)
  ),
  Pattern = c(
    # Main strain
    "CAST >> B6 (strong)", "CAST > B6 (moderate)",
    "B6 >> CAST (strong)", "B6 > CAST (moderate)",
    "Consistent only (no interactions)",
    # Main domain
    "Increasing posterior (PM lowest)", "Decreasing posterior (PM highest)",
    "MAX peak (non-monotone)", "High domain specificity (range > 1.0)",
    # Main time
    "Expression increases over time", "Expression decreases over time",
    # Domain × Time
    "MAX slope differs from PM", "POST slope differs from PM",
    "Both MAX & POST slopes differ from PM",
    # Domain × Strain
    "PM strongest strain effect", "MAX strongest strain effect",
    "POST strongest strain effect", "Opposite directions across domains",
    "High domain specificity (range > 1.0)",
    "Static (no time interactions)", "Time-affected (main time effect)",
    "Time-varying (three-way interaction)",
    # Time × Strain
    "Diverging (CAST faster)", "Diverging (B6 faster)",
    "Temporal dynamics: Uniform across domains",
    "Temporal dynamics: Domain-affected (parallel slopes)",
    # Three-way
    "Complex (differ in initial state & dynamics)",
    "Primarily initial differences",
    "Primarily dynamic differences"
  ),
  N_Genes = c(
    length(genes_CAST_high_strong), length(genes_CAST_high_mod),
    length(genes_B6_high_strong), length(genes_B6_high_mod),
    length(genes_main_only),
    length(genes_domain_increasing_posterior), length(genes_domain_decreasing_posterior),
    length(genes_domain_MAX_peak), length(genes_domain_high_specificity),
    length(genes_time_increasing), length(genes_time_decreasing),
    length(genes_domain_time_MAX_diverges), length(genes_domain_time_POST_diverges),
    length(genes_domain_time_both_diverge),
    length(genes_PM_strongest_strain), length(genes_MAX_strongest_strain),
    length(genes_POST_strongest_strain), length(genes_opposite_directions),
    length(genes_high_domain_specificity),
    length(genes_domain_strain_static), length(genes_domain_strain_time_affected),
    nrow(domain_strain_patterns |> filter(str_detect(time_status, "Three-way"))),
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
    "Domain sig, no interactions; expression higher toward posterior",
    "Domain sig, no interactions; expression higher toward anterior",
    "Domain sig, no interactions; MAX domain shows peak expression",
    "Domain sig, no interactions; LFC range > 1.0 across domains",
    "Time sig, no interactions; expression increases over development",
    "Time sig, no interactions; expression decreases over development",
    "Domain×Time sig; MAX time slope differs from PM",
    "Domain×Time sig; POST time slope differs from PM",
    "Domain×Time sig; both MAX and POST slopes differ from PM",
    "Domain×Strain sig; PM has strongest |strain effect|",
    "Domain×Strain sig; MAX has strongest |strain effect|",
    "Domain×Strain sig; POST has strongest |strain effect|",
    "Domain×Strain sig; strain effect direction varies across domains",
    "Domain×Strain sig; strain effect range > 1.0 log2FC across domains",
    "Domain×Strain sig; no time/3-way interaction",
    "Domain×Strain sig; has time×strain interaction",
    "Domain×Strain sig; has three-way interaction (pattern changes over time)",
    "Time×Strain sig; positive slope (CAST increases faster)",
    "Time×Strain sig; negative slope (B6 increases faster)",
    "Time×Strain sig; no domain/3-way interaction",
    "Time×Strain sig; has domain interaction (parallel slopes)",
    "Three-way: Domains differ in both t=0 state and slope",
    "Three-way: Domains differ mainly at t=0, similar slopes",
    "Three-way: Similar at t=0, divergent slopes"
  )
)

print(pattern_summary)


# ===================================================================
# PART 9: GENE ASSIGNMENT TABLE
# ===================================================================

all_classified_genes <- strain_effects_classified |>
  filter(any_strain_sig | domain_sig | time_sig) |>
  select(gene_id, symbol,
         main_strain_lfc, main_strain_pval, main_strain_pattern,
         has_any_strain_interact,
         strain_domain_int_sig, strain_time_int_sig, threeway_int_sig,
         passes_domain_filter, passes_time_filter, passes_domain_time_filter)

gene_assignments <- all_classified_genes |>
  mutate(
    # Main strain
    CAST_high_strong  = gene_id %in% genes_CAST_high_strong,
    CAST_high_mod     = gene_id %in% genes_CAST_high_mod,
    B6_high_strong    = gene_id %in% genes_B6_high_strong,
    B6_high_mod       = gene_id %in% genes_B6_high_mod,
    Main_strain_only  = gene_id %in% genes_main_only,
    
    # Main domain
    Domain_increasing_posterior = gene_id %in% genes_domain_increasing_posterior,
    Domain_decreasing_posterior = gene_id %in% genes_domain_decreasing_posterior,
    Domain_MAX_peak             = gene_id %in% genes_domain_MAX_peak,
    Domain_high_specificity     = gene_id %in% genes_domain_high_specificity,
    
    # Main time
    Time_increasing = gene_id %in% genes_time_increasing,
    Time_decreasing = gene_id %in% genes_time_decreasing,
    
    # Domain × Time interaction
    DomainTime_MAX_diverges  = gene_id %in% genes_domain_time_MAX_diverges,
    DomainTime_POST_diverges = gene_id %in% genes_domain_time_POST_diverges,
    DomainTime_both_diverge  = gene_id %in% genes_domain_time_both_diverge,
    
    # Domain × Strain interaction
    PM_strongest_strain           = gene_id %in% genes_PM_strongest_strain,
    MAX_strongest_strain          = gene_id %in% genes_MAX_strongest_strain,
    POST_strongest_strain         = gene_id %in% genes_POST_strongest_strain,
    Opposite_directions           = gene_id %in% genes_opposite_directions,
    High_domain_strain_specificity = gene_id %in% genes_high_domain_specificity,
    Domain_strain_static          = gene_id %in% genes_domain_strain_static,
    Domain_strain_time_affected   = gene_id %in% genes_domain_strain_time_affected,
    
    # Time × Strain interaction
    Diverging_CAST_faster = gene_id %in% genes_diverging_CAST_faster,
    Diverging_B6_faster   = gene_id %in% genes_diverging_B6_faster,
    Time_strain_uniform         = gene_id %in% genes_time_uniform,
    Time_strain_domain_affected = gene_id %in% genes_time_domain_affected,
    
    # Three-way interaction
    ThreeWay_complex = gene_id %in% genes_threeway_complex,
    ThreeWay_initial = gene_id %in% genes_threeway_initial,
    ThreeWay_dynamic = gene_id %in% genes_threeway_dynamic
  )


# ===================================================================
# PART 10: SAVE OUTPUTS
# ===================================================================

write_csv(pattern_summary,       here("strain_pattern_summary.csv"))
write_csv(main_patterns,         here("main_strain_effect_patterns.csv"))
write_csv(domain_main_patterns,  here("domain_main_effect_patterns.csv"))
write_csv(time_main_patterns,    here("time_main_effect_patterns.csv"))
write_csv(domain_time_patterns,  here("domain_time_interaction_patterns.csv"))
write_csv(domain_strain_patterns,here("domain_strain_patterns.csv"))
write_csv(time_strain_patterns,  here("time_strain_patterns.csv"))
write_csv(threeway_patterns,     here("threeway_patterns.csv"))
write_csv(gene_assignments,      here("strain_gene_pattern_assignments.csv"))
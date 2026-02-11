# ===================================================================
# STRAIN GENE PATTERN CLASSIFICATION
# Uses objects created in Quarto document: strain_effects_classified
# Run this AFTER the DESeq2 analysis chunks in the Quarto document
# ===================================================================

library(tidyverse)
library(here)

# ===================================================================
# PART 1: MAIN STRAIN EFFECT PATTERNS
# These are genes with a consistent strain effect across domains and
# time — no significant domain × strain or time × strain interactions.
# Uses strain_effects_filtered (which excludes interaction genes).
# ===================================================================

main_patterns <- strain_effects_filtered |>
  select(gene_id, symbol, main_strain_lfc, main_strain_pval, main_strain_pattern) |>
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

# ===================================================================
# PART 2: DOMAIN × STRAIN INTERACTION PATTERNS
# These genes have strain effects that differ across palate domains.
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
         max_domain_effect, max_domain_direction, domain_consistency) |>
  mutate(
    detailed_pattern = paste(max_domain_direction, "|", domain_consistency)
  ) |>
  arrange(desc(abs(domain_int_MAX_lfc) + abs(domain_int_POST_lfc)))

genes_PM_biggest_effect <- domain_patterns |>
  filter(max_domain_effect == "PM") |> pull(gene_id)
genes_MAX_biggest_effect <- domain_patterns |>
  filter(max_domain_effect == "MAX") |> pull(gene_id)
genes_POST_biggest_effect <- domain_patterns |>
  filter(max_domain_effect == "POST") |> pull(gene_id)
genes_opposite_directions <- domain_patterns |>
  filter(domain_consistency == "Variable") |> pull(gene_id)

# ===================================================================
# PART 3: TIME × STRAIN INTERACTION PATTERNS
# These genes have strain effects that change over developmental time.
# Source: strain_effects_classified with significant time LRT.
# ===================================================================

time_patterns <- strain_effects_classified |>
  filter(any_strain_sig & time_int_sig & time_int_pass) |>
  select(gene_id, symbol, time_int_slope, time_int_slope_lfcSE, time_int_slope_pval, time_int_pval,
         main_strain_lfc, time_pattern) |>
  mutate(
    trajectory_pattern = case_when(
      time_int_slope > 0.2 ~ "Strong divergence (CAST faster)",
      time_int_slope > 0.1 ~ "Moderate divergence (CAST faster)",
      time_int_slope < -0.2 ~ "Strong divergence (B6 faster)",
      time_int_slope < -0.1 ~ "Moderate divergence (B6 faster)",
      TRUE ~ "Weak/parallel"
    )
  ) |>
  arrange(desc(abs(time_int_slope)))

genes_diverging_CAST_faster <- time_patterns |>
  filter(time_int_slope > 0.1) |> pull(gene_id)
genes_diverging_B6_faster <- time_patterns |>
  filter(time_int_slope < -0.1) |> pull(gene_id)

# ===================================================================
# PART 4: THREE-WAY INTERACTION PATTERNS
# These genes have domain-specific temporal strain dynamics.
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
    slope_PM = time_int_slope,
    slope_MAX = time_int_slope + threeway_MAX_slope,
    slope_POST = time_int_slope + threeway_POST_slope,
    
    max_slope_diff_domain = case_when(
      abs(slope_PM - mean(c(slope_MAX, slope_POST))) > 
        abs(slope_MAX - mean(c(slope_PM, slope_POST))) &
        abs(slope_PM - mean(c(slope_MAX, slope_POST))) > 
        abs(slope_POST - mean(c(slope_PM, slope_MAX))) ~ "PM",
      abs(slope_MAX - mean(c(slope_PM, slope_POST))) > 
        abs(slope_POST - mean(c(slope_PM, slope_MAX))) ~ "MAX",
      TRUE ~ "POST"
    ),
    
    slope_consistency = case_when(
      sign(slope_PM) == sign(slope_MAX) & sign(slope_MAX) == sign(slope_POST) ~ "Consistent",
      TRUE ~ "Variable"
    ),
    
    pattern = paste(max_slope_diff_domain, "diverges most |", slope_consistency)
  ) |>
  arrange(desc(abs(threeway_MAX_slope) + abs(threeway_POST_slope)))

genes_PM_biggest_slope_diff <- threeway_patterns |>
  filter(max_slope_diff_domain == "PM") |> pull(gene_id)
genes_MAX_biggest_slope_diff <- threeway_patterns |>
  filter(max_slope_diff_domain == "MAX") |> pull(gene_id)
genes_POST_biggest_slope_diff <- threeway_patterns |>
  filter(max_slope_diff_domain == "POST") |> pull(gene_id)
genes_consistent_slopes <- threeway_patterns |>
  filter(slope_consistency == "Consistent") |> pull(gene_id)

# ===================================================================
# PART 5: DOMAIN SPECIFICITY PATTERNS
# Which domains show significant AND well-estimated strain effects?
# Source: strain_effects_classified with significant domain LRT.
# ===================================================================

domain_specificity <- strain_effects_classified |>
  filter(any_strain_sig & domain_int_sig) |>
  mutate(
    effect_PM = strain_effect_PM_pass,
    effect_MAX = strain_effect_MAX_pass,
    effect_POST = strain_effect_POST_pass,
    
    specificity_pattern = case_when(
      effect_PM & !effect_MAX & !effect_POST ~ "PM only",
      !effect_PM & effect_MAX & !effect_POST ~ "MAX only",
      !effect_PM & !effect_MAX & effect_POST ~ "POST only",
      effect_PM & effect_MAX & !effect_POST ~ "PM & MAX",
      effect_PM & !effect_MAX & effect_POST ~ "PM & POST",
      !effect_PM & effect_MAX & effect_POST ~ "MAX & POST",
      effect_PM & effect_MAX & effect_POST ~ "All domains",
      TRUE ~ "None (filtered)"
    )
  ) |>
  select(gene_id, symbol, 
         strain_effect_PM, strain_effect_MAX, strain_effect_POST,
         effect_PM, effect_MAX, effect_POST, specificity_pattern)

genes_strain_PM_only <- domain_specificity |>
  filter(specificity_pattern == "PM only") |> pull(gene_id)
genes_strain_MAX_only <- domain_specificity |>
  filter(specificity_pattern == "MAX only") |> pull(gene_id)
genes_strain_POST_only <- domain_specificity |>
  filter(specificity_pattern == "POST only") |> pull(gene_id)
genes_strain_all_domains <- domain_specificity |>
  filter(specificity_pattern == "All domains") |> pull(gene_id)

# ===================================================================
# PART 6: SUMMARY TABLE
# ===================================================================

pattern_summary <- tibble(
  Category = c(
    rep("Main Strain Effect", 4),
    rep("Domain × Strain", 4),
    rep("Time × Strain", 2),
    rep("Three-Way Interaction", 4),
    rep("Domain Specificity", 4)
  ),
  Pattern = c(
    "CAST >> B6 (strong)", "CAST > B6 (moderate)",
    "B6 >> CAST (strong)", "B6 > CAST (moderate)",
    "Biggest effect in PM", "Biggest effect in MAX", 
    "Biggest effect in POST", "Opposite directions",
    "Diverging (CAST faster)", "Diverging (B6 faster)",
    "PM diverges most", "MAX diverges most",
    "POST diverges most", "Consistent slopes",
    "PM only", "MAX only", "POST only", "All domains"
  ),
  N_Genes = c(
    length(genes_CAST_high_strong), length(genes_CAST_high_mod),
    length(genes_B6_high_strong), length(genes_B6_high_mod),
    length(genes_PM_biggest_effect), length(genes_MAX_biggest_effect),
    length(genes_POST_biggest_effect), length(genes_opposite_directions),
    length(genes_diverging_CAST_faster), length(genes_diverging_B6_faster),
    length(genes_PM_biggest_slope_diff), length(genes_MAX_biggest_slope_diff),
    length(genes_POST_biggest_slope_diff), length(genes_consistent_slopes),
    length(genes_strain_PM_only), length(genes_strain_MAX_only),
    length(genes_strain_POST_only), length(genes_strain_all_domains)
  ),
  Description = c(
    "Consistent strain effect, LFC > 1", "Consistent strain effect, 0.5 < LFC ≤ 1",
    "Consistent strain effect, LFC < -1", "Consistent strain effect, -1 ≤ LFC < -0.5",
    "Domain int. sig, largest |effect| in PM", "Domain int. sig, largest |effect| in MAX",
    "Domain int. sig, largest |effect| in POST", "Domain int. sig, effect direction varies",
    "Time int. sig, positive slope", "Time int. sig, negative slope",
    "3-way sig, PM most different slope", "3-way sig, MAX most different slope",
    "3-way sig, POST most different slope", "3-way sig, all slopes same direction",
    "Domain int. sig, strain effect only in PM", "Domain int. sig, strain effect only in MAX",
    "Domain int. sig, strain effect only in POST", "Domain int. sig, strain effect in all"
  )
)

print(pattern_summary)

# ===================================================================
# PART 7: GENE ASSIGNMENT TABLE
# One row per gene with significant strain effect of any kind.
# ===================================================================

# All genes with any significant strain effect
all_strain_genes <- strain_effects_classified |>
  filter(any_strain_sig) |>
  select(gene_id, symbol)

gene_assignments <- all_strain_genes |>
  mutate(
    # Main effects (consistent, no interactions)
    CAST_high_strong = gene_id %in% genes_CAST_high_strong,
    CAST_high_mod = gene_id %in% genes_CAST_high_mod,
    B6_high_strong = gene_id %in% genes_B6_high_strong,
    B6_high_mod = gene_id %in% genes_B6_high_mod,
    # Domain interaction patterns
    PM_biggest_effect = gene_id %in% genes_PM_biggest_effect,
    MAX_biggest_effect = gene_id %in% genes_MAX_biggest_effect,
    POST_biggest_effect = gene_id %in% genes_POST_biggest_effect,
    Opposite_directions = gene_id %in% genes_opposite_directions,
    # Time interaction patterns
    Diverging_CAST_faster = gene_id %in% genes_diverging_CAST_faster,
    Diverging_B6_faster = gene_id %in% genes_diverging_B6_faster,
    # Three-way patterns
    ThreeWay_PM = gene_id %in% genes_PM_biggest_slope_diff,
    ThreeWay_MAX = gene_id %in% genes_MAX_biggest_slope_diff,
    ThreeWay_POST = gene_id %in% genes_POST_biggest_slope_diff,
    ThreeWay_consistent = gene_id %in% genes_consistent_slopes,
    # Domain specificity
    Strain_PM_only = gene_id %in% genes_strain_PM_only,
    Strain_MAX_only = gene_id %in% genes_strain_MAX_only,
    Strain_POST_only = gene_id %in% genes_strain_POST_only,
    Strain_all_domains = gene_id %in% genes_strain_all_domains
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
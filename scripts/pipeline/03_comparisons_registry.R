# 03_comparisons_registry.R — define the 42 main and 30 sensitivity comparisons
# Usage: Rscript scripts/pipeline/03_comparisons_registry.R
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", ".")   # default: the current directory
#                                                 # on; set PBL_ROOT to run elsewhere.
out_d <- file.path(root, "data/02_processed"); rep_d <- file.path(root, "doc/05_reports")

map <- read.delim(file.path(root, "data/01_metadata/sample_mapping.tsv"), stringsAsFactors = FALSE)
gens <- c("Wildtype", "HSF1.KD", "MSN24.KO", "Double.KDKO")
conds <- c("25_0", "37_10", "37_30", "42_10", "42_30")
stopifnot(setequal(unique(map$knockout), gens))

rows <- list()
add <- function(comp_id, family, g, g_ref, t, t_ref, s, s_ref, mutant_vs_wt_set = FALSE, pairs_3cond_set = FALSE, note = "")
  rows[[length(rows) + 1]] <<- data.frame(comp_id, family, genotype = g, genotype_ref = g_ref,
    temperature = t, temperature_ref = t_ref, time = s, time_ref = s_ref,
    n_rep_per_arm = 3, main42_candidate_a = mutant_vs_wt_set, main42_candidate_b = pairs_3cond_set, note)

# 1) time-series: 4 genotypes x 2 temps x 2 times vs same-genotype 25_0  -> 16 (main)
for (g in gens) for (t in c(37, 42)) for (s in c(10, 30))
  add(sprintf("%s_%d_%d_vs_25_0", g, t, s), "time_series", g, g, t, 25, s, 0,
      note = "main analysis, 42 comparisons")

# 2) genotype: all 6 pairs x 5 conditions -> 30 (sensitivity; two definitions of the 18-comparison main set are evaluated,
#    so both definitions are recorded)
# direction convention: mutant-vs-WT keeps WT as reference; mutant-mutant pairs alphabetical
pairs <- combn(gens, 2, simplify = FALSE)
for (cond in conds) for (p in pairs) {
  if ("Wildtype" %in% p) { a <- p[p != "Wildtype"]; b <- "Wildtype" } else { a <- p[1]; b <- p[2] }
  cond_t <- as.integer(sub("_.*", "", cond)); cond_s <- as.integer(sub(".*_", "", cond))
  mutant_vs_wt_set <- (b == "Wildtype" && cond != "25_0")                      # 3 mutants x 4 stress = 12
  pairs_3cond_set <- (cond %in% c("42_10", "42_30", "37_30"))                 # 6 pairs x 3 conditions = 18
  add(sprintf("%s_vs_%s_at_%s", a, b, cond), "genotype", a, b, cond_t, cond_t, cond_s, cond_s,
      mutant_vs_wt_set, pairs_3cond_set,
      note = if (mutant_vs_wt_set) "candidate: mutants-vs-WT interpretation" else if (pairs_3cond_set) "candidate: 6pairs-x-3cond interpretation" else "")
}

# 3) temperature: 42 vs 37 at same genotype & time -> 8 (main)
for (g in gens) for (s in c(10, 30))
  add(sprintf("%s_42_%d_vs_37_%d", g, s, s), "temperature", g, g, 42, 37, s, s,
      note = "main analysis, 42 comparisons")

reg <- do.call(rbind, rows)
reg <- reg[order(reg$family, reg$comp_id), ]
fwrite(reg, file.path(out_d, "comparisons.csv"))

L <- c("# Comparison registry (scripts/pipeline/03_comparisons_registry.R)", "",
       sprintf("- registered: %d comparisons (main analysis and sensitivity; all executed comparisons are defined in this registry)", nrow(reg)),
       sprintf("- time_series: 16 comparisons (main analysis)"), sprintf("- temperature: 8 comparisons (main analysis)"),
       "- genotype: 30 all-pairs comparisons as a sensitivity set, with two candidate subsets marked:",
       "  - `main42_candidate_a`: 3 mutants vs wild type x 4 stress conditions = 12 comparisons",
       "  - `main42_candidate_b`: 6 pairs x {42_10, 42_30, 37_30} = 18 comparisons",
       "- the sensitivity analysis covers all 30 genotype comparisons.",
       "- direction: in comp_id `A_vs_B`, B is the reference; within the genotype family test and reference follow alphabetical order.",
       "- prefilter: `rowSums(count >= 10) >= x`, x set in the batch run.")
writeLines(L, file.path(rep_d, "comparisons_registry_2026-09-04.md"))
cat(sprintf("registered %d comparisons -> data/02_processed/comparisons.csv\n", nrow(reg)))
cat(sprintf("families: %s\n", paste(names(table(reg$family)), table(reg$family), sep = "=", collapse = ", ")))
cat(sprintf("main-42 definition 1 sum=%d, definition 2 sum=%d\n", sum(reg$main42_candidate_a), sum(reg$main42_candidate_b)))

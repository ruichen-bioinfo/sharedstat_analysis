# p1_orientation.R -- re-expresses the S. paradoxus x S. kudriavzevii (SpxSk) fits with the first-named allele
# as reference: one reference arm per hybrid in every P1 quantity, as stated in the manuscript.
#
# WHY. 02_fay_multimethod.R sets design$Allele <- factor(design$Allele), so the reference level of every fit is
# the alphabetically first allele. That is the first-named allele in ScxSk, ScxSp, ScxSu, SkxSu and SpxSu but not in
# SpxSk, where the fits use Sk: eff[time] is Sk's 15 min response and eff[inter] is Time15:AlleleSp = Sp - Sk.
# P1's h and kappa (94, 102d, 109: substr(g, 1, 2)) use Sp. This helper re-expresses the SpxSk fits with Sp as the
# reference, exactly, so every P1 quantity uses one arm:
#   model A: time' = time + inter (Sp's 15 min response), inter' = -inter (Sk - Sp), allele' = -allele
#   model B: inter' = -inter, allele' = -allele (no time term under model B)
# stat, pval and padj are returned as stored, not re-expressed (padj is not invariant to the reference level for every
# test; see below). Other groups are returned untouched. The alphabetical crossspecies scripts (08, 09, 10, 11 default,
# 12_artefact_corrected_modifier, 14, 56) do not use this helper.
#
# Not used by the P1 analysis scripts. The effects re-express exactly, but padj does not (limma's voom-weighted
# multi-coefficient F is not invariant to the reparametrised `time` test), so P1 reads the refit
# fay_multimethod_p1order.tsv.gz written by 02b_p1order_refit.R. Kept only as the documented check that refit effects
# equal the re-expressed ones (check b of 02b_p1order_refit.R).
# Usage: source("scripts/crossspecies/p1_orientation.R"); X <- p1_orient(fread(... fay_multimethod.tsv.gz ...))
p1_orient <- function(X, focal = "SpxSk") {
  stopifnot(data.table::is.data.table(X),
            all(c("gene", "term", "eff", "condition", "group", "model", "method") %in% names(X)))
  X <- data.table::copy(X)
  by <- c("gene", "condition", "group", "model", "method")
  I <- X[group == focal & model == "A" & term == "inter", c(by, "eff"), with = FALSE]
  stopifnot(!anyDuplicated(I[, by, with = FALSE]))
  ## every model-A time row of the focal group must have its inter partner, or time' would be left unconverted
  stopifnot(X[group == focal & model == "A" & term == "time", .N] == nrow(I))
  data.table::setnames(I, "eff", "eff_inter")
  I[, term := "time"]
  n_hit <- nrow(X[I, on = c(by, "term"), nomatch = NULL])
  stopifnot(n_hit == nrow(I))
  X[I, on = c(by, "term"), eff := eff + i.eff_inter]   # model A: Sp's response = Sk's + (Sp - Sk), original inter used
  X[group == focal & term %in% c("inter", "allele"), eff := -eff]
  X[]
}

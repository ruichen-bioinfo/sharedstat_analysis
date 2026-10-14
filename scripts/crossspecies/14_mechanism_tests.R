# 14_mechanism_tests.R — three tests that the mechanistic explanation has to pass.
#
# Starting point. The current account is that HSR regulatory divergence has two
# separable components: a coherent global gain on the whole heat-induced program (invisible
# to per-gene FDR; ScxSp is the pure case) and gene-specific rewiring (what Fay counted;
# SpxSu is the pure case), and that the two are orthogonal across crosses (Spearman -0.600).
# The thermal-niche account is falsified on three instruments. That is a description. To be
# a mechanism it has to say where the gain acts, and that claim has to be falsifiable.
#
# The claim: a coherent multiplicative rescaling of an entire regulon is what divergence at a
# shared upstream node looks like -- the activator's abundance, activation threshold, or the
# chaperone-titration feedback that sets its activity -- as opposed to divergence in the
# individual promoters, which can only produce gene-specific changes. Three consequences:
#
#   (T1) Completeness of the species decomposition. The four strains are a 2x2: thermal class
#        (Sc,Sp thermophilic; Sk,Su cryophilic) crossed with provenance (Sc,Su domesticated
#        wine; Sp,Sk wild oak), one strain per cell. So the three residual df of the centred
#        additive model are exactly three orthogonal contrasts -- thermal, provenance, and
#        their interaction -- and nothing else. Script 12 section 3 tests only the first of the
#        three and finds it significantly in the wrong direction. The other two need no new
#        model, just a reparameterisation of a fit that already exists. If provenance or the interaction carries more signal than
#        thermal class does, that is a substantive redirection, not a footnote.
#
#   (T2) The gain should be regulon-specific. If the gain sits on the Hsf1 node, the gain
#        slope estimated within Hsf1-dependent genes should exceed the slope within
#        Msn2/4-dependent genes. If the gain is instead something global and non-regulatory
#        (RNA composition, normalisation, growth rate), it should be the same in both. Note
#        this is not the test of scripts 06-07: that one asks whether the per-gene
#        divergence rate differs by module and finds no difference -- which is the prediction
#        for the rewiring component. The two together are a double dissociation, and the "no"
#        in scripts 06-07 is what makes a "yes" here meaningful rather than a generic module effect.
#
#   (T3) The gain should actually be multiplicative. "Global multiplicative rescaling" is a
#        strong functional form: eff_inter proportional to the shared response, one constant
#        for all genes. It is falsified if the relation bends -- if it is driven by the
#        strongly-induced tail, or saturates. Tested with a grouping estimator: bin genes on
#        the independent magnitude S (other crosses, co-pooled partner excluded), then take
#        per-bin means. Selection on S is independent of both noises, so per-bin means of the
#        focal cross are unbiased and the artefact of script 11 does not enter -- the
#        binned slope is a consistent estimator in its own right, so it also cross-checks the
#        IV slope by an entirely different route.
#
# Usage: Rscript scripts/crossspecies/14_mechanism_tests.R [--nperm=200]

suppressPackageStartupMessages({ library(data.table) })

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
outd <- file.path(root, "data/02_processed/crossspecies")
HYB <- c("ScxSk","ScxSp","ScxSu","SkxSu","SpxSk","SpxSu")
SPECIES <- c("Sc","Sk","Sp","Su")
METHODS <- c("DESeq2","edgeR","limma")
ALPHA <- 0.05
partner <- c(ScxSp = "SkxSu", SkxSu = "ScxSp", ScxSu = "SpxSk",
             SpxSk = "ScxSu", ScxSk = "SpxSu", SpxSu = "ScxSk")
args <- commandArgs(trailingOnly = TRUE)
np <- grep("^--nperm=", args, value = TRUE)
NPERM <- if (length(np)) as.integer(sub("^--nperm=", "", np[1])) else 200
say <- function(...) cat(sprintf(...), "\n", sep = "")
set.seed(20260905)
say("=== 14_mechanism_tests | %s | nperm=%d ===", format(Sys.time()), NPERM)

X <- fread(cmd = paste("zcat", shQuote(file.path(outd, "fay_multimethod.tsv.gz"))))
X <- X[condition == "heat" & group %in% HYB & model == "A"]
GIV <- fread(file.path(outd, "gain_iv.tsv"))
Rtab <- X[term == "time", .(gene, group, method, R = eff, padj_t = padj)]
IItab <- X[term == "inter", .(gene, group, method, ei = eff, padj_i = padj)]

## the 2x2: thermal class x provenance, one strain per cell (Fay et al. 2023 strain table)
CLS <- data.table(sp = SPECIES,
                  thermal  = c("thermo","cryo","thermo","cryo"),
                  source   = c("dom","wild","wild","dom"),
                  strain   = c("T73","CR85","N17","BMV58"))
say("\nthe four strains form a complete 2x2 with one strain per cell:")
say("%-6s %-8s %-6s %-6s", "sp", "strain", "thermal", "source")
for (i in 1:4) say("%-6s %-8s %-6s %-6s", CLS$sp[i], CLS$strain[i], CLS$thermal[i], CLS$source[i])

## three mutually orthogonal contrasts on (Sc, Sk, Sp, Su), spanning the 3 residual df
CON <- rbind(thermal = c(-1, 1, -1, 1)/2,       # cryo  - thermo
             source  = c( 1,-1, -1, 1)/2,       # dom   - wild
             interact= c( 1, 1, -1,-1)/2)       # (dom-wild in thermo) - (dom-wild in cryo)
colnames(CON) <- SPECIES
say("\northogonality check (CON %%*%% t(CON), should be diagonal):")
print(round(CON %*% t(CON), 6))
say("this is a complete reparameterisation: 3 orthogonal contrasts, 3 residual df, nothing left over")

## ------------------------------------------------- shared machinery from script 12 ----
Dm <- matrix(0, 6, 4, dimnames = list(HYB, SPECIES))
for (g in HYB) { a <- sort(unlist(strsplit(g, "x"))); Dm[g, a[1]] <- -1; Dm[g, a[2]] <- 1 }
Dc <- Dm %*% (diag(4) - matrix(1/4, 4, 4))
Gi <- MASS::ginv(Dc)

build <- function(mt) {
  W  <- dcast(IItab[method == mt, .(gene, group, ei)], gene ~ group, value.var = "ei")
  Rw <- dcast(Rtab[method == mt], gene ~ group, value.var = "R")
  M  <- merge(W, Rw, by = "gene", suffixes = c("", ".R"))
  Y  <- as.matrix(M[, ..HYB]); RR <- as.matrix(M[, paste0(HYB, ".R"), with = FALSE])
  keep <- complete.cases(Y) & complete.cases(RR) & rowSums(Y^2) > 0
  a <- GIV[model == "A" & method == mt, .(artefact = median(artefact)), by = group]
  lam <- a$artefact[match(HYB, a$group)]
  list(gene = M$gene[keep], R = RR[keep, , drop = FALSE],
       corr = Y[keep, , drop = FALSE] - RR[keep, , drop = FALSE] * rep(lam, each = sum(keep)))
}

## ============================ (T1) the full 2x2 decomposition ==========================
say("\n=== (T1) complete 2x2 decomposition of the corrected species modifiers ===")
say("each contrast's per-gene value, median over genes with R2 > 0.8. Null = shuffle the six")
say("cross labels within gene, refit, recompute -- the same null as script 12 section 1.")
say("%-7s %-9s %7s %9s %8s %11s %11s", "method","contrast","n","median","frac>0","null_95pct","p_perm")
resT1 <- list()
r2_of <- function(Y) { H <- Dc %*% Gi; rss <- rowSums((Y - Y %*% t(H))^2)
                       tss <- rowSums(Y^2); ifelse(tss > 0, 1 - rss/tss, NA_real_) }
for (mt in METHODS) {
  B <- build(mt); Y <- B$corr
  r2 <- r2_of(Y); sel <- which(!is.na(r2) & r2 > 0.8)
  Mh <- Y %*% t(Gi); colnames(Mh) <- SPECIES
  CV <- Mh[sel, , drop = FALSE] %*% t(CON)            # genes x 3 contrasts
  ## null: permute cross labels within gene, refit, recompute contrasts on the same genes
  nullabs <- matrix(NA_real_, NPERM, 3, dimnames = list(NULL, rownames(CON)))
  for (k in seq_len(NPERM)) {
    Yp <- Y[, sample(6), drop = FALSE]
    Mp <- (Yp %*% t(Gi)); colnames(Mp) <- SPECIES
    nullabs[k, ] <- apply(abs(Mp[sel, , drop = FALSE] %*% t(CON)), 2, median)
  }
  for (j in 1:3) {
    v <- CV[, j]; obs <- median(abs(v))
    pp <- (1 + sum(nullabs[, j] >= obs)) / (1 + NPERM)
    bt <- binom.test(sum(v > 0), length(v), 0.5)
    resT1[[length(resT1)+1]] <- data.table(method = mt, contrast = rownames(CON)[j],
      n = length(v), median = median(v), abs_median = obs, frac_pos = mean(v > 0),
      null_95 = quantile(nullabs[, j], 0.95), p_perm = pp, p_sign = bt$p.value)
    say("%-7s %-9s %7d %+9.4f %7.1f%% %11.4f %11.4g", mt, rownames(CON)[j], length(v),
        median(v), 100*mean(v > 0), quantile(nullabs[, j], 0.95), pp)
  }
}
T1 <- rbindlist(resT1)
fwrite(T1, file.path(outd, "species_2x2_decomposition.tsv"), sep = "\t")
say("\nranking of the three contrasts by |median| (median over methods):")
rk <- T1[, .(abs_med = median(abs_median), md = median(median)), by = contrast][order(-abs_med)]
for (i in seq_len(nrow(rk))) say("  %-9s |%.4f|  (signed %+.4f)", rk$contrast[i], rk$abs_med[i], rk$md[i])

## ============================ (T1b) same three contrasts on the IV gain slopes =========
say("\n=== (T1b) the same three contrasts on the six IV gain slopes (exact permutation) ===")
say("6 slopes, 4 species scalars, 3 residual df. All 720 relabellings enumerated, so the")
say("p-values are exact rather than Monte Carlo.")
perms <- do.call(rbind, lapply(1:720, function(i) NULL))
allperm <- function(v) { n <- length(v)
  if (n == 1) return(matrix(v, 1)); out <- NULL
  for (i in seq_len(n)) { sub <- allperm(v[-i])
    out <- rbind(out, cbind(v[i], sub)) }; out }
PM <- allperm(1:6)
say("%-7s %-9s %10s %12s %10s", "method","contrast","observed","null_95_abs","p_exact")
resT1b <- list()
for (mt in METHODS) {
  sl <- GIV[model == "A" & method == mt, .(iv = median(slope_iv)), by = group]
  v <- sl$iv[match(HYB, sl$group)]
  obs <- as.numeric(CON %*% Gi %*% v)
  nl <- t(apply(PM, 1, function(p) as.numeric(CON %*% Gi %*% v[p])))
  for (j in 1:3) {
    pe <- mean(abs(nl[, j]) >= abs(obs[j]))
    resT1b[[length(resT1b)+1]] <- data.table(method = mt, contrast = rownames(CON)[j],
      observed = obs[j], null95 = quantile(abs(nl[, j]), 0.95), p_exact = pe)
    say("%-7s %-9s %+10.4f %12.4f %10.4g", mt, rownames(CON)[j], obs[j],
        quantile(abs(nl[, j]), 0.95), pe)
  }
}
fwrite(rbindlist(resT1b), file.path(outd, "gain_2x2_decomposition.tsv"), sep = "\t")

## ============================ (T2) is the gain regulon-specific? ========================
say("\n=== (T2) gain slope within the Hsf1 regulon vs within the Msn2/4 regulon ===")
MOD  <- fread(file.path(outd, "yeast_modules_37_30.tsv"))
HOGS <- fread(file.path(outd, "hog_to_scer.tsv"))
mkey <- function(ym) merge(HOGS[!is.na(sys) & in_yeast == TRUE, .(gene, sys)],
                           MOD[method == ym, .(sys = gene, module)], by = "sys")[, .(gene, module)]
key <- mkey("DESeq2")
key[, reg := fifelse(module %in% c("Hsf1","both"), "Hsf1",
             fifelse(module == "Msn24", "Msn24", NA_character_))]
key <- key[!is.na(reg)]
say("regulon sizes after HOG bridging: Hsf1(+both) %d, Msn2/4 %d",
    key[reg == "Hsf1", .N], key[reg == "Msn24", .N])

## the instrument: shared response over the other crosses, co-pooled partner excluded
instrument <- function(g, mt) {
  keep <- setdiff(HYB, c(g, partner[[g]]))
  Rtab[group %in% keep & method == mt,
       .(S = mean(R, na.rm = TRUE), nS = sum(is.finite(R)),
         nsig = sum(!is.na(padj_t) & padj_t < ALPHA)), by = gene][nS >= 3]
}
say("\nper cross (n>=40 per regulon required):")
say("%-8s %-7s %-7s %6s %10s %10s", "cross","method","regulon","n","slope_iv","slope_ols")
resT2 <- list()
for (g in HYB) for (mt in METHODS) {
  D <- merge(merge(IItab[group == g & method == mt, .(gene, ei)],
                   Rtab[group == g & method == mt, .(gene, R, padj_t)], by = "gene"),
             instrument(g, mt), by = "gene")
  D <- merge(D, key, by = "gene")
  D <- D[!is.na(padj_t) & padj_t < ALPHA & is.finite(ei) & is.finite(R) & is.finite(S)]
  for (rg in c("Hsf1","Msn24")) {
    s <- D[reg == rg]
    if (nrow(s) < 40) next
    iv <- cov(s$ei, s$S) / cov(s$R, s$S)
    ols <- unname(coef(lm(ei ~ R, s))[2])
    resT2[[length(resT2)+1]] <- data.table(group = g, method = mt, regulon = rg,
                                           n = nrow(s), slope_iv = iv, slope_ols = ols)
    say("%-8s %-7s %-7s %6d %+10.4f %+10.4f", g, mt, rg, nrow(s), iv, ols)
  }
}
T2 <- rbindlist(resT2)
fwrite(T2, file.path(outd, "gain_by_regulon.tsv"), sep = "\t")
if (nrow(T2)) {
  W <- dcast(T2, group + method ~ regulon, value.var = "slope_iv")
  W <- W[is.finite(Hsf1) & is.finite(Msn24)]
  say("\npaired within cross x method, n=%d cells:", nrow(W))
  say("  Hsf1 median   %+.4f", median(W$Hsf1))
  say("  Msn2/4 median %+.4f", median(W$Msn24))
  say("  paired difference (Hsf1 - Msn2/4) median %+.4f", median(W$Hsf1 - W$Msn24))
  say("  cells with Hsf1 gain steeper (more negative): %d of %d", sum(W$Hsf1 < W$Msn24), nrow(W))
  if (nrow(W) >= 4) {
    wt <- suppressWarnings(wilcox.test(W$Hsf1, W$Msn24, paired = TRUE))
    bt <- binom.test(sum(W$Hsf1 < W$Msn24), nrow(W), 0.5)
    say("  paired Wilcoxon p = %.4g | sign test p = %.4g", wt$p.value, bt$p.value)
  }
  say("\nprediction was: Hsf1 steeper if the gain sits on the Hsf1 node; equal if the gain is")
  say("a non-regulatory global factor. Section 6.5.3 found the per-gene divergence RATE does")
  say("not differ by module, so a difference here would be a double dissociation.")
}

## ============================ (T3) is the gain actually multiplicative? ================
say("\n=== (T3) functional form: grouping estimator on independent magnitude bins ===")
say("bins are defined on S (other crosses only), so per-bin means of the focal cross are")
say("unbiased. A multiplicative gain predicts per-bin mean(ei) proportional to mean(R) with")
say("one constant; curvature falsifies it.")
NB <- 10
say("\n%-8s %-7s %6s %10s %10s %10s %10s", "cross","method","bins","slope_bin","slope_iv",
    "R2_linear","p_curve")
resT3 <- list(); binlong <- list()
for (g in HYB) for (mt in METHODS) {
  D <- merge(merge(IItab[group == g & method == mt, .(gene, ei)],
                   Rtab[group == g & method == mt, .(gene, R, padj_t)], by = "gene"),
             instrument(g, mt)[nsig >= 2], by = "gene")
  D <- D[is.finite(ei) & is.finite(R) & is.finite(S)]
  if (nrow(D) < 500) next
  D[, bin := cut(S, breaks = quantile(S, seq(0, 1, length.out = NB + 1), na.rm = TRUE),
                 include.lowest = TRUE, labels = FALSE)]
  bt <- D[!is.na(bin), .(mR = mean(R), mE = mean(ei), n = .N), by = bin][order(bin)][n >= 30]
  if (nrow(bt) < 6) next
  lin <- lm(mE ~ 0 + mR, bt, weights = bt$n)
  qd  <- lm(mE ~ 0 + mR + I(mR^2), bt, weights = bt$n)
  pc  <- anova(lin, qd)$`Pr(>F)`[2]
  ivv <- GIV[model == "A" & method == mt & group == g, median(slope_iv)]
  resT3[[length(resT3)+1]] <- data.table(group = g, method = mt, bins = nrow(bt),
    slope_bin = unname(coef(lin)[1]), slope_iv = ivv,
    r2_lin = summary(lin)$r.squared, p_curve = pc)
  binlong[[length(binlong)+1]] <- data.table(group = g, method = mt, bt)
  say("%-8s %-7s %6d %+10.4f %+10.4f %10.3f %10.4g", g, mt, nrow(bt),
      coef(lin)[1], ivv, summary(lin)$r.squared, pc)
}
T3 <- rbindlist(resT3)
fwrite(T3, file.path(outd, "gain_functional_form.tsv"), sep = "\t")
fwrite(rbindlist(binlong), file.path(outd, "gain_binned_means.tsv"), sep = "\t")
if (nrow(T3)) {
  say("\nlinear-through-origin fit: median R2 %.3f over %d cells", median(T3$r2_lin), nrow(T3))
  say("curvature term significant at 0.05 in %d of %d cells", sum(T3$p_curve < 0.05, na.rm = TRUE), nrow(T3))
  say("agreement of the two independent estimators (binned vs IV): Spearman %.3f, median |diff| %.4f",
      cor(T3$slope_bin, T3$slope_iv, method = "spearman"), median(abs(T3$slope_bin - T3$slope_iv)))
}

say("\nwrote species_2x2_decomposition.tsv, gain_2x2_decomposition.tsv, gain_by_regulon.tsv,")
say("      gain_functional_form.tsv, gain_binned_means.tsv")
say("=== done ===")

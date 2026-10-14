## 87_allele_rho.R -- helper sourced by 87a/87b. Inter-gene correlation on the data allele-contrast quantities
## are built from: per hybrid, per-library allele difference
## (a2 - a1 in the same library), 8 libraries (4 timepoints x 2 blocks), residuals from a timepoint-means fit (4 residual df); interGeneCorrelation for
## rho-bar, and mean pairwise r^2 - 1/4 for the bias-corrected mean rho^2; each averaged over the six hybrids.
## Same normalisation as 04/82 (per-hybrid TMM, logCPM, prior.count 2). Library-level shifts cancel in the difference.
suppressPackageStartupMessages({ library(limma); library(edgeR) })
.AR <- local({
  root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>"); ext <- file.path(root, "data/03_external/fay2023")
  dat <- read.csv(file.path(ext, "htseq.ortho.csv")); design <- read.csv(file.path(ext, "design.csv"))
  one <- dat$one == 1; cts <- as.matrix(dat[one, -(1:14)]); rownames(cts) <- dat$HOG[one]
  D <- data.table::as.data.table(design); D[, col := Name][, num := as.integer(sub("^[A-Za-z]+", "", Name))]
  HYB <- c("ScxSk","ScxSp","ScxSu","SkxSu","SpxSk","SpxSu")
  H <- D[Strain %in% HYB & Condition == "heat"][, Time := as.integer(Time)]
  lapply(setNames(HYB, HYB), function(hy) {
    hs <- H[Strain == hy]; y <- calcNormFactors(DGEList(counts = cts[, hs$col]), method = "TMM"); lc <- cpm(y, log = TRUE, prior.count = 2)
    L1 <- hs[Allele == substr(hy, 1, 2)][order(num)]; L2 <- hs[Allele == substr(hy, 4, 5)][order(num)]
    stopifnot(identical(L1$num, L2$num))
    Dm <- lc[, L2$col] - lc[, L1$col]; des <- model.matrix(~ 0 + factor(L1$Time))
    list(Dm = Dm, des = des, R = Dm - Dm %*% des %*% solve(crossprod(des)) %*% t(des))
  })
})
.AR_CACHE <- new.env()
allele_rho <- function(hogs) {
  key <- paste(length(hogs), sum(utf8ToInt(paste(sort(hogs), collapse = ""))) %% 1e9)
  if (!is.null(.AR_CACHE[[key]])) return(.AR_CACHE[[key]])
  set.seed(871)
  v <- sapply(.AR, function(x) { g <- intersect(hogs, rownames(x$Dm)); if (length(g) < 5) return(c(NA, NA))
    if (length(g) > 1500) g <- sample(g, 1500)
    rho <- interGeneCorrelation(x$Dm[g, ], x$des)$correlation
    pr <- cbind(sample(g, 8000, TRUE), sample(g, 8000, TRUE)); pr <- pr[pr[, 1] != pr[, 2], , drop = FALSE]
    r <- vapply(seq_len(nrow(pr)), function(k) cor(x$R[pr[k, 1], ], x$R[pr[k, 2], ]), 0)
    c(max(0, rho), max(0, mean(r^2, na.rm = TRUE) - 1 / 4)) })
  out <- c(rho = mean(v[1, ], na.rm = TRUE), rho2 = mean(v[2, ], na.rm = TRUE)); .AR_CACHE[[key]] <- out; out
}

#!/usr/bin/env Rscript
## 94_p1_anchor_diagnosis.R -- diagnosis of the reference-arm reliability. The per-hybrid ranks of the measured
## reference-arm h_arm (91) do not match the ranks implied by Table 1 under IV validity
## (naive/gamma_IV): e.g. ScxSp measured 0.79 but implied ~0.94; SpxSk measured 0.95 but
## implied ~0.84. This discrepancy was specified for evaluation before the analysis was run.
## This diagnoses it
## once: recompute h_arm on Table 1's own gene set (model A; padj<0.05 filter, and the
## no-filter set), 15 min only, per hybrid, and compare with naive/gamma_IV. The analysis is
## evaluated without conditioning on agreement.
suppressPackageStartupMessages({ library(data.table); library(edgeR) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
source(file.path(root, "pkg/sharedstat/R/reliability.R"))
ext<-file.path(root,"data/03_external/fay2023"); csd<-file.path(root,"data/02_processed/crossspecies")
HYB<-c("ScxSk","ScxSp","ScxSu","SkxSu","SpxSk","SpxSu"); METHODS<-c("DESeq2","edgeR","limma")
ALPHA<-0.05; PRIOR<-2; MINCPM<-1
partner<-c(ScxSp="SkxSu",SkxSu="ScxSp",ScxSu="SpxSk",SpxSk="ScxSu",ScxSk="SpxSu",SpxSu="ScxSk")
say <- function(...) cat(sprintf(...), "\n", sep = "")

## --- Table 1 gene set: the genes 11_gain_iv.R keeps, per hybrid (union over methods) ---
X<-fread(cmd=paste("zcat",shQuote(file.path(csd,"fay_multimethod_p1order.tsv.gz"))))   # SpxSk refitted with the first-named allele (Sp) as reference
X<-X[condition=="heat" & group%in%HYB]
Rtab<-X[model=="A" & term=="time",.(gene,group,method,R=eff,padj_t=padj)]
## genes with reference-arm padj<0.05 (the filter of 11_gain_iv.R), and all finite (no filter)
tbl1_genes <- function(g, filtered=TRUE){
  d<-Rtab[group==g & is.finite(R)]
  if(filtered) d<-d[!is.na(padj_t) & padj_t<ALPHA]
  unique(d$gene)
}
## implied reference-arm h under IV validity = naive_gain / gain_iv = (1+slope_ols)/(1+slope_iv)
gi<-fread(file.path(csd,"gain_iv_p1order.tsv"))[model=="A"]   # first-named allele as reference
med<-function(g,col) median(gi[group==g,get(col)])
implied<-sapply(HYB,function(g){ (1+med(g,"slope_ols"))/(1+med(g,"slope_iv")) })

## --- per-hybrid logCPM, 15 min only, reference-arm h on Table 1's gene set ---
D<-as.data.table(read.csv(file.path(ext,"design.csv"))); D[,col:=Name]
D[,num:=as.integer(sub("^[A-Za-z]+","",Name))]; D[,block:=fifelse(num<=32L,"block1","block2")]
Hd<-D[Strain%in%HYB & Condition=="heat"]; Hd[,Time:=as.integer(Time)]
dat<-read.csv(file.path(ext,"htseq.ortho.csv")); one<-dat$one==1
cts<-as.matrix(dat[one,-(1:14)]); rownames(cts)<-dat$HOG[one]
HB<-fread(file.path(csd,"hog_to_scer.tsv"))
## Table 1 genes (fay_multimethod) are already HOG ids, matching htseq.ortho rownames directly.

h_arm_on <- function(g, genes_hog, tp=15L){
  hs<-Hd[Strain==g][order(Allele,block,Time)]; a1<-substr(g,1,2)
  m<-cts[,hs$col,drop=FALSE]; y<-calcNormFactors(DGEList(counts=m),method="TMM")
  lcpm<-cpm(y,log=TRUE,prior.count=PRIOR)
  gv<-function(bk,t){cc<-hs[Allele==a1&block==bk&Time==t,col];stopifnot(length(cc)==1L);lcpm[,cc]}
  rA1<-gv("block1",tp)-gv("block1",0L); rA2<-gv("block2",tp)-gv("block2",0L)
  hog<-intersect(genes_hog, rownames(lcpm))
  b1<-rA1[hog]; b2<-rA2[hog]; ok<-is.finite(b1)&is.finite(b2)
  if(sum(ok)<10) return(c(h=NA,n=sum(ok)))
  r<-reliability_from_blocks(b1[ok],b2[ok],center=TRUE); c(h=r$h,n=r$n_pairs)
}

say("=== anchor diagnosis: reference-arm h on Table 1's OWN gene set, 15 min, per hybrid ===")
say("%-6s %10s %12s %12s %10s","hybrid","implied","h_arm(filt)","h_arm(nofilt)","n_filt")
res<-list()
for(g in HYB){
  gf<-tbl1_genes(g,TRUE); gn<-tbl1_genes(g,FALSE)
  hf<-h_arm_on(g,gf); hn<-h_arm_on(g,gn)
  res[[g]]<-data.table(hybrid=g,implied=implied[[g]],h_filt=hf["h"],h_nofilt=hn["h"],n_filt=hf["n"])
  say("%-6s %10.4f %12.4f %12.4f %10d",g,implied[[g]],hf["h"],hn["h"],as.integer(hf["n"]))
}
R<-rbindlist(res)
say("\n=== rank comparison (1 = smallest) ===")
rk<-function(v) rank(v,ties.method="first")
say("implied ranks     : %s", paste(rk(R$implied),collapse=" "))
say("h_arm(filt) ranks : %s", paste(rk(R$h_filt),collapse=" "))
say("h_arm(nofilt)ranks: %s", paste(rk(R$h_nofilt),collapse=" "))
say("Spearman implied vs h_arm(filt): %.3f", cor(rk(R$implied),rk(R$h_filt)))
say("Spearman implied vs h_arm(nofilt): %.3f", cor(rk(R$implied),rk(R$h_nofilt)))
fwrite(R,file.path(root,"data/02_processed/hsf1/p1_anchor_diagnosis.tsv"),sep="\t")
say("=== done (evaluated without conditioning on agreement) ===")

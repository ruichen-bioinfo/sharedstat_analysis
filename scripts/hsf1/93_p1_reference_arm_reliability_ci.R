#!/usr/bin/env Rscript
## 93_p1_reference_arm_reliability_ci.R -- interval supplement to script 91 (reference-arm reliability).
## Reference-arm h bootstrap CI + per-library inter-gene correlation VIF, done exactly as
## script 82: 8 libraries per hybrid (4 times x 2 blocks), residuals from a
## time-means fit (4 residual d.f.), using the first-named allele's per-library logCPM.
## Then reversal-regime counts under point / rho-bar-VIF / rho^2-VIF, beside the DoD counts.
suppressPackageStartupMessages({ library(data.table); library(edgeR); library(limma) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
source(file.path(root, "pkg/sharedstat/R/reliability.R"))
ext<-file.path(root,"data/03_external/fay2023"); csd<-file.path(root,"data/02_processed/crossspecies")
hsd<-file.path(root,"data/02_processed/hsf1")
ANCHOR<-"37_30"; HYB<-c("ScxSp","ScxSk","ScxSu","SpxSk","SpxSu","SkxSu"); TP<-c(15L,30L,60L)
PRIOR<-2; MINCPM<-1; NB<-400; set.seed(20260928)
say <- function(...) cat(sprintf(...), "\n", sep = "")

dat<-read.csv(file.path(ext,"htseq.ortho.csv")); design<-read.csv(file.path(ext,"design.csv"))
one<-dat$one==1; cts<-as.matrix(dat[one,-(1:14)]); rownames(cts)<-dat$HOG[one]
D<-as.data.table(design); D[,col:=Name]; D[,num:=as.integer(sub("^[A-Za-z]+","",Name))]
D[,block:=fifelse(num<=32L,"block1","block2")]
H<-D[Strain%in%HYB & Condition=="heat"]; H[,Time:=as.integer(Time)]
POL<-data.table(hybrid=HYB,a1=substr(HYB,1,2),a2=substr(HYB,4,5))
HB<-fread(file.path(csd,"hog_to_scer.tsv"))
CS<-fread(file.path(hsd,sprintf("hsf1_core_set_%s.tsv",ANCHOR)))
CL<-fread(file.path(hsd,sprintf("hsf1_epistasis_classes_%s.tsv",ANCHOR)))
UC<-dcast(CL[,.(gene,method,class)],gene~method,value.var="class")
UC[,unanimous:=fifelse(DESeq2==edgeR&edgeR==limma,DESeq2,NA_character_)]
map<-function(g){h<-HB[Name%in%g|sys%in%g,unique(gene)];h[h%in%rownames(cts)]}
STR<-list(Hsf1_core=CS[set=="Hsf1_core_unanimous",gene],Hsf1_joint=CS[set=="Hsf1_joint_candidate",gene],
  Msn24_only=UC[unanimous=="Msn24_only",gene],TF_independent=UC[unanimous=="TF_independent",gene],
  other_induced=UC[unanimous%in%c("Msn24_only","TF_independent"),gene],not_induced=UC[unanimous=="not_induced",gene])
HOGS<-lapply(STR,map)

rows<-list()
for(g in HYB){
  hs<-H[Strain==g][order(Allele,block,Time)]; a1<-POL[hybrid==g,a1]; a2<-POL[hybrid==g,a2]
  m<-cts[,hs$col,drop=FALSE]; y<-calcNormFactors(DGEList(counts=m),method="TMM")
  lcpm<-cpm(y,log=TRUE,prior.count=PRIOR); craw<-cpm(y)
  k1<-rowSums(craw[,hs[Allele==a1,col],drop=FALSE]>=MINCPM)>=4L
  k2<-rowSums(craw[,hs[Allele==a2,col],drop=FALSE]>=MINCPM)>=4L; keep<-k1&k2
  gv<-function(al,bk,t){cc<-hs[Allele==al&block==bk&Time==t,col];stopifnot(length(cc)==1L);lcpm[,cc]}
  ## --- per-library reference-arm logCPM (8 libraries), as in script 82 ---
  ## The reference-arm reliability is of the a1 heat response measured across libraries;
  ## the relevant correlation is that of a1's per-library logCPM, residualised on time
  ## (4 groups incl. t=0, giving 4 residual d.f. over 8 libraries), exactly as 82 does for the
  ## allele-difference contrast Dm.
  L1 <- hs[Allele==a1][order(num)]                       # 8 libraries of the reference allele
  stopifnot(nrow(L1)==8L)
  A1m <- lcpm[, L1$col]; colnames(A1m) <- paste0("lib",L1$num)
  des <- model.matrix(~ 0 + factor(L1$Time))
  Rres <- A1m - A1m %*% des %*% solve(crossprod(des)) %*% t(des)   # residuals, 4 d.f.
  for(t in TP){
    rA1<-gv(a1,"block1",t)-gv(a1,"block1",0L); rA2<-gv(a1,"block2",t)-gv(a1,"block2",0L)
    for(st in names(HOGS)){
      gg<-intersect(HOGS[[st]],rownames(lcpm)[keep]); if(length(gg)<10) next
      b1<-rA1[gg]; b2<-rA2[gg]; ok<-is.finite(b1)&is.finite(b2); b1<-b1[ok]; b2<-b2[ok]; ggk<-gg[ok]; n<-length(b1)
      if(n<10) next
      r<-reliability_from_blocks(b1,b2,center=TRUE)
      hb<-replicate(NB,{k<-sample(n,replace=TRUE);reliability_from_blocks(b1[k],b2[k],center=TRUE)$h})
      ## rho-bar: mean inter-gene correlation of the a1 per-library residuals (limma), as 82
      gs<-if(n>2000) sample(ggk,2000) else ggk
      rho<-tryCatch(interGeneCorrelation(A1m[gs,,drop=FALSE], des)$correlation, error=function(e) NA_real_)
      ## rho^2 (bias-corrected mean r^2) form, as 82: pairs of residual rows, mean r^2 minus 1/df_pairs baseline (1/4)
      pr<-matrix(sample(gs,2*min(20000,choose(length(gs),2)),replace=TRUE),ncol=2); pr<-pr[pr[,1]!=pr[,2],]
      rr<-vapply(seq_len(nrow(pr)),function(k) cor(Rres[pr[k,1],],Rres[pr[k,2],]),0)
      rho2_bc<-max(0, mean(rr^2,na.rm=TRUE)-1/4)
      rows[[length(rows)+1L]]<-data.table(stratum=st,hybrid=g,time=t,n_genes=n,h=r$h,
        lo=unname(quantile(hb,0.025)),hi=unname(quantile(hb,0.975)),
        rho=rho, rho2_bc=rho2_bc)
    }
  }
}
B<-rbindlist(rows)
B[,rho_pos:=pmax(rho,0)]; B[is.na(rho_pos),rho_pos:=0]
B[,vif :=1+(n_genes-1)*rho_pos]
B[,vif2:=1+(n_genes-1)*rho2_bc]
B[,hi_adj :=pmin(1,h+sqrt(vif )*(hi-h))]
B[,hi_adj2:=pmin(1,h+sqrt(vif2)*(hi-h))]
fwrite(B,file.path(hsd,"p1_reference_arm_h_adjusted.tsv"),sep="\t")
E<-B[is.finite(h)&h>0]
say("=== reference-arm h reversal counts: point vs correlation-adjusted (per-library) ===")
say("estimable cells: %d (6 strata incl. other_induced x 6 hybrids x 3 tp)", nrow(E))
say("inter-gene rho-bar (a1 per-lib residuals): median %.4f (range %.4f-%.4f); VIF median %.1f",
    median(E$rho,na.rm=TRUE),min(E$rho,na.rm=TRUE),max(E$rho,na.rm=TRUE),median(E$vif))
say("rho^2 (bias-corrected) median %.4f; VIF2 median %.1f", median(E$rho2_bc), median(E$vif2))
say("\n%-6s %8s %12s %14s %14s %14s","gamma","1/gamma","point","interval(no adj)","rho-bar VIF","rho^2 VIF")
for(g in c(1.1,1.2,1.3,1.5,2.0))
  say("%-6.1f %8.4f %12d %14d %14d %14d",g,1/g,
      sum(E$h<1/g), sum(E$hi<1/g), sum(E$hi_adj<1/g), sum(E$hi_adj2<1/g))
say("(interval counts = whole CI below 1/gamma; the adjusted forms widen the CI by sqrt(VIF))")
say("=== done ===")

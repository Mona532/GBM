# ============================================================
# tcga_ilc_survival.R — ILC-NK signature survival in TCGA-GBM
# ============================================================
# Uses top10 DEGs from ILC-NK program as gene signature,
# scores TCGA-GBM samples, finds optimal cutoff, plots KM curve.
# Expression: UCSC Xena HiSeqV2 (user-provided file)
# Clinical: downloaded from UCSC Xena
# ============================================================

suppressMessages({
  library(survival)
  library(survminer)
  library(maxstat)
})

# ---- 1. Read TCGA expression ----
expr_file <- "C:/Users/qingy/Downloads/TCGA.GBM.sampleMap_HiSeqV2_PANCAN.gz"
expr <- read.table(gzfile(expr_file), header=TRUE, sep="\t", row.names=1, check.names=FALSE)
cat(sprintf("Expression: %d genes x %d samples\n", nrow(expr), ncol(expr)))

# ---- 2. Top10 ILC-NK DEGs (from posdev_all_DEGs.csv) ----
ilc_genes <- c("BTBD17","LINC02798","PNCK","CRMP1","FABP7",
               "ADCYAP1R1","CX3CL1","EGFR","LUZP2","F3")
present <- intersect(ilc_genes, rownames(expr))
cat(sprintf("ILC genes present: %d/%d\n", length(present), length(ilc_genes)))
if (length(present) < 5) stop("Too few ILC genes in TCGA expression matrix")

# ---- 3. Download TCGA-GBM clinical data from UCSC Xena ----
clinical_url <- "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.GBM.sampleMap%2FGBM_clinicalMatrix"
clinical_file <- tempfile(fileext=".tsv")
download.file(clinical_url, clinical_file, mode="wb")
clinical <- read.table(clinical_file, header=TRUE, sep="\t", check.names=FALSE)
cat(sprintf("Clinical: %d patients\n", nrow(clinical)))

# ---- 4. Compute ILC signature score per sample ----
ilc_score <- colMeans(expr[present, , drop=FALSE], na.rm=TRUE)
score_df <- data.frame(
  sample = names(ilc_score),
  ilc_score = ilc_score,
  # map sample ID to patient ID (TCGA-XX-XXXX-XX -> TCGA-XX-XXXX)
  patient = sub("^(TCGA-[0-9]{2}-[0-9]{4}).*", "\\1", names(ilc_score)),
  stringsAsFactors = FALSE
)
# average across technical replicates (different -01/-02 suffix)
score_patient <- aggregate(ilc_score ~ patient, data=score_df, FUN=mean)

# ---- 5. Merge with survival data ----
clinical$patient <- sub("^(TCGA-[0-9]{2}-[0-9]{4}).*", "\\1", clinical$sampleID)
surv_df <- merge(score_patient, clinical, by="patient")
cat(sprintf("Merged: %d patients with expression and survival\n", nrow(surv_df)))

# ---- 6. Build OS and PFS ----
surv_df$OS.time <- as.numeric(surv_df$CDE_survival_time)
surv_df$OS      <- ifelse(surv_df$CDE_vital_status=="DECEASED", 1, 0)
surv_df$PFS.time <- ifelse(!is.na(surv_df$days_to_tumor_recurrence) & surv_df$days_to_tumor_recurrence>0,
                            surv_df$days_to_tumor_recurrence, surv_df$days_to_last_followup)
surv_df$PFS <- ifelse(!is.na(surv_df$days_to_tumor_recurrence) & surv_df$days_to_tumor_recurrence>0, 1, 0)
surv_df <- surv_df[!is.na(surv_df$OS.time) & surv_df$OS.time>0, ]
cat(sprintf("Filtered: %d patients, OS events=%d, PFS events=%d\n",
            nrow(surv_df), sum(surv_df$OS), sum(surv_df$PFS,na.rm=TRUE)))

# ---- 7. KM for OS and PFS using optimal cutoff ----
for(endpoint in c("OS","PFS")){
  tvar<-paste0(endpoint,".time");evar<-endpoint
  s<-surv_df[!is.na(surv_df[[tvar]]) & surv_df[[tvar]]>0,]
  cut<-maxstat.test(Surv(s[[tvar]],s[[evar]])~ilc_score,data=s,smethod="LogRank",pmethod="Lau94")
  best<-cut$estimate;s$group<-ifelse(s$ilc_score>best,"ILC-high","ILC-low")
  cat(sprintf("%s: cutoff=%.3f high=%d low=%d\n",endpoint,best,sum(s$group=="ILC-high"),sum(s$group=="ILC-low")))
  fit<-survfit(Surv(s[[tvar]],s[[evar]])~group,data=s)
  pval<-surv_pvalue(fit)$pval;cat(sprintf("  log-rank p=%.2e\n",pval))
  pdf(file.path("E:/GBM/results",paste0("fig_tcga_ilc_",endpoint,"_km.pdf")),width=4.5,height=4.5)
  p<-ggsurvplot(fit,data=s,pval=TRUE,pval.size=3.5,palette=c("#E74C3C","#3498DB"),
    legend.title="ILC signature",legend.labs=c("High","Low"),
    xlab=ifelse(endpoint=="OS","Overall survival (days)","Progression-free survival (days)"),
    ylab="Survival probability",ggtheme=theme_classic(base_size=11,base_family="sans"),
    risk.table=TRUE,risk.table.height=0.25,break.time.by=500,conf.int=FALSE)
  print(p,newpage=FALSE);dev.off()
  jpeg(file.path("E:/GBM/results",paste0("fig_tcga_ilc_",endpoint,"_km.jpg")),width=4.5,height=4.5,units="in",res=400,quality=95)
  print(p,newpage=FALSE);dev.off()
}
cat("Done\n")
                     smethod="LogRank", pmethod="Lau94")
best_cut <- cut$estimate
cat(sprintf("Optimal cutoff: %.3f\n", best_cut))
surv_df$group <- ifelse(surv_df$ilc_score > best_cut, "ILC-high", "ILC-low")

# ---- 7. Kaplan-Meier ----
fit <- survfit(Surv(OS.time, OS) ~ group, data=surv_df)
pval <- surv_pvalue(fit)$pval

pdf(file.path("E:/GBM/results", "fig_tcga_ilc_km.pdf"), width=4.5, height=4.5)
p <- ggsurvplot(fit, data=surv_df,
                pval=TRUE, pval.size=3.5,
                palette=c("#E74C3C","#3498DB"),
                legend.title="ILC signature", legend.labs=c("High","Low"),
                xlab="Overall survival (days)", ylab="Survival probability",
                ggtheme=theme_classic(base_size=11, base_family="sans"),
                risk.table=TRUE, risk.table.height=0.25,
                break.time.by=500, conf.int=FALSE)
print(p, newpage=FALSE)
dev.off()

jpeg(file.path("E:/GBM/results", "fig_tcga_ilc_km.jpg"), width=4.5, height=4.5, units="in", res=400, quality=95)
print(p, newpage=FALSE)
dev.off()

cat(sprintf("Saved KM plot. High n=%d, Low n=%d, log-rank p=%.2e\n",
            sum(surv_df$group=="ILC-high"), sum(surv_df$group=="ILC-low"), pval))
cat("Done\n")

suppressMessages({library(survival);library(survminer);library(maxstat)})

# ---- Read CGGA data ----
expr<-read.table("E:/GBM/data/cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt",header=TRUE,sep="\t",row.names=1,check.names=FALSE)
clinical<-read.table("E:/GBM/data/cgga/CGGA.mRNAseq_693_clinical.20200506.txt",header=TRUE,sep="\t",check.names=FALSE)
cat(sprintf("Expr: %d genes x %d samples\n",nrow(expr),ncol(expr)))
cat(sprintf("Clinical: %d patients\n",nrow(clinical)))

# ---- Filter GBM only ----
clinical$Grade<-trimws(clinical$Grade)
gbm<-clinical[grepl("GBM|glioblastoma|IV|4",clinical$Histology,ignore.case=TRUE) |
              grepl("GBM|glioblastoma|IV|4",clinical$Grade,ignore.case=TRUE),]
cat(sprintf("GBM patients: %d\n",nrow(gbm)))

# ---- ILC signature (top10 ILC-NK DEGs) ----
ilc_genes<-c("BTBD17","LINC02798","PNCK","CRMP1","FABP7",
             "ADCYAP1R1","CX3CL1","EGFR","LUZP2","F3")
present<-intersect(ilc_genes,rownames(expr))
cat(sprintf("ILC genes present: %d/%d\n",length(present),length(ilc_genes)))

# ---- Score GBM samples ----
expr_gbm<-expr[,intersect(colnames(expr),gbm$CGGA_ID),drop=FALSE]
ilc_score<-colMeans(expr_gbm[present,,drop=FALSE],na.rm=TRUE)
score_df<-data.frame(CGGA_ID=names(ilc_score),ilc_score=ilc_score,stringsAsFactors=FALSE)
surv_df<-merge(score_df,gbm,by="CGGA_ID")

# ---- Survival ----
surv_df$OS.time<-as.numeric(surv_df$OS)
surv_df$OS<-as.numeric(surv_df[["Censor (alive=0; dead=1)"]])
surv_df<-surv_df[!is.na(surv_df$OS.time)&surv_df$OS.time>0,]
cat(sprintf("GBM with survival: %d, events=%d\n",nrow(surv_df),sum(surv_df$OS)))

# ---- KM with optimal cutoff ----
cut<-maxstat.test(Surv(OS.time,OS)~ilc_score,data=surv_df,smethod="LogRank",pmethod="Lau94")
best<-cut$estimate;surv_df$group<-ifelse(surv_df$ilc_score>best,"ILC-high","ILC-low")
cat(sprintf("Optimal cutoff: %.3f, high=%d low=%d\n",best,sum(surv_df$group=="ILC-high"),sum(surv_df$group=="ILC-low")))
fit<-survfit(Surv(OS.time,OS)~group,data=surv_df)
pval<-surv_pvalue(fit)$pval;cat(sprintf("log-rank p=%.2e\n",pval))

root<-"E:/GBM/results"
pdf(file.path(root,"fig_cgga_ilc_OS_km.pdf"),width=4.5,height=4.5)
p<-ggsurvplot(fit,data=surv_df,pval=TRUE,pval.size=3.5,palette=c("#E74C3C","#3498DB"),
  legend.title="ILC signature",legend.labs=c("High","Low"),
  xlab="Overall survival (days)",ylab="Survival probability",
  ggtheme=theme_classic(base_size=11,base_family="sans"),
  risk.table=TRUE,risk.table.height=0.25,break.time.by=500,conf.int=FALSE)
print(p,newpage=FALSE);dev.off()
jpeg(file.path(root,"fig_cgga_ilc_OS_km.jpg"),width=4.5,height=4.5,units="in",res=400,quality=95)
print(p,newpage=FALSE);dev.off()
cat("Done\n")

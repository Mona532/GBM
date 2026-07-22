# Export hires images from Seurat RDS for verification
suppressPackageStartupMessages({library(Seurat);library(png)})
root<-"E:/GBM/results/visium_rds"
out <-file.path(root,"hires_export")
dir.create(out,showWarnings=FALSE)
files<-list.files(root,pattern="\\.rds$",recursive=TRUE,full.names=TRUE)
cat(sprintf("Found %d RDS files\n",length(files)))
for(rds in files){
  sid<-sub("\\.rds$","",basename(rds));se<-readRDS(rds)
  if(length(se@images)>0&&!is.null(se@images[[1]]@image)){
    img<-se@images[[1]]@image
    writePNG(img,file.path(out,paste0(sid,"_hires.png")))
    cat(sprintf("%s: %dx%d\n",sid,ncol(img),nrow(img)))
  }else cat(sprintf("%s: no image\n",sid))
}
cat("Done\n")

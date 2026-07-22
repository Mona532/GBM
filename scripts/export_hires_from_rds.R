# Copy original hires images from spatial/ directories
root<-"E:/GBM/results/visium_rds"
out <-file.path(root,"hires_export")
dir.create(out,showWarnings=FALSE)
dirs<-list.files(root,full.names=TRUE)
for(d in dirs){
  if(!file.exists(file.path(d,"spatial","tissue_hires_image.png"))) next
  sid<-basename(d)
  file.copy(file.path(d,"spatial","tissue_hires_image.png"),
            file.path(out,paste0(sid,"_hires.png")),overwrite=TRUE)
}
cat(sprintf("Copied %d hires images\n",length(list.files(out))))

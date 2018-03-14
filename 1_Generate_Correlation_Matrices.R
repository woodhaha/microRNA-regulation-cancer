library("tictoc")
library("Hmisc")
tic("total runtime")
tic("data load")
types<-c("COAD","HNSC","OV","UCEC","CESC")
temp<-lapply(types,function(i){
  setwd("/DDP/cancer-ext/data-paired/")
  assign(paste0(i,"_rna"),read.csv(paste0(i,"-rna.csv"),header = TRUE)[,-1],envir = globalenv())
  assign(paste0(i,"_mir"),read.csv(paste0(i,"-mir.csv"),header = TRUE)[,-1],envir = globalenv())
})

index.34b<-sapply(types,function(i){
  which(get(paste0(i,"_mir"))[,1]== "hsa-mir-34b")
})
rm(temp)
toc()
tic("Generating correlation matrix and file write")
write.table("",file = "miR-34b-corr-fpkm.txt",quote=F,append = FALSE,sep = ",",col.names = F,row.names = F)
lapply(1:length(types), function(i){
  tic(types[i])
  name<-types[i]
  sub<-index.34b[i]
  x<-t(get(paste0(name,"_mir"))[sub,-1])
  y<-get(paste0(name,"_rna"))[,1]
  corr<-sapply(1:nrow(get(paste0(name,"_rna"))[,-1]),function(j){
    COR<-rcorr(x,(as.numeric(get(paste0(name,"_rna"))[j,-1])),type="spearman")
    COR$r[1,2]
  })
  pvals<-sapply(1:nrow(get(paste0(name,"_rna"))[,-1]),function(j){
    COR<-rcorr(x,t(get(paste0(name,"_rna"))[j,-1]),type="spearman")
    COR$P[1,2]
  })
  z<-cbind(as.character(y),corr,pvals)
  colnames(z)<-c("gene","correlation","p-value")
  write.table(paste0("#",name,collapse = ","),file = "miR-34b-corr-fpkm.txt",append= T,sep = ",",quote=F,col.names = F,row.names = F)
  z<-z[which(pvals<0.05),]
  write.table(z,file = "miR-34b-corr-fpkm.txt",append = T,quote=F,sep=",")
  toc()
})
toc()
toc()
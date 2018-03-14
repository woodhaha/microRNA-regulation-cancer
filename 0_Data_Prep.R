
# library import ----------------------------------------------------------
library("tictoc")
library("jsonlite")
sink("data-paired/logs.txt",append = TRUE)
print(Sys.time())
rm(list=ls()) #clearing workspace
tic("total") #run time is approx 3900s
types<-c("COAD","HNSC","OV","UCEC","CESC")

barcode_filter<-20
#looping over cancer data types

lapply(types,function(type){
  tic(i)
  print(paste0("cancer type: ",type))

# metadata import ---------------------------------------------------------

  metadata<-fromJSON(paste0("data/",type,"-metadata.json"))
  print("metadata contents")
  print(colnames(metadata))#use associated_entities
  x<- as.character(sapply(metadata$associated_entities,function(i){
    i[3]
  }))
  x<-substring(x,1,barcode_filter)
  meta<-cbind(metadata[,c(2,3,5)],x)
  y<-which(meta$data_format=="TXT")
  meta[y,2]<-substring(meta[y,2],1,45)#Should not need changing
  x<-c(x[y],x[-y])
  meta<-rbind(meta[y,],meta[-y,])
  
  
  print(" filename splicing check")
  print(sum(sapply(meta$data_format,function(i){substring(i,nchar(i)-3)=="TXT"}))==length(y))
  print(" sort and contents check")
  pint(identical(1:length(y),which(meta$data_format=="TXT")))
  flush.console()

# unpaired data import ----------------------------------------------------
  
  rna<-lapply(meta[1:length(y),2], function(i){
    read.table(file=paste("data/",type,"/rna/",i,sep = "",collapse = ""),sep="\t")
  })
  mir<-lapply(meta[(length(y)+1):length(x),2],function(i){
    read.table(file=paste("data/",type,"/mir/",i,sep = "",collapse = ""),sep="\t",header = TRUE)
  })
  
  z<-Reduce(cbind,rna)
  rna<-z[,c(1,seq(2,ncol(z),2))]
  z<-Reduce(cbind,mir)
  mir<-z[,c(1,seq(3,ncol(z),4))]
  rm(z)
  colnames(rna)[-1] <- x[1:length(y)]
  colnames(mir)[-1] <- x[(length(y)+1):length(x)]
  raw_rna<-rna
  raw_mir<-mir

# pairing of samples ------------------------------------------------------

  n<-c(sum(meta$data_format=="TXT"),length(x))
  y<-x[which(duplicated(x))]
  z1<-sapply(y,function(i){
    which(x[1:n[1]]==i)
  })
  z2<-sapply(y,function(i){
    which(x[(n[1]+1):n[2]]==i)
  })
  print(cbind(z1,z2))
  temp<-unique(c(which(as.integer(sapply(z1,length))==0),which(as.integer(sapply(z2,length))==0)))
  if(!identical(temp,integer())){y<-y[-temp]} else{print("No unpaired double IDs")}
  z<-(lapply(y,function(i){
    which(x==i)
  }))
  z<-z[which(sapply(z, length)==2)]
  z<-matrix(unlist(z,recursive = TRUE),ncol = 2,byrow = TRUE)
  rna<-raw_rna[,c(1,1+z[,1])]
  mir<-raw_mir[,c(1,(z[,2]-(n[1]-1)))]
  rm(x,y,z)
  print("pairing check")
  print(identical(colnames(mir[,-1]),colnames(rna[,-1])))
  ref<-read.csv(file="ref.mapped.csv",header = TRUE)
  rna<-rna[as.integer(ref[,1]),]
  print("ID mapping check")
  print(identical(as.character(ref[,2]),substring(as.character(rna[,1]),1,15)))
  rna[,1]<-ref[,3]
  flush.console()
  
# filtering data ----------------------------------------------------------

  x<-sapply(1:nrow(rna),function(i){
    sum(rna[i,-1]==0)
  })
  rna<-rna[-which(x>(ncol(rna[,-1])/4)),]
  x<-sapply(1:nrow(mir),function(i){
    sum(mir[i,-1]==0)
  })
  mir<-mir[-which(x>(ncol(mir[,-1])/2)),]

# file write --------------------------------------------------------------

  write.csv(rna,file = paste0("data-paired/",type,"-rna.csv"))
  write.csv(mir,file = paste0("data-paired/",type,"-mir.csv"))
  
  toc()
  flush.console()
})
toc()
print(Sys.time())
rm(list=ls())
sink()

# importing libraries ----------------------------------------------------------
library("tictoc") #timing functions tic() and toc()
library("jsonlite") #fromJSON() required to load json files
sink("logs/0_Data_Prep.txt",split=TRUE) #printing output to file and console
print(Sys.time()) #printing runtime
rm(list=ls()) #clearing workspace
tic("total") #run time is approx 3900s
types<-readLines("types")

barcode_filter<-20 # splicing TCGA barcode before plate/aliquot

#looping over cancer data types
lapply(types,function(type){
  tic(type)
  print(paste0("cancer type: ",type))

# metadata import ---------------------------------------------------------

  metadata<-fromJSON(paste0("fpkm/data/",type,"_metadata.json"))
  print("metadata contents")
  print(colnames(metadata))#use associated_entities for barcode
  x<- as.character(sapply(metadata$associated_entities,function(i){
    i[3] #barcode is element 3 of associated_entities
  }))
  x<-substring(x,1,barcode_filter) #splicing barcode 
  meta<-cbind(metadata[,c(2,3,5)],x) #filtering metadata columns to 
  y<-which(meta$data_format=="TXT") # only miRNA quantification files end with .txt in metadata
  meta[y,2]<-substring(meta[y,2],1,45)#filename splice to get rid of .gz compression extension
  x<-c(x[y],x[-y]) #ordering barcode
  meta<-rbind(meta[y,],meta[-y,])	#ordering metadata and adding spliced barcode
  
  
  print(" filename splicing check")
  print(sum(sapply(meta$data_format,function(i){substring(i,nchar(i)-3)=="TXT"}))==length(y))
  print(" sort and contents check")
  print(identical(1:length(y),which(meta$data_format=="TXT")))
  flush.console()

# unpaired data import ----------------------------------------------------
  
  rna<-lapply(meta[1:length(y),2], function(i){
    read.table(file=paste("fpkm/data/",type,"/rna/",i,sep = "",collapse = ""),sep="\t")
  })
  mir<-lapply(meta[(length(y)+1):length(x),2],function(i){
    read.table(file=paste("fpkm/data/",type,"/mir/",i,sep = "",collapse = ""),sep="\t",header = TRUE)
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

  write.csv(rna,file = paste0("fpkm/data-paired/",type,"-rna.csv"))
  write.csv(mir,file = paste0("fpkm/data-paired/",type,"-mir.csv"))
  
  print(toc())
  })
toc()
print(Sys.time())
rm(list=ls())
sink()
library("tictoc")
library("Hmisc")

sink("logs/1_Generate_Correlation_Matrices.txt", split = TRUE, append = FALSE)
test <- c(TRUE,FALSE)[2]
Sys.time()
tic("total runtime")

# data import -------------------------------------------------------------

tic("data load")
types <- readLines("types")
setwd("fpkm/data-paired/")
invisible(lapply(types, function(i) {
  assign(paste0(i, "_rna"), read.csv(paste0(i, "-rna.csv"), header = TRUE)[, -1], envir = globalenv())
  assign(paste0(i, "_mir"), read.csv(paste0(i, "-mir.csv"), header = TRUE)[, -1], envir = globalenv())
}))
toc()
print("data imported into environment")
setwd("../..")

# common genes and miRs across cancers ------------------------------------

mir <-
  Reduce(
    intersect,
    list(
      as.character(CESC_mir[, 1]),
      as.character(HNSC_mir[, 1]),
      as.character(OV_mir[, 1]),
      as.character(UCEC_mir[, 1])
    )
  )
rna <-
  Reduce(
    intersect,
    list(
      as.character(CESC_rna[, 1]),
      as.character(HNSC_rna[, 1]),
      as.character(OV_rna[, 1]),
      as.character(UCEC_rna[, 1])
    )
  )

print("common genes and microRNAs through cancers calculated")

# correlating mir vs genes ------------------------------------------------

tic("Generating correlation matrix and file write")
print("starting correlation")

if(test){
  types <- types[1] # remove after testing
  rna <- rna[1:5] #remove comment for testing
}

invisible(lapply(1:length(types), function(i) {
  t <- proc.time()
  print(types[i])
  write.table(paste0(c("gene",mir),collapse="\t"),
              file = paste0("fpkm/correlations/",types[i],"miR-gene-corr-spearman.txt"),
              quote = F,
              append = FALSE,
              sep = "\t",
              col.names = F,
              row.names = F
  )
  write.table(paste0(c("gene",mir),collapse="\t"),
              file = paste0("fpkm/correlations/",types[i],"miR-gene-pvalue-spearman.txt"),
              quote = F,
              append = FALSE,
              sep = "\t",
              col.names = F,
              row.names = F
  )
  x <- get(paste0(types[i], "_mir"))
  #ordering row of mir and rna matrices
  x <-
    x[sapply(mir, function(j) {
      which(as.character(x[, 1]) == j)
    }),]
  y <- get(paste0(types[i], "_rna"))
  y <- y[sapply(rna, function(j) {
    which(as.character(y[, 1]) == j)[1] #some elements are repeated after ensembl to gene name conversion 
  }),]
  #starting correlation loop  
  sapply(1:(as.integer(length(rna)/1000)+1), function(j) {
      
      COR <- rcorr(t(y[(j*1000-999):(j*1000), -1]), t(x[, -1]), type = "spearman")
      rownames(COR$r)[1:1000]<-rna[(j*1000-999):(j*1000)]
      rownames(COR$P)[1:1000]<-rna[(j*1000-999):(j*1000)]
      write.table(COR$r[1:1000,1001:nrow(COR$r)],
                  file = paste0("fpkm/correlations/",types[i],"miR-gene-corr-spearman.txt"),
                  quote = F,
                  append = TRUE,
                  sep = "\t",
                  col.names = FALSE,
                  row.names = TRUE
      )
      write.table(COR$P[1:1000,1001:nrow(COR$P)],
                  file = paste0("fpkm/correlations/",types[i],"miR-gene-pvalue-spearman.txt"),
                  quote = F,
                  append = TRUE,
                  sep = "\t",
                  col.names = FALSE,
                  row.names = TRUE
      )      
      cat(paste0(min(length(rna),j*1000),"/",length(rna),"  "))
  })
  cat("\n")
  print(proc.time() - t)
}))

toc()
rm(list=ls())
library(reshape2)
tic("melt")
types <- readLines("types")
invisible(lapply(types, function(i) {
  assign(paste0(i, "_corr"), read.table(paste0("fpkm/correlations/",i,"miR-gene-corr-spearman.txt"),sep = "\t", header = TRUE), envir = globalenv())
  assign(paste0(i, "_pvalue"), read.csv(paste0("fpkm/correlations/",i,"miR-gene-pvalue-spearman.txt"),sep = "\t", header = TRUE), envir = globalenv())
}))
if(Reduce(all,lapply(types,function(i){
  identical(get(paste0(i,"_corr"))[,1],get(paste0(i,"_pvalue"))[,1])
}))){
  temp<-lapply(types,function(i){
    cbind(melt(get(paste0(i,"_corr"))),melt(get(paste0(i,"_pvalue")))[,3],i)
  })
  temp<-Reduce(rbind,temp)
  colnames(temp)<-c("gene","microRNA","correlation","pvalue","cancer")
  temp<-temp[which(as.numeric(temp$pvalue)< 0.001),]
  write.table(temp,"generated_files/correlations-spearman.tsv",row.names = F, col.names = T, quote = F, sep="\t" )
  
}
print("melt done")

toc()
toc()
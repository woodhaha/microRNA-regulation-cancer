rm(list=ls())

test <- F

progress <- function (x, max = 100) {
  percent <- x / max * 100
  cat(sprintf('\r[%-50s] %d%%',
              paste(rep('=', percent / 2), collapse = ''),
              floor(percent)))
  if (x == max)
    cat('\n')
}
# import data -------------------------------------------------------------
print(Sys.time())

# setwd("D:/DDP/scripts/")

library(readr)
library(dplyr)
library(tictoc)

tic("file reads")
correlations<-read_tsv("correlations-spearman.tsv")
correlations$microRNA<-gsub("\\.","-",correlations$microRNA)
correlations$microRNA<-gsub("mir","miR",correlations$microRNA)
correlations$pvalue<-p.adjust(correlations$pvalue, method = "bonferroni")
correlations <- filter(correlations, pvalue < 0.001, correlation < 0)

if(test){
  correlations <- filter(correlations, microRNA=="hsa-miR-34a" |  microRNA=="hsa-miR-34b" |  microRNA=="hsa-miR-34c")
  print("WARNING: TEST RUN")
}

MIRDIP<- read_csv("database_downloads/MIRDIP_MICRORNA_C-Version4.0.txt",
                  col_names = readLines("database_downloads/MIRDIP_Column_Names.txt")[-1])
MIRDIP$MICRORNA<-sub("-3p","",MIRDIP$MICRORNA)
MIRDIP$MICRORNA<-sub("-5p","",MIRDIP$MICRORNA)
MIRDIP<-MIRDIP[(MIRDIP$GENE_SYMBOL %in% unique(correlations$gene))&(MIRDIP$MICRORNA %in% unique(correlations$microRNA)),]
MIRDIP<-filter(MIRDIP,SCORE_CLASS=="E" | SCORE_CLASS=="G" | SCORE_CLASS=="F")
#MIRDIP SCORE_CLASS are B Bad, P Poor, F Favorable, G Good, E Excellent

CLASH<- read_tsv("database_downloads/CLASH.txt", skip= 30)[,c(2,6)]
CLASH$microRNA_name<-sub("_microRNA","",substring(CLASH$microRNA_name,22))
CLASH$microRNA_name<-sub("-5p","",CLASH$microRNA_name)
CLASH$microRNA_name<-sub("-3p","",CLASH$microRNA_name)
CLASH$microRNA_name<-sub("\\*","",CLASH$microRNA_name)
CLASH$mRNA_name <- sub("_mRNA","",substring(CLASH$mRNA_name,33))

toc()

# mapping prediction/confirmation onto correlation ------------------------
tic("mapping")

tic("CLASH")
print(Sys.time())
CLASH_C <-sapply(1:nrow(correlations),function(i){
  progress(i,nrow(correlations))
  nrow(filter(CLASH, microRNA_name == correlations$microRNA[i] & mRNA_name == correlations$gene[i])) != 0
})
toc()

tic("MIRDIP")
print(Sys.time())
mirDIP_P <-sapply(1:nrow(correlations),function(i){
  progress(i,nrow(correlations))
  nrow(filter(MIRDIP, MICRORNA == correlations$microRNA[i] & GENE_SYMBOL == correlations$gene[i])) != 0
})
toc()

correlations<- mutate(correlations, CLASH_C , mirDIP_P , duplicate=duplicated(paste0(correlations$gene,correlations$microRNA)))
paste0(as.integer(nrow(filter(correlations, mirDIP_P == TRUE| CLASH_C ==TRUE | duplicate == TRUE ))/(0.01*nrow(correlations))
),"% confirmed or predicted")
write_tsv(correlations,"generated_files/correlations-mapped.tsv")
write_csv(filter(correlations, mirDIP_P==TRUE | CLASH_C==TRUE | duplicate == TRUE),"generated_files/correlations_filtered.csv")
temp<- t(sapply(unique(correlations$microRNA),function(i){
  df<- filter(correlations, microRNA==i)
  (c(sum(unlist(df[,6]))/(0.01*nrow(df)),sum(unlist(df[,7]))/(0.01*nrow(df))))
}))
colnames(temp)<-c("CLASH_C","mirDIP_P")
write.csv(temp,"generated_files/percentage_predictions_confirmed.csv")
rm(temp)

toc()

# Stats on predictions ----------------------------------------------------

temp <- sapply(seq(min(correlations$correlation),max(correlations$correlation),0.01),function(i){
  nrow(filter(correlations, duplicate==TRUE & correlation < i & mirDIP_P== TRUE))*100/nrow(filter(correlations,duplicate==TRUE &  correlation < i))
})
ggplot(data.frame(x=seq(min(correlations$correlation),max(correlations$correlation),0.01),y=temp),aes(x=x,y=y))+
  geom_area()+theme_bw()+ scale_x_reverse() + ylim(0,100) +
  ggtitle("variation of percentage predicted vs total correlated using spearman") + xlab("correlation cutoff")+
  ylab("percentage of correlated genes predicted")

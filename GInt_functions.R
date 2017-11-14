#get data from hicpro/hiccompare
#install.packages("Amelia")
#sample code:
#library(Amelia)
#library(data.table)
library(InteractionSet)
library(readr)
library(HiCcompare)
#missmap(allhiccontacts)
# allhiccontacts_na_removed<-na.omit(allhiccontacts)
# 
# region1hic <- GRanges(allhiccontacts$chr1[1:10],
#                       IRanges(as.numeric(allhiccontacts$start1)[1:10], as.numeric(allhiccontacts$end1)[1:10]))
# region2hic <- GRanges(allhiccontacts$chr2[1:10],
#                       IRanges(as.numeric(allhiccontacts$start2)[1:10], as.numeric(allhiccontacts$end2)[1:10]))
# hicinteraction<-GInteractions(region1hic,region2hic)
# region2hic <- GRanges(allhiccontacts$chr1,ranges =  IRanges(as.numeric(allhiccontacts$start1), as.numeric(allhiccontacts$end1)))
#all.regions <- GRanges("chrA", IRanges(0:9*10+1, 1:10*10))
#index.1 <- c(1,5,10)
#index.2 <- c(3,2,6)
#region.1 <- all.regions[index.1]
#region.2 <- all.regions[index.2]
#gi <- GInteractions(region.1, region.2)

#allhiccontacts_readr<-read_csv("W:/dalgleishjl/straw/tohiccompare/allhiccontacts.csv")
#allhiccontacts_fread<-fread("W:/dalgleishjl/straw/tohiccompare/allhiccontacts.csv")
#head(allhiccontacts)
#class(allhiccontacts_readr)
# makeGenomicInteractionsFromHiCcompare<-function(hiccontactdf,n=nrow(hiccontactdf))
# {
#   region1hic <- GRanges(hiccontactdf$chr1[1:n],
#                         IRanges(as.numeric(hiccontactdf$start1)[1:n], as.numeric(hiccontactdf$end1)[1:n]))
#   region2hic <- GRanges(hiccontactdf$chr2[1:n],
#                         IRanges(as.numeric(hiccontactdf$start2)[1:n], as.numeric(hiccontactdf$end2)[1:n]))
#   hicinteraction<-GInteractions(region1hic,region2hic)
# }
#make_InteractionSet(hic.table = allhiccontacts[1:10,])
#n<-10
#SHORT TEST CODE (enable to validate method), W is equivalent to /data/dalgleishjl/
# allhiccontacts<-read_csv("W:/dalgleishjl/straw/tohiccompare/allhiccontacts.csv")
# hicinteraction_test<-makeGenomicInteractionsFromHiCcompare(hiccontactdf = na.omit(allhiccontacts_readr))

makeGenomicInteractionsFromHiCcompare<-function(hiccontactdf,n=nrow(hiccontactdf))
{
  GInteractions( GRanges(hiccontactdf$chr1[1:n],
                        IRanges(as.numeric(hiccontactdf$start1)[1:n], as.numeric(hiccontactdf$end1)[1:n])),
  GRanges(hiccontactdf$chr2[1:n],
                        IRanges(as.numeric(hiccontactdf$start2)[1:n], as.numeric(hiccontactdf$end2)[1:n])))
}
#From Sean's code
makeGenomicInteractionsFromHiCPro <- function(matrixfile,bedfile) {
  beddf = rtracklayer::import(bedfile)
  matrixdf = read_tsv(matrixfile,col_names=FALSE)
  colnames(matrixdf) = c('bini','binj','count')
  GenomicInteractions(beddf[matrixdf[['bini']]],beddf[matrixdf[['binj']]],matrixdf[['count']])
}

#adding p-values:
GInteractions( GRanges(hiccontactdf$chr1[1:n],
                       IRanges(as.numeric(hiccontactdf$start1)[1:n], as.numeric(hiccontactdf$end1)[1:n])),
               GRanges(hiccontactdf$chr2[1:n],
                       IRanges(as.numeric(hiccontactdf$start2)[1:n], as.numeric(hiccontactdf$end2)[1:n])))
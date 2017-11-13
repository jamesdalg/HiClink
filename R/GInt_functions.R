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
library(InteractionSet)
library(readr)
library(HiCcompare)

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
# GInteractions( GRanges(hiccontactdf$chr1[1:n],
#                        IRanges(as.numeric(hiccontactdf$start1)[1:n], as.numeric(hiccontactdf$end1)[1:n])),
#                GRanges(hiccontactdf$chr2[1:n],
#                        IRanges(as.numeric(hiccontactdf$start2)[1:n], as.numeric(hiccontactdf$end2)[1:n])))



plotSushiArcPlotsWithGenes<-function(chrom,chromstart,chromend,window,txbedpe,controlbedpe,chiptrack1){
  par(mar=c(1,1,1,1))
  par(mfrow=c(4,1))
  txbedpe_window<-txbedpe[txbedpe$chr1==chrom & txbedpe$start1 >= (chromstart-window) & txbedpe$end1 <= (chromend+window),]
  Neg.bedpe_window<-Neg.bedpe[Neg.bedpe$chr1==chrom & Neg.bedpe$start1 >= (chromstart-window) & Neg.bedpe$end1 <= (chromend+window),]
  
  pbpeNeg = plotBedpe(bedpedata = Neg.bedpe_window,chrom = chrom,chromstart = (chromstart-window),chromend = (chromend+window),
                      heights = Neg.bedpe_window$score,plottype="loops",color="red")
  labelgenome(chrom, chromstart,chromend,n=3,scale="Mb")
  
  pbpePg = plotBedpe(bedpedata = txbedpe_window,chrom = chrom,chromstart = (chromstart-window),chromend = (chromend+window),
                     heights = txbedpe_window$score,plottype="loops",color="blue")
  labelgenome(chrom, (chromstart-window),(chromend+window),n=3,scale="Mb")
  legend("topright",inset =0.01,legend=c("T47DPg","T47DNeg"),
         col=SushiColors(3)(3),pch=19,bty='n',text.font=2)
  ensemblgenes_window<-ensemblgenes_patch_removed[ensemblgenes_patch_removed$chrom==chrom & ensemblgenes_patch_removed$start >= (chromstart - window)  & ensemblgenes_patch_removed$end <= (chromend + window),]
  
  pg = plotGenes(ensemblgenes_window,chrom,(chromstart-window),(chromend+window) ,
                 types="exon",maxrows=1,bheight=0.2,
                 plotgenetype="arrow",bentline=FALSE,
                 labeloffset=.4,fontsize=1.2,arrowlength = 0.025,
                 labeltext=TRUE)
  plotBedgraph(signal = chiptrack1,chrom=chrom,chromstart=(chromstart-window),chromend = (chromend+window),transparency = 0.50,color="green")
  
}
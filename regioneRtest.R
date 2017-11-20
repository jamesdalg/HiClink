#library(c(regioneR,BioCStyle,knitr,BSgenome.Hsapiens.UCSC.hg19.masked,BSgenome.Hsapiens.UCSC.hg19,testthat))
library(regioneR)
library(BiocStyle)
library(knitr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
library(testthat)
library(foreach)
library(readr)
library(data.table)
library(devtools)
library(plyr)
library(data.table)

options(error = utils::recover)
#install.packages(c("BiocStyle", "knitr", "BSgenome.Hsapiens.UCSC.hg19.masked", "testthat"))
hg19genmask<-getGenomeAndMask(genome = BSgenome.Hsapiens.UCSC.hg19,mask =BSgenome.Hsapiens.UCSC.hg19.masked )
#need to install the masked version.
#hg19genmask$mas
#get Granges1 & 2, run permutation test
#these tests will probably have to be run on biowulf, ultimately. Use small samples of each GRange.
#production analysis code (not functionalized yet)
allhiccontacts_readr<-read_csv("W:/dalgleishjl/straw/tohiccompare/allhiccontacts_withoutheaders.csv")
hiccontactdf<-na.omit(allhiccontacts_readr)
n<-nrow(hiccontactdf)
hiccontacts_with_metadata<-GenomicInteractions( GRanges(hiccontactdf$chr1[1:n],
                                                        IRanges(as.numeric(hiccontactdf$start1)[1:n], as.numeric(hiccontactdf$end1)[1:n])),
                                                GRanges(hiccontactdf$chr2[1:n],
                                                        IRanges(as.numeric(hiccontactdf$start2)[1:n], as.numeric(hiccontactdf$end2)[1:n])),...=as.data.frame(hiccontactdf[,7:ncol(hiccontactdf)]))
setwd("W:/dalgleishjl/HiClink/HiClink/")
devtools::load_all(".")
#get interaction data
hicgenint_small<-makeGenomicInteractionsFromHiCcompare(na.omit(read_csv("W:/dalgleishjl/straw/tohiccompare/allhiccontacts.csv",n_max=10000)),includemetadata = TRUE)
hicgenint_full<-makeGenomicInteractionsFromHiCcompare(na.omit(read_csv("W:/dalgleishjl/straw/tohiccompare/allhiccontacts_withoutheaders.csv")),includemetadata = TRUE)
hicgenint_full_sig<-hicgenint_full[hicgenint_full@elementMetadata$....p.value<0.05,]
#get chipseq data
SATB1.Pgr.MACS.summits.bed<-data.table::fread("W:/dalgleishjl/chipseq/swarmoutput/T47-D_SATB1_ChIP_Pgr_08.20.09T-47D_SATB1-10.02.09_summits.bed",sep="\t",col.names = c("chrom","chromstart","chromend","name","score"))
SATB1.Pg.MACS.summits.bed<-data.table::fread("W:/dalgleishjl/chipseq/swarmoutput/T47D-HiC-Pro-SATB1-1T47D-HiC-Neg-SATB1-1_summits.bed",sep="\t",col.names = c("chrom","chromstart","chromend","name","score"))
SATB1.Pgr.MACS.summits.bed.GR<-GenomicRanges::GRanges(SATB1.Pgr.MACS.summits.bed)
SATB1.Pg.MACS.summits.bed.GR<-GenomicRanges::GRanges(SATB1.Pg.MACS.summits.bed)
chromosomestokeep<-paste0("chr",c(1:22,"X"))
#the above is needed to fix the error in keepseqlevels... invalid seqlevels: chrM (or chrY) errors in the filterChromosomes function.
SATB1.Pgr.MACS.summits.bed.GR.filtered<-filterChromosomes(SATB1.Pgr.MACS.summits.bed.GR,organism = "hg",chr.type = "canonical",keep.chr =intersect(seqlevels(SATB1.Pgr.MACS.summits.bed.GR),chromosomestokeep))
SATB1.Pg.MACS.summits.bed.GR.filtered<-filterChromosomes(SATB1.Pg.MACS.summits.bed.GR,organism = "hg",chr.type = "canonical",keep.chr =intersect(seqlevels(SATB1.Pg.MACS.summits.bed.GR),chromosomestokeep))
save.image("data_before_permtest.RData")
#need to convert the genome to granges, as well as the mask
BSgenome.Hsapiens.UCSC.hg19
chippermtest<-permTest(A=SATB1.Pgr.MACS.summits.bed.GR.filtered,B=SATB1.Pg.MACS.summits.bed.GR.filtered,ntimes=1000,
                              randomize.function=randomizeRegions,
                              evaluate.function=numOverlaps, count.once=TRUE,
                              genome="hg19", mc.set.seed=FALSE, mc.cores=1)
#divide into genomic bins, test the first anchors, then test the second anchors in each bin.
#anchorOne(hicgenint)
#the hicresults file MUST BE FILTERED OR the permutation test needs to use statistics from HiCcompare as a metric.
#overlap % could also be used-- see overlapRegions(type"binA",get.pctA=T)
firstranges<-anchorOne(hicgenint)
currentchrom<-"chr10"
currentbinstart<-0
currentbinwidth<-10000
firstranges[seqnames(firstranges)==currentchrom,]
SATB1.Prg.chip.perm.test<-permTest(
A=GenomicRanges::restrict(x=SATB1.Pgr.MACS.summits.bed.GR.filtered,start=currentbinstart,end = currentbinstart+currentbinwidth),
B=GenomicRanges::restrict(x=anchorOne(hicgenint_full_sig),start=currentbinstart,end = currentbinstart+currentbinwidth)
,ntimes=100,
                       randomize.function=randomizeRegions,
                       evaluate.function=numOverlaps, count.once=TRUE,
                       genome="hg19", mc.set.seed=FALSE, mc.cores=1)
save.image("databeforepermtestwithchipseqandanchors.RData")
ptm <- proc.time()
SATB1.Prg.chip.perm.test.anchor1.100perm<-permTest(A=SATB1.Pgr.MACS.summits.bed.GR.filtered,B=anchorTwo(hicgenint_full_sig),ntimes=100,
                                                  randomize.function=randomizeRegions,
                                                  evaluate.function=numOverlaps, count.once=TRUE,
                                                  genome="hg19", mc.set.seed=FALSE, mc.cores=1)
proc.time() - ptm
SATB1.Prg.chip.perm.test.anchor2.100perm<-permTest(A=SATB1.Pgr.MACS.summits.bed.GR.filtered,B=anchorTwo(hicgenint_full_sig),ntimes=100,
                                                  randomize.function=randomizeRegions,
                                                  evaluate.function=numOverlaps, count.once=TRUE,
                                                  genome="hg19", mc.set.seed=FALSE, mc.cores=1)
proc.time() - ptm
SATB1.Prg.chip.perm.test.anchor1.1kperm<-permTest(A=SATB1.Pgr.MACS.summits.bed.GR.filtered,B=anchorTwo(hicgenint_full_sig),ntimes=1000,
                                                   randomize.function=randomizeRegions,
                                                   evaluate.function=numOverlaps, count.once=TRUE,
                                                   genome="hg19", mc.set.seed=FALSE, mc.cores=1)
proc.time() - ptm
SATB1.Prg.chip.perm.test.anchor2.1kperm<-permTest(A=SATB1.Pgr.MACS.summits.bed.GR.filtered,B=anchorTwo(hicgenint_full_sig),ntimes=1000,
                                                   randomize.function=randomizeRegions,
                                                   evaluate.function=numOverlaps, count.once=TRUE,
                                                   genome="hg19", mc.set.seed=FALSE, mc.cores=1)
proc.time() - ptm
SATB1.Prg.chip.perm.test.anchor1.10kperm<-permTest(A=SATB1.Pgr.MACS.summits.bed.GR.filtered,B=anchorTwo(hicgenint_full_sig),ntimes=10000,
                                           randomize.function=randomizeRegions,
                                           evaluate.function=numOverlaps, count.once=TRUE,
                                           genome="hg19", mc.set.seed=FALSE, mc.cores=1)

SATB1.Prg.chip.perm.test.anchor2.10kperm<-permTest(A=SATB1.Pgr.MACS.summits.bed.GR.filtered,B=anchorTwo(hicgenint_full_sig),ntimes=10000,
                       randomize.function=randomizeRegions,
                       evaluate.function=numOverlaps, count.once=TRUE,
                       genome="hg19", mc.set.seed=FALSE, mc.cores=1)
hicgenint_full_sig_sample<-sample(hicgenint_full_sig,5000)
getBinwiseOverlapPercentsForGenome<-function(hiccontacts_with_metadata,chipseqGRanges,binwidth)
{
  #iterate across a list of chromosomes (for loop)
 #iterate across a chromosome (for loop)

  #what to do about the fact that there are two contacts?
  #One could simply throw them in together, but are first contacts different than second contacts?
  #conversely, there could be data where patterns cannot be seen without looking at both.
  #Doing both will let the user decide.
  #do both contacts separately, just anchorOne,just anchorTwo, or both together as a combined GRange object.
  
  #put all of this as an output in a compound S4 object at the end
  
}
#make first a small function that does binwise across a single chromosome, across a single set of anchors, defaulting to anchorOne.
generateHiCbinsubset<-function(hiccontacts_with_metadata,chrom=NULL,start=0,binwidth=NULL,anchor="first")
{
  GenomicRanges::restrict(x=firstanchor[seqnames(InteractionSet::anchors(hiccontacts_with_metadata)[,anchor])==chrom,],start=start,end = start+binwidth)
}
getSignificantInteractions<-function(hiccontacts_with_metadata,p.value_threshold=0.05)
{
  hiccontacts_with_metadata[hiccontacts_with_metadata@elementMetadata$....p.value<0.05,]
}
chipseqGRanges<-SATB1.Pgr.MACS.summits.bed.GR.filtered

#library(AneuFinder)
#install.packages("AneuFinder")
chipseqGRanges<-sample(SATB1.Pgr.MACS.summits.bed.GR.filtered,100)
hicgenint_full_sig_sample<-sample(hicgenint_full_sig,100)
hiccontacts_with_metadata<-hicgenint_full_sig_sample
getBinwiseOverlapStats<-function(hiccontacts_with_metadata,chipseqGRanges,binwidth,chroms="ALL")
{
  #testcase #1
#  binstart<-0
#  binwidth<-1e7
 # chrom<-"chr11"
  #end test data
  #a better way to access anchors: anchors(genomicInteractions)$first
  #we will only consider intrachromosomal interactions for the purposes of this function. It would be simple to rearrange things to allow for intrachromosomal (i.e. where anchor1's chromosome is not equal to anchor2's chromosome)
  if (chroms=="ALL")
  {  
  chromset<-intersect(intersect(seqnames(anchorOne(hiccontacts_with_metadata)),seqnames(chipseqGRanges)),intersect(seqnames(anchorTwo(hiccontacts_with_metadata)),seqnames(chipseqGRanges)))
  } else
  {
    chromset<-intersect(chroms,
    intersect(intersect(seqnames(anchorOne(hiccontacts_with_metadata)),seqnames(chipseqGRanges)),intersect(seqnames(anchorTwo(hiccontacts_with_metadata)),seqnames(chipseqGRanges))))
    
  }
  
  genomestats<-foreach(chromindex=1:length(chromset)) %dopar%
{
  print(chromset[chromindex])
  chrom<-chromset[chromindex]
  anchoroneranges_in_chrom<-anchors(hiccontacts_with_metadata)$first[seqnames(anchors(hiccontacts_with_metadata)$first)==chrom,]
  anchortworanges_in_chrom<-anchors(hiccontacts_with_metadata)$second[seqnames(anchors(hiccontacts_with_metadata)$second)==chrom,]
#  anchortworanges<-ranges(anchorOne(hiccontacts_with_metadata))
  #anchoroneranges_in_chrom<-anchoroneranges[seqnames(hiccontacts_with_metadata@elementMetadata),]
  #anchortworanges_in_chrom<-ranges(anchorOne(hiccontacts_with_metadata))
  chipseqGRanges[seqnames(chipseqGRanges)==chrom,]
  maxpos<-max(max(anchoroneranges_in_chrom@ranges@start+anchoroneranges_in_chrom@ranges@width),max(anchortworanges_in_chrom@ranges@start+anchortworanges_in_chrom@ranges@width))
  maxbinstart<-floor(maxpos/binwidth)*binwidth
  maxbinend<-ceiling(maxpos/binwidth)*binwidth
  #get the size of the chromosome (or the maximum position of the end of the range, rounded down to the nearest binwidth for the start)
  binstarts<-c(seq(from=0, to=maxbinstart,by=binwidth))
#it only gathers stats for those regions that have overlaps. Others throw an error and pass handling takes over.
  chromstats<-foreach::foreach(binindex=1:(length(binstarts)-1),.combine=rbind,.errorhandling = "pass") %dopar%
  {
    print(binstarts[binindex])
  firstanchor<-anchorOne(hiccontacts_with_metadata)
  secondanchor<-anchorOne(hiccontacts_with_metadata)
  hiccontacts_with_metadata_anchorOne_currentbin_subset<-GenomicRanges::restrict(x=firstanchor[seqnames(firstanchor)==chrom,],start=binstarts[binindex],end = binstarts[binindex]+binwidth)
  hiccontacts_with_metadata_anchorTwo_currentbin_subset<-GenomicRanges::restrict(x=secondanchor[seqnames(secondanchor)==chrom,],start=binstarts[binindex],end = binstarts[binindex]+binwidth)
  chipseqGRanges_currentbin_subset<-GenomicRanges::restrict(x=chipseqGRanges[seqnames(chipseqGRanges)==chrom,],start=binstarts[binindex],end = binstarts[binindex]+binwidth)
  
  anchorOneOverlaps<-regioneR::overlapRegions(hiccontacts_with_metadata_anchorOne_currentbin_subset,chipseqGRanges_currentbin_subset,get.pctA = T,get.pctB = T,get.bases = T,type="any")
  anchorTwoOverlaps<-regioneR::overlapRegions(hiccontacts_with_metadata_anchorTwo_currentbin_subset,chipseqGRanges_currentbin_subset,get.pctA = T,get.pctB = T,get.bases = T,type="any")
 #regioneR::overlapRegions(A=hiccontacts_with_metadata_bin_subset,B=chipseqGRanges,type="any",get.bases=T) 
  
  #tabulate all overlaps for both anchors and combine
  labeledAnchorOneOverlaps<-cbind(anchorOneOverlaps,rep(1,nrow(anchorOneOverlaps)))
  names(labeledAnchorOneOverlaps)[length(labeledAnchorOneOverlaps)]<-"anchor"
  labeledAnchorTwoOverlaps<-cbind(anchorTwoOverlaps,rep(2,nrow(anchorTwoOverlaps)))
  names(labeledAnchorTwoOverlaps)[length(labeledAnchorTwoOverlaps)]<-"anchor"
  allBinOverlaps<-rbind(labeledAnchorOneOverlaps,labeledAnchorTwoOverlaps)
  # if (nrow(allBinOverlaps)==0)
  # #{
  # #  Instead of adding a stats line with NA for a bin that has no overlaps, the program bypasses this region and moves forward.
  # 
  # }
  if (nrow(allBinOverlaps)!=0)
  {
  #colMeans(allBinOverlaps)
  bin_jaccard<-nrow(allBinOverlaps)/(nrow(hiccontacts_with_metadata_anchorOne_currentbin_subset)+nrow(hiccontacts_with_metadata_anchorOne_currentbin_subset))
  binstats<-c(chrom,binstarts[binindex],binstarts[binindex]+binwidth,binwidth,matrixStats::colMeans2(as.matrix(allBinOverlaps[,c("ov.bases","pct.basesA","pct.basesB")])),bin_jaccard)
  names(binstats)<-c("chrom","binstart","binend","binwidth","ov.bases.mean","pct.overlappingbasesHiC.mean","pct.overlappingbaseschipseq.mean","bin_jaccard")
  binstats
  } else {
    #stop(paste0("no bin overlaps in region",chrom,":",binstarts[binindex],"-",binstarts[binindex]+binwidth))
    binstats<-NULL
  }
    
}  #need to calculate jaccard index for the bin.
}  
#if (vebose==true)
# output<-list("overlaps"=allBinOverlaps)
#  output$overlaps
#  
}
registerDoMC(8)
ptm <- proc.time()
genomestats_100sample<-getBinwiseOverlapStats(hiccontacts_with_metadata = sample(hicgenint_full_sig,100),chipseqGRanges = sample(SATB1.Pgr.MACS.summits.bed.GR.filtered,100),binwidth=1e7)
proc.time() - ptm
genomestats_5000sample<-getBinwiseOverlapStats(hiccontacts_with_metadata = sample(hicgenint_full_sig,5000),chipseqGRanges = sample(SATB1.Pgr.MACS.summits.bed.GR.filtered,5000),binwidth=1e7)
save.image("overlapfunction_genome_wide_complete.RData")
genomestats_1e7complete<-getBinwiseOverlapStats(hiccontacts_with_metadata = hicgenint_full_sig,chipseqGRanges = SATB1.Pgr.MACS.summits.bed.GR.filtered,binwidth=1e7)
genomestats_1e7complete_df<-ldply(genomestats_1e7complete,data.frame)
library(plyr)
genomestats_1e7complete<-getBinwiseOverlapStats(hiccontacts_with_metadata = hicgenint_full_sig,chipseqGRanges = SATB1.Pgr.MACS.summits.bed.GR.filtered,binwidth=1e7)
genomestats_1e7complete_df<-ldply(genomestats_1e7complete,data.frame)
fwrite(genomestats_1e7complete_df,"genomestats_1e7complete_SATB1.Pgr_Pg_Neg_HiCcompare_0.05sig.csv")
genomestats_1e6complete<-getBinwiseOverlapStats(hiccontacts_with_metadata = hicgenint_full_sig,chipseqGRanges = SATB1.Pgr.MACS.summits.bed.GR.filtered,binwidth=1e6)
genomestats_1e6complete_df<-ldply(genomestats_1e6complete,data.frame)
fwrite(genomestats_1e6complete_df,"genomestats_1e6complete_SATB1.Pgr_Pg_Neg_HiCcompare_0.05sig.csv")
genomestats_1e5complete<-ldply(getBinwiseOverlapStats(hiccontacts_with_metadata = hicgenint_full_sig,chipseqGRanges = SATB1.Pgr.MACS.summits.bed.GR.filtered,binwidth=1e5),data.frame)
fwrite(genomestats_1e5complete_df,"genomestats_1e5complete_SATB1.Pgr_Pg_Neg_HiCcompare_0.05sig.csv")
genomestats_1e4complete<-ldply(getBinwiseOverlapStats(hiccontacts_with_metadata = hicgenint_full_sig,chipseqGRanges = SATB1.Pgr.MACS.summits.bed.GR.filtered,binwidth=1e4),data.frame)
fwrite(genomestats_1e4complete_df,"genomestats_1e5complete_SATB1.Pgr_Pg_Neg_HiCcompare_0.05sig.csv")
genomestats_1e3complete<-ldply(getBinwiseOverlapStats(hiccontacts_with_metadata = hicgenint_full_sig,chipseqGRanges = SATB1.Pgr.MACS.summits.bed.GR.filtered,binwidth=1e3),data.frame)
fwrite(genomestats_1e3complete_df,"genomestats_1e5complete_SATB1.Pgr_Pg_Neg_HiCcompare_0.05sig.csv")

#-----------------------------------------------------------------------------------------------------
#EXAMPLE1 code for regioneR (p.29)


set.seed(12345)
cpgHMM <- toGRanges("http://www.haowulab.org/software/makeCGI/model-based-cpg-islands-hg19.txt")
promoters <- toGRanges("http://gattaca.imppc.org/regioner/data/UCSC.promoters.hg19.bed")
cpgHMM_filtered <- filterChromosomes(cpgHMM, organism="hg", chr.type="canonical")
promoters_filtered <- filterChromosomes(promoters, organism="hg", chr.type="canonical")
cpgHMM_2K <- sample(cpgHMM_filtered, 2000)
pt <- overlapPermTest(cpgHMM_2K, promoters_filtered, ntimes=1000, genome="hg19", count.once=TRUE)
pt
mean(pt$permuted)
plot(pt) 

#EXAMPLE2 (final, p.32)
set.seed(12345)
download.file("http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsSydhHepg2Rad21IggrabUniPk.narrowPeak.gz", "Rad21.gz")
download.file("http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsUwHepg2CtcfUniPk.narrowPeak.gz", "Ctcf.gz")
HepG2_Rad21 <- toGRanges(gzfile("Rad21.gz"), header=FALSE)
HepG2_Ctcf <- toGRanges(gzfile("Ctcf.gz"), header=FALSE)
promoters <- toGRanges("http://gattaca.imppc.org/regioner/data/UCSC.promoters.hg19.bed")
promoters <- filterChromosomes(promoters, organism="hg19")
HepG2_Rad21_5K <- sample(HepG2_Rad21, 5000)
pt_Rad21_5k_vs_Ctcf <- permTest(A=HepG2_Rad21_5K, B=HepG2_Ctcf, ntimes=1000,
                                randomize.function=circularRandomizeRegions,
                                evaluate.function=numOverlaps, count.once=TRUE,
                                genome="hg19", mc.set.seed=FALSE, mc.cores=1)
pt_Rad21_5k_vs_Prom <- permTest(A=HepG2_Rad21_5K, B=promoters, ntimes=1000,
                                randomize.function=circularRandomizeRegions,
                                evaluate.function=numOverlaps, count.once=TRUE,
                                genome="hg19", mc.set.seed=FALSE, mc.cores=1)
pt_Rad21_5k_vs_Ctcf
pt_Rad21_5k_vs_Prom
plot(pt_Rad21_5k_vs_Ctcf, main="Rad21_5K vs CTCF")
plot(pt_Rad21_5k_vs_Prom, main="Rad21_5K vs Promoters")

lz_Rad21_vs_Ctcf_1 <- localZScore(A=HepG2_Rad21_5K, B=HepG2_Ctcf,
                                  pt=pt_Rad21_5k_vs_Ctcf,
                                  window=1000, step=50, count.once=TRUE)
lz_Rad21_vs_Prom_1 <- localZScore(A=HepG2_Rad21_5K, B=promoters,
                                  pt=pt_Rad21_5k_vs_Prom,
                                  window=1000, step=50, count.once=TRUE)
plot(lz_Rad21_vs_Ctcf_1, main="Rad21_5k vs CTCF (1Kbp)")
plot(lz_Rad21_vs_Prom_1, main="Rad21_5k vs promoters (1Kbp)")

#NOT RUN
lz_Rad21_vs_Ctcf_2 <- localZScore(A=HepG2_Rad21_5K, B=HepG2_Ctcf,
                                  pt=pt_Rad21_5k_vs_Ctcf,
                                  window=10000, step=500, count.once=TRUE)
lz_Rad21_vs_Prom_2 <- localZScore(A=HepG2_Rad21_5K, B=promoters,
                                  pt=pt_Rad21_5k_vs_Prom,
                                  window=10000, step=500, count.once=TRUE)
plot(lz_Rad21_vs_Ctcf_2, main="Rad21_5k vs CTCF (10Kbp)")
plot(lz_Rad21_vs_Prom_2, main="Rad21_5k vs promoters (10Kbp)")
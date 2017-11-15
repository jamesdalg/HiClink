#library(c(regioneR,BioCStyle,knitr,BSgenome.Hsapiens.UCSC.hg19.masked,BSgenome.Hsapiens.UCSC.hg19,testthat))
library(regioneR)
library(BiocStyle)
library(knitr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
library(testthat)
#install.packages(c("BiocStyle", "knitr", "BSgenome.Hsapiens.UCSC.hg19.masked", "testthat"))
hg19genmask<-getGenomeAndMask(genome = BSgenome.Hsapiens.UCSC.hg19,mask =BSgenome.Hsapiens.UCSC.hg19.masked )
#need to install the masked version.
hg19genmask$mas
#get Granges1 & 2, run permutation test
#these tests will probably have to be run on biowulf, ultimately. Use small samples of each GRange.
hiccontactdf<-na.omit(allhiccontacts_readr)
n<-nrow(hiccontactdf)
hiccontacts_with_metadata<-GenomicInteractions( GRanges(hiccontactdf$chr1[1:n],
                                                        IRanges(as.numeric(hiccontactdf$start1)[1:n], as.numeric(hiccontactdf$end1)[1:n])),
                                                GRanges(hiccontactdf$chr2[1:n],
                                                        IRanges(as.numeric(hiccontactdf$start2)[1:n], as.numeric(hiccontactdf$end2)[1:n])),...=as.data.frame(hiccontactdf[,7:ncol(hiccontactdf)]))


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
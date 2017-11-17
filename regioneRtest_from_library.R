library(HiClink)

#test case #2-- multiple chromosomes
guessmax<-as.numeric(strsplit(system(paste("wc -l","/data/CCRBioinfo/dalgleishjl/straw/tohiccompare/allhiccontacts.csv"),intern=T)," ")[[1]][1])
#strsplit(guessmax," ")[[1]][1]
full_hiccontacts_table_without_headers<-na.omit(read_csv("/data/CCRBioinfo/dalgleishjl/straw/tohiccompare/allhiccontacts_withoutheaders.csv"))
full_hiccontacts_table_with_headers<-na.omit(read_csv("/data/CCRBioinfo/dalgleishjl/straw/tohiccompare/archive/allhiccontacts.csv"))
library(compare)
#comparison<-compare(as.data.frame(full_hiccontacts_table),as.data.frame(full_hiccontacts_table_with_headers))
hicgenint_full<-makeGenomicInteractionsFromHiCcompare(na.omit(read_csv("/data/CCRBioinfo/dalgleishjl/straw/tohiccompare/allhiccontacts_withoutheaders.csv")),includemetadata = TRUE)
hicgenint_full_sig<-hicgenint_full[hicgenint_full@elementMetadata$....p.value<0.05,]
SATB1.Pgr.MACS.summits.bed<-data.table::fread("/data/CCRBioinfo/dalgleishjl/chipseq/swarmoutput/T47-D_SATB1_ChIP_Pgr_08.20.09T-47D_SATB1-10.02.09_summits.bed",sep="\t",col.names = c("chrom","chromstart","chromend","name","score"))

SATB1.Pgr.MACS.summits.bed.GR<-GenomicRanges::GRanges(SATB1.Pgr.MACS.summits.bed)
chromosomestokeep<-paste0("chr",c(1:22,"X"))
SATB1.Pgr.MACS.summits.bed.GR.filtered<-filterChromosomes(SATB1.Pgr.MACS.summits.bed.GR,organism = "hg",chr.type = "canonical",keep.chr =intersect(seqlevels(SATB1.Pgr.MACS.summits.bed.GR),chromosomestokeep))
registerDoMC()
library(SNOW)
genomestats_1e7complete<-getBinwiseOverlapStats(hiccontacts_with_metadata = hicgenint_full_sig,chipseqGRanges = SATB1.Pgr.MACS.summits.bed.GR.filtered,binwidth=1e7)
genomestats_1e7complete_df<-ldply(genomestats_1e7complete,data.frame)
fwrite(genomestats_1e7complete_df,"genomestats_1e7complete_SATB1.Pgr_Pg_Neg_HiCcompare_0.05sig.csv")
#additional analysis code to run:
genomestats_1e6complete<-getBinwiseOverlapStats(hiccontacts_with_metadata = hicgenint_full_sig,chipseqGRanges = SATB1.Pgr.MACS.summits.bed.GR.filtered,binwidth=1e6)
genomestats_1e6complete_df<-ldply(genomestats_1e6complete,data.frame)
fwrite(genomestats_1e6complete_df,"genomestats_1e6complete_SATB1.Pgr_Pg_Neg_HiCcompare_0.05sig.csv")
genomestats_1e5complete<-ldply(getBinwiseOverlapStats(hiccontacts_with_metadata = hicgenint_full_sig,chipseqGRanges = SATB1.Pgr.MACS.summits.bed.GR.filtered,binwidth=1e5),data.frame)
fwrite(genomestats_1e5complete_df,"genomestats_1e5complete_SATB1.Pgr_Pg_Neg_HiCcompare_0.05sig.csv")
genomestats_1e4complete<-ldply(getBinwiseOverlapStats(hiccontacts_with_metadata = hicgenint_full_sig,chipseqGRanges = SATB1.Pgr.MACS.summits.bed.GR.filtered,binwidth=1e4),data.frame)
fwrite(genomestats_1e4complete_df,"genomestats_1e5complete_SATB1.Pgr_Pg_Neg_HiCcompare_0.05sig.csv")
genomestats_5e3complete<-ldply(getBinwiseOverlapStats(hiccontacts_with_metadata = hicgenint_full_sig,chipseqGRanges = SATB1.Pgr.MACS.summits.bed.GR.filtered,binwidth=5e3),data.frame)
fwrite(genomestats_1e4complete_df,"genomestats_1e5complete_SATB1.Pgr_Pg_Neg_HiCcompare_0.05sig.csv")

genomestats_1e3complete<-ldply(getBinwiseOverlapStats(hiccontacts_with_metadata = hicgenint_full_sig,chipseqGRanges = SATB1.Pgr.MACS.summits.bed.GR.filtered,binwidth=1e3),data.frame)
fwrite(genomestats_1e3complete_df,"genomestats_1e5complete_SATB1.Pgr_Pg_Neg_HiCcompare_0.05sig.csv")

#end analysis code

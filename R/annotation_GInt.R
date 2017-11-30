#reference: https://github.com/kbroman/pkg_primer/blob/gh-pages/example/stage3/R/plot_crayons.R
#' Annotates HiCcompare interaction data output.
#' 
#' @param filename filename of the HiCcompare output, in csv format.If NULL, then the provided dataframe is used instead.
#'
#' @param df dataframe containing HiCcompare output, containing fields "chr1"        "start1"      "end1"        "chr2"        "start2"      "end2"        "IF1"         "IF2"         "D"  "M"           "adj.IF1"     "adj.IF2"     "adj.M"       "mc"          "A"           "p.value"     "fold.change"
#' @param is Interactionset containing the above. Please note: only one of the three options may be used.
#' @return annotatedInteractionSet An annotated interaction set, containing HGNC symbols
#' @author James L. T. Dalgleish
#' @seealso \code{\link[GIChipSeqCompare]}
#' @keywords chipseq HiC HiChip GInteractions
#' 
#' @examples 
#' afunction()
#' 
#' @export
#' @importFrom readr read_csv
#' @importFrom GenomicInteractions
#' 
#' 
library(data.table)
library(GenomicRanges)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicFeatures)
csv_input<-data.table::fread("W:/dalgleishjl/hiclink/HiClink/output/Copy of genomestats_1e4complete_SATB1.Pgr_Pg_Neg_HiCcompare_0.05sig.csv")
tenkbhicbins<-GenomicRanges::makeGRangesFromDataFrame(csv_input,keep.extra.columns = T)
hg19.refseq.db<-GenomicFeatures::makeTxDbFromUCSC(genome="hg19", table="refGene")
refseq.genes = genes(hg19.refseq.db)
refseq.transcripts = transcriptsBy(hg19.refseq.db, by="gene")
hg19_refseq_promoters <- unique(unlist(promoters(refseq.transcripts, 2500,2500)))
mart = biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
genes <- biomaRt::getBM(attributes = c("hgnc_symbol", "refseq_mrna"), filter = "refseq_mrna",
                        values = hg19_refseq_promoters$tx_name, mart = mart)
genetxtable<-na.omit(genes[match(hg19_refseq_promoters$tx_name, genes$refseq_mrna),])
hg19_refseq_promoters$geneSymbol <- genes$hgnc_symbol[match(hg19_refseq_promoters$tx_name, genes$refseq_mrna)]
names(hg19_refseq_promoters)<-hg19_refseq_promoters$geneSymbol
na.symbol <- is.na(names(hg19_refseq_promoters))
names(hg19_refseq_promoters)[na.symbol] <- hg19_refseq_promoters$tx_name[na.symbol]

save.image("annotation_GInt.RData")
overlaps<-GenomicRanges::findOverlaps(hg19_refseq_promoters,tenkbhicbins)
#overlaps@to
mcols(tenkbhicbins)$geneSymbol<-unlist(lapply(mcols(hg19_refseq_promoters)$geneSymbol,paste,collapse=","))

#tenkbhicbins_inpromoter<-tenkbhicbins %over% promoters(hg19.refseq.db,upstream=1.5e5,downstream=2e3)
tenkbhicbins_inpromoter2<-overlapsAny(tenkbhicbins,promoters(hg19.refseq.db))
promoters_intenkbhicbins<-overlapsAny(promoters(hg19.refseq.db),tenkbhicbins)
tenkbhicbins_inpromoter_subset<-tenkbhicbins[tenkbhicbins_inpromoter,]
tenkbhicbin_hg19promoter_subset<-promoters(hg19.refseq.db)[promoters_intenkbhicbins,]
#mcols(inpromoter_subset)$geneName<-
tenkbhicbin_hg19promoter_subset
#table(tenkbhicbins_inpromoter == tenkbhicbins_inpromoter2) #they are identical.
#small use case: 
promoters(hg19.refseq.db)$tx_name[subjectHits(findOverlaps(tenkbhicbins, promoters(hg19.refseq.db),maxgap=1e5))]

#creating a dataframe:
hits_with_txnames<-cbind(as.data.frame(findOverlaps(tenkbhicbins, promoters(hg19.refseq.db),maxgap=1e5)),promoters(hg19.refseq.db)$tx_name[subjectHits(findOverlaps(tenkbhicbins, promoters(hg19.refseq.db),maxgap=1e5))])
hits_with_genenames<-cbind(
  as.data.frame(findOverlaps(tenkbhicbins, hg19_refseq_promoters,maxgap=1e5)),
  hg19_refseq_promoters$geneSymbol[subjectHits(findOverlaps(tenkbhicbins, hg19_refseq_promoters,maxgap=1e5))],
  hg19_refseq_promoters$tx_name[subjectHits(findOverlaps(tenkbhicbins, hg19_refseq_promoters,maxgap=1e5))])

#need to collapse down the txnames.
names(hits_with_txnames)[3]<-"txname"
names(hits_with_genenames)[3]<-"geneSymbol"
names(hits_with_genenames)[4]<-"txname"
library(plyr)
hits_with_txnames_collapsed<-ddply(hits_with_txnames, .(queryHits), summarize, txname=paste(txname, collapse="|"))
mcols(tenkbhicbins)$txname[hits_with_txnames_collapsed$queryHits]<-hits_with_txnames_collapsed$txname
hits_with_genenames_collapsed<-ddply(hits_with_genenames, .(queryHits), summarize, geneSymbol=paste(geneSymbol, collapse="|"))
mcols(tenkbhicbins)$geneSymbol[hits_with_genenames_collapsed$queryHits]<-hits_with_genenames_collapsed$geneSymbol

save.image("working_annotations_with_txname_and_geneSymbol.RData")
#tenkbhicbins[queryHits(findOverlaps(tenkbhicbins, promoters(hg19.refseq.db),maxgap=1e5))]
#tenkbhicbins[queryHits(findOverlaps(tenkbhicbins, promoters(hg19.refseq.db),maxgap=1e5))]


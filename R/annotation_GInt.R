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
#' 
#' 
#' 
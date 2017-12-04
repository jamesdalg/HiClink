library(data.table)
library(readr)
#NegallValidPairs<-fread("Z:/hicpro/fanHiC10112017/complete/output_sbatch_multithreaded_lscratch_new_params/hic_results/data/T47D-HiChip-Neg/T47D-HiChip-Neg_allValidPairs")
#PgallValidPairs<-fread("Z:/hicpro/fanHiC10112017/complete/output_sbatch_multithreaded_lscratch_new_params/hic_results/data/T47D-HiChip-Pg/T47D-HiChip-Pg_allValidPairs")
setwd("/data/CCRBioinfo/dalgleishjl/quasar/")

#fwrite("T47D-HiChip-Pg_allValidPairs.hifive_write.table.raw",x=PgallValidPairs[,2:7],sep="\t",row.names = F,col.names = F)
#fwrite("T47D-HiChip-Neg_allValidPairs.hifive_write.table.raw",x=NegallValidPairs[,2:7],sep="\t",row.names = F,col.names = F)

convertValidPairsToHiFiveRaw<-function(validpairsfilename)
{
 Validpairsdt<-fread(validpairsfilename) 
 fwrite(paste0(validpairsfilename,".hifive.raw"),x=Validpairsdt[,2:7],sep="\t",row.names = F,col.names = F)
 #print(System(paste0("ls -l ",paste0(validpairsfilename,".hifive.raw"))))
}
#test case #2 literature
convertValidPairsToHiFiveRaw("/data/CCRBioinfo/dalgleishjl/hicpro/output_sbatch_literature_hichip/hic_results/data/SRR3467175/SRR3467175_allValidPairs")
convertValidPairsToHiFiveRaw("/data/CCRBioinfo/dalgleishjl/hicpro/output_sbatch_literature_hichip/hic_results/data/SRR3467176/SRR3467176_allValidPairs")

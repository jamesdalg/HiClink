#http://wiki.wubrowse.org/Long-range
#from hiccompare data, we could simply take a single score and add it.
#the format is below:
#Each line indicates an interaction event involving two regions from the genome. 3 fields are present for each line: coordinate of region 1, coordinate of region 2, score.
#
#
#
library(readr)
library(foreach)
hicgenint_full<-read_csv("W:/dalgleishjl/straw/tohiccompare/allhiccontacts_withoutheaders.csv",n_max=5000)
hicgenint_full_sig<-hicgenint_full[hicgenint_full$p.value<0.05,]
#
#-log(hicgenint_full_sig$p.value)
#
#hicgenint_full_sig$
  
  convertToWubrowse<-function(hiccontactsdataframe=NULL,pvaluecutoff=0.05,format="HiCcompare",score=foldchange,trackformat=T)
  {
    if(format=="HiCcompare")
    {
      hiccontactsdataframe_sig<-hiccontactsdataframe[hiccontactsdataframe$p.value<pvaluecutoff,]
      i<-1
      outputdf<-foreach(i=1:nrow(hiccontactsdataframe_sig),.combine="rbind") %do%
      {
        if(trackformat==F)
          {
        outputrow<-cbind(
          paste(c(hiccontactsdataframe_sig$chr1[i],hiccontactsdataframe_sig$start1[i],hiccontactsdataframe_sig$end1[i]),collapse=","),
          paste(c(hiccontactsdataframe_sig$chr2[i],hiccontactsdataframe_sig$start2[i],hiccontactsdataframe_sig$end2[i]),collapse=","),
          hiccontactsdataframe_sig$fold.change[i])
        outputrow
        }

        # if(trackformat==T)
        # {
        #   # chr1   \t   111   \t   222   \t   chr2:333-444,55   \t   1   \t   .
        #   # chr2   \t   333   \t   444   \t   chr1:111-222,55   \t   2   \t   .
        #   
        #   outputrowpair<-rbind(c(hiccontactsdataframe_sig$chr1[i],hiccontactsdataframe_sig$start1[i],hiccontactsdataframe_sig$end1[i],
        #                paste0(hiccontactsdataframe_sig$chr2[i],":",hiccontactsdataframe_sig$start2[i],"-",
        #                       hiccontactsdataframe_sig$end2[i],",",hiccontactsdataframe_sig$fold.change[i])
        #                ,i,"."),
        #                c(hiccontactsdataframe_sig$chr2[i],hiccontactsdataframe_sig$start2[i],hiccontactsdataframe_sig$end2[i],
        #                  paste0(hiccontactsdataframe_sig$chr1[i],":",hiccontactsdataframe_sig$start1[i],"-",
        #                         hiccontactsdataframe_sig$end1[i],",",hiccontactsdataframe_sig$fold.change[i])
        #                  ,nrow(hiccontactsdataframe_sig)+i,"."))
        # outputrowpair
        #   }
      }
    }
  }
test_wubrowse_output<-convertToWubrowse(hiccontactsdataframe = hicgenint_full)
#library(data.table)
setwd("W:/dalgleishjl/HiClink/HiClink/")
write.table(test_wubrowse_output,paste0(getwd(),"/output/test_wubrowse_output.tsv"),sep="\t",quote=F,row.names=F,col.names = F)
nontrackformat_output<-convertToWubrowse(hiccontactsdataframe = hicgenint_full,trackformat = F)
write.table(outputdf,"test_wubrowse_nontrack.txt",sep="\t",quote = F,row.names = F,col.names = F)

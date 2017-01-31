library(RTCGAToolbox)

cases=getFirehoseDatasets()  #To get all the TCGA cases

gisticDate = getFirehoseAnalyzeDates(last=3)

cases_filtered=vector()
TCGA_data=list()
for ( i in seq(1,length(cases))) if(nchar(cases[i]) <= 4) {cases_filtered[i]=cases[i]}
cases_filtered=cases_filtered[!is.na(cases_filtered)]
case_filtered_lower=tolower(cases_filtered)
#TCGA_Mutation_Frequency=list()
for (i in 1:35)
{
  TCGA_data[[case_filtered_lower[i]]]=getFirehoseData(dataset=cases_filtered[i], runDate=gisticDate[1], gistic2_Date=NULL,Clinic=FALSE, RNAseq_Gene=FALSE, mRNA_Array=FALSE, Mutation=TRUE);
  #TCGA_Mutation_Frequency[[paste(case_filtered_lower[i],'clinical',sep="_")]]=getData(TCGA_data[[i]], "GISTIC")
  #write.table(TCGA_Mutation_Frequency[[i]],file=paste(case_filtered_lower[i],"tsv",sep="."),row.names = FALSE)
  #TCGA_data[[i]]=NULL
    
}


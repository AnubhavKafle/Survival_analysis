library(survival)

plotting = function(Cancer_Expression, clinical_data)
{
  event = c() #To have all the events coded in 0 or 1
  
  days_to_death = c() #days to death, if not, then follow-up data
  
  expression_rna = t(apply(Cancer_Expression,1,as.numeric)) #TO eliminate factor problems
  
  total_gene_list = colnames(Cancer_Expression) #Obtaining all the genes in the sample
  
  rownames(expression_rna) <- as.vector(sapply(rownames(Cancer_Expression), function(x) tolower(paste(unlist(strsplit(x,"\\."))[[1]],unlist(strsplit(x,"\\."))[[2]],unlist(strsplit(x,"\\."))[[3]],sep="."))))#To keep all the rows in small letters
  
  colnames(expression_rna) = colnames(Cancer_Expression)  ##############
  
  expression_rna_countMorethanZero = expression_rna[,which(apply(expression_rna, 2, median) >= 0.1)]  ####to filter out gene that have counts less than 0.1 in 50 or more of the sample.
  
  expression_rna_log2 = log2(t(apply(t(apply(expression_rna_countMorethanZero,1,as.numeric)), 1, function(i) i+1 ))) ### Log2 transform of the filtered gene counts
  
  colnames(expression_rna_log2) = colnames(expression_rna_countMorethanZero) #################
  
  match_index = which(rownames(clinical_data)%in% rownames(expression_rna_log2)) #Check which patient ID are matching with the tumor data
  
  clinical_match=as.data.frame(clinical_data[match_index,]) #filter the clinical match. There are Normal patients as well so less mathches
  
  for ( i in 1:length(as.character(clinical_match$days_to_death))) #Take follow-day where death days are not available
  {
    days_to_death[i] = ifelse(is.na(as.character(clinical_match$days_to_death)[i]),as.numeric(as.character(clinical_match$days_to_last_followup))[i],as.numeric(as.character(clinical_match$days_to_death))[i])
    
    event[i] = as.numeric(as.character(clinical_match$vital_status))[i]
  }
  
  rm(clinical_data)
  
  clinical_survival = as.data.frame(cbind(days_to_death,event)) #Create a data frame for survival fit analysis
  
  rownames(clinical_survival)= rownames(clinical_match)#Matching rownames
  
  expression_rna_ordered_log2_median = expression_rna_log2[rownames(clinical_survival),] #Creating expression matrix ordered according to rownames of clinical_match
  
  rm(expression_rna,expression_rna_countMorethanZero,match_index,expression_rna_log2,event, days_to_death,clinical_match)
  
  plot(survfit(Surv(clinical_survival$days_to_death,clinical_survival$event) ~ as.character(expression_rna_ordered_log2_median[,grep(paste('^',genes,sep=""), colnames(expression_rna_ordered_log2_median))] > median(expression_rna_ordered_log2_median[,grep(paste('^',genes,sep=""), colnames(expression_rna_ordered_log2_median))]))),col= c("blue","red"), main=paste(cancer,"Cancer Specific Survival Plot"))
  
  mtext(paste("(p=",signif(p_val, 2),")"),3)
  
  legend(x="bottomleft",legend=paste(genes,c("ct high","ct low")),col=c("red","blue"),lwd=1)
  
  return(NULL)
  
}
  
  Cancer_Expression = t(read.table(list.files("~/Desktop/labRotation1_AnubhavK/gene_expression/",pattern=cancer,full.names = T),header=T,sep="\t",row.names=1))[,-1]
  
  clinical_data = t(read.table(list.files("~/Desktop/labRotation1_AnubhavK/Clinical_Firehose/",pattern=cancer,full.names = T),header=T,sep="\t",row.names = 1))[,-1]
  
  plotting(Cancer_Expression, clinical_data)
  #################################################################################
  #Variables
  ################
  genes = "GLA"
  cancer = "UVM"
  p_val = 0.00000147611925505053
  #################################
  #################################
  #################################
  rm(list=ls())
  dev.off()
  

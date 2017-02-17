TCGA_Samples = c("BRCA","GBM","OV","LUAD","UCEC","KIRC","HNSC","LGG","THCA","LUSC","PRAD","STAD","SKCM","COAD","BLCA","CESC","KIRP","SARC","LAML","ESCA","PAAD","PCPG","READ","TGCT","THYM","ACC","UVM","DLBC","UCS","CHOL") #All TCGA Samples without MESO

TCGA_Samples = TCGA_Samples[order(TCGA_Samples)]

Mutation_survival = function(ACC, clinical_data)
{
  unique_genes = as.character(unique(ACC$Hugo_Symbol))
  unique_genes = unique_genes[order(unique_genes)]
  genes_acc = list()
  
  for(i in unique_genes)
  {
    genes_acc[[genes]] = unique(as.character(ACC[rownames(ACC[which(ACC$Hugo_Symbol==genes),]),3]))
  }
  
  clinical_match = as.data.frame(clinical_data)
  
  event = c() #To have all the events coded in 0 or 1
  
  days_to_death = c() #days to death, if not, then follow-up data
  
  ##Preparing clinical survival data###
  for ( i in 1:length(as.character(clinical_match$days_to_death))) #Take follow-day where death days are not available
  {
    days_to_death[i] = ifelse(is.na(as.character(clinical_match$days_to_death)[i]),as.numeric(as.character(clinical_match$days_to_last_followup))[i],as.numeric(as.character(clinical_match$days_to_death))[i])
    
    event[i] = as.numeric(as.character(clinical_match$vital_status))[i]
  }
  
  clinical_survival = as.data.frame(cbind(days_to_death,event))
  
  rownames(clinical_survival)= rownames(clinical_match)
  
  ####################
  survival_mutation_matrix = list()
  
  #status = c()
  
  for(j in names(genes_acc))
  {
    status = c()
    
    for(i in rownames(clinical_survival))
    {
      status = c(status, ifelse( i %in% genes_acc[[j]][], 1, 0)) 
    }
    
    survival_mutation_matrix[[j]] = cbind(status)
    rownames(survival_mutation_matrix[[j]]) = rownames(clinical_survival)
  }
  
  results_coxph = list() #Creating a result list
  
  if ( sum(survival_mutation_matrix[[genes]][]) > 0 )
      
    {  
      plot(survfit(Surv(clinical_survival$days_to_death,clinical_survival$event) ~ as.factor(survival_mutation_matrix[[genes]][])),col= c("blue","red"), main=paste(cancer,"Cancer Specific Survival Plot for mutation data"))
          
      mtext(paste("(p=",signif(p_val, 2),")"),3)
    
      legend(x="bottomleft",legend=paste(genes,c("Mut yes","Mut No")),col=c("red","blue"),lwd=1)
    
    
  }
  return(NULL)
  
}
  
  ACC = read.table(list.files("~/Desktop/labRotation1_AnubhavK/Mutation_data/Mutation_data_Formatted",pattern=cancer,full.names = T),header=T,sep="\t")
  
  clinical_data = t(read.table(list.files("~/Desktop/labRotation1_AnubhavK/Clinical_Firehose",pattern=cancer,full.names = T),header=T,sep="\t",row.names = 1))[,-1]
  
  Mutation_survival(ACC,clinical_data)
  
  #################################################################################
  #Variables
  ################
  genes = "DISP2"
  cancer = "OV"
  p_val = 0.000000529721621678192
  #################################
  #################################
  #################################
  rm(list=ls())
  dev.off()
  
  
  

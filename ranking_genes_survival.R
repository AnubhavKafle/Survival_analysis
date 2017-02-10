library(VennDiagram)

TCGA_Samples = c("BRCA","GBM","OV","LUAD","UCEC","KIRC","HNSC","LGG","THCA","LUSC","PRAD","STAD","SKCM","COAD","BLCA","CESC","KIRP","SARC","LAML","ESCA","PCPG","READ","TGCT","THYM","ACC","UVM","DLBC","UCS","CHOL","KICH","LIHC") #All TCGA Samples without MESO and PAAD

TCGA_Samples = TCGA_Samples[order(TCGA_Samples)]

venny = function(Coxph, mutation, Differential)
{
  
  coxph_cancer = as.character(Coxph_significant$X)
  
  mutation_cancer = rownames(mutation_significant)

  overlap1 = calculate.overlap(x = list("coxph_cancer" = coxph_cancer,"mutation_cancer" = mutation_cancer))
  
  common_Expre_Mut_genes = overlap1[[3]]
  
  expression_survival = c()
  mutation_survival = c()
    
  for ( i in common_Expre_Mut_genes)
    
  {
    expression_survival = c(expression_survival, Coxph_significant[Coxph_significant$X==i,2])
    mutation_survival = c(mutation_survival, mutation_significant[i,])
    
  }
  
  common_expression_survival = cbind.data.frame(common_Expre_Mut_genes,expression_survival)
  common_mutation_survival = cbind.data.frame(common_Expre_Mut_genes, mutation_survival)
  common_mutation_survival = common_mutation_survival[order(common_mutation_survival$mutation_survival),]
  common_expression_survival = common_expression_survival[order(common_expression_survival$expression_survival),]

  write.csv(common_expression_survival, file = paste(colnames(Coxph_significant)[2],"ranked_significant_genes_expressionSurvival.csv",sep="")
  write.csv(common_mutation_survival, file = paste(colnames(Coxph_significant)[2],"ranked_significant_genes_MutationSurvival.csv",sep="")
  return(NULL)
  }

#output_results = list()

for ( i in TCGA_Samples)
{
  Differential_significant = NULL
  
  Coxph_significant = read.csv(list.files("/media/pathway/48F918113AE39451/New_results/Significant_coxph",pattern=paste("Formatted_",i,sep=""),full.names = T))
  
  mutation_significant = read.table(list.files("/media/pathway/48F918113AE39451/New_results/Significant_genes_mutation_survival",pattern=paste("Significant_",i,sep=""),full.names = T),sep="\t", row.names = 1)
  
  #if (length(list.files("C:/Users/Anubhav/Desktop/Differential_analysis/New_results/Significant_genes_Differential_analysis",pattern=paste("Significant_",i,sep=""),full.names = T)))Differential_significant = read.table(list.files("C:/Users/Anubhav/Desktop/Differential_analysis/New_results/Significant_genes_Differential_analysis",pattern=paste("Significant_",i,sep=""),full.names = T),header=T,sep="\t",row.names =1)
  
}



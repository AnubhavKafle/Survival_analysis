library(VennDiagram)

TCGA_Samples = c("BRCA","GBM","OV","LUAD","UCEC","KIRC","HNSC","LGG","THCA","LUSC","PRAD","STAD","SKCM","COAD","BLCA","CESC","KIRP","SARC","LAML","ESCA","PAAD","PCPG","READ","TGCT","THYM","ACC","UVM","DLBC","UCS","CHOL","KICH","LIHC") #All TCGA Samples without MESO

TCGA_Samples = TCGA_Samples[order(TCGA_Samples)]

venny = function(Coxph, mutation, Differential)
{
  common_in_all_genes = c()
  
  
  
  coxph_cancer = as.character(Coxph$X)
  
  mutation_cancer = rownames(mutation)
  
  Differential_cancer = rownames(Differential)
  

  overlap1 = calculate.overlap(x = list("coxph_cancer" = coxph_cancer,"mutation_cancer" = mutation_cancer))
    
  common_Expre_Mut_genes = overlap1[[3]]
    
  if(length(Differential_cancer > 0)){
  overlap2 = calculate.overlap(x = list("Differential_cancer" = Differential_cancer, "Expression_mutation_overlap" = common_Expre_Mut_genes))
    
  common_in_all_genes = overlap2[[3]]
  
  }
  
  if(length(common_in_all_genes) > 0) return (list(length(coxph_cancer),length(mutation_cancer), length(common_Expre_Mut_genes),length(common_in_all_genes)))
  else return(list(length(coxph_cancer),length(mutation_cancer), length(common_Expre_Mut_genes),NA))
  
}

output_results = list()
  
Expression_survival = vector()
Mutation_survival = vector()
Common_Expression_Mutation_survival = vector()
Common_in_all = vector()
  for ( i in TCGA_Samples)
  {
    Differential_significant = NULL
    
    Coxph_significant = read.csv(list.files("C:/Users/Anubhav/Desktop/Differential_analysis/New_results/Significant_coxph",pattern=paste("Formatted_",i,sep=""),full.names = T))
    
    if(length(list.files("C:/Users/Anubhav/Desktop/Differential_analysis/New_results/Significant_genes_mutation_survival",pattern=paste("Significant_",i,sep=""),full.names = T)))mutation_significant = read.table(list.files("C:/Users/Anubhav/Desktop/Differential_analysis/New_results/Significant_genes_mutation_survival",pattern=paste("Significant_",i,sep=""),full.names = T),sep="\t", row.names = 1)
    
    if (length(list.files("C:/Users/Anubhav/Desktop/Differential_analysis/New_results/Significant_genes_Differential_analysis",pattern=paste("Significant_",i,sep=""),full.names = T)))Differential_significant = read.table(list.files("C:/Users/Anubhav/Desktop/Differential_analysis/New_results/Significant_genes_Differential_analysis",pattern=paste("Significant_",i,sep=""),full.names = T),header=T,sep="\t",row.names =1)
    
    results = as.matrix(venny(Coxph_significant, mutation_significant, Differential_significant))
    
    Expression_survival = c(Expression_survival, results[[1]])
    Mutation_survival = c(Mutation_survival, results[[2]])
    Common_Expression_Mutation_survival = c(Common_Expression_Mutation_survival, results[[3]])
    Common_in_all = c(Common_in_all, results[[4]])
  
  }
table = cbind(Expression_survival, Mutation_survival, Common_Expression_Mutation_survival, Common_in_all)
rownames(table) = TCGA_Samples
write.csv(table, file = "results_summary.csv")
  

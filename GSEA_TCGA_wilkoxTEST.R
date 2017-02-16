TCGA_samples = c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

GSEA_test = function(expression_survival_common, mutation_survial_common)
{
  featured_genes = as.character(expression_survival_common$common_Expre_Mut_genes) #genes with significant value
  
  pw_genesNrs = list()
  pw_genes = list()
  for ( pw in pw_names ) {
    pathway_genes = pwgenes[[pw]] #genes in a particular pathway
    matched_pathway_genes = pathway_genes %in% featured_genes # Matching genes in a particular pathway  
    pw_genes[[pw]] = unique(pathway_genes[matched_pathway_genes])#Filtering genes that dont match  
    pw_genesNrs[[pw]] = length(pathway_genes[matched_pathway_genes]) # maintaining a record of no of genes left
  }
  
  W_pValue = list()
  for ( pw in pwnames ) {
    reduced_genes = pw_genes[[pw]]
    pw_ranks = apply(cbind(match(reduced_genes, expression_survival_common$common_Expre_Mut_genes), match(reduced_genes, mutation_survial_common$common_Expre_Mut_genes)),1,sum)
    if(!(length(pwranks) == 0)) {
      wilcox_test = wilcox.test(pw_ranks,(1:length(expression_survival_common$common_Expre_Mut_genes))[-pwranks],
                                alternative="less")
      W_pValue[count] = wilcox_test$p.value
    } else { W_pValue[count] = 1 }
    ###### WILCOX ENDE
  }
  
  
}

# split up info into pwnames, description and gene lists of same length
load("fwdpwanalyses/pwgenes.RData")
pathway_details=read.csv(file="fwdpwanalyses/pathways.csv")
pw_names = as.character(pathway_details$pw)
pw_description = as.character(pathway_details$desc)
rm(pathway_details)

for ( i in TCGA_Samples)
{
  
  Expression_survival_pVal = read.csv(list.files("~/Desktop/labRotation1_AnubhavK/New_results_pendrive_backup/significant_survival_ranked_Common_genes/",pattern=i,full.names = T)[1])
  
  Mutation_survival_pVal = read.csv(list.files("~/Desktop/labRotation1_AnubhavK/New_results_pendrive_backup/significant_survival_ranked_Common_genes/",pattern=i,full.names = T)[2])

  wilcox_test_GSEA = GSEA_test(Expression_survival_pVal, Mutation_survival_pVal)
    
}



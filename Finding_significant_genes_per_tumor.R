cox_pvalues = read.table("TCGA_Cox_Pval_Allgenes.txt", sep= "\t", header = T, row.names = 1)

cox_hazardRatio = read.table("TCGA_Cox_Hazardratio.txt", sep="\t", header = T , row.names = 1)

cox_pvalues[!is.finite(as.matrix(cox_pvalues))] <- 0

cox_hazardRatio[!is.finite(as.matrix(cox_hazardRatio))] <- 0

genes_significantP_greaterHazard = list()

for ( i in 1:length(colnames(cox_hazardRatio)))
{
  for(j in 1:length(rownames(cox_hazardRatio)))
  {
    if ( cox_pvalues[j,i] <= 0.05 & cox_hazardRatio[j,i] >= 1)
    {
      genes_significantP_greaterHazard[[rownames(cox_pvalues)[j]]] = c(cox_pvalues[j,i], cox_hazardRatio[j,i])
    }
  }
  significant_genes = as.matrix(genes_significantP_greaterHazard)
  
  write.table(significant_genes, file = paste(colnames(cox_hazardRatio)[i],".txt"), sep = "\t")
}

Expression_survival = read.table("Genes_with_only_Pvalues_new_coxph.txt", header = T, sep = "\t", row.names = 1)

#head(Expression_survival[1:5,1:5])

#head(colnames(Expression_survival))

for(j in 1:length(colnames(Expression_survival)))
{
  cancer = Expression_survival[which(Expression_survival[,j] <= 0.05),]
  
  write.table(cancer[,j,drop = FALSE], file = paste(colnames(Expression_survival)[j],"_Significant.txt"), sep = "\t", row.names = T)
              
}

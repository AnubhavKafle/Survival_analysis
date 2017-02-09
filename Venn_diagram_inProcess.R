library(VennDiagram)

grid.newpage()
require(gridExtra)
#grid.arrange(gTree(children=ss), top="BRCA")

TCGA_Samples = c("BRCA","GBM","OV","LUAD","UCEC","KIRC","HNSC","LGG","THCA","LUSC","PRAD","STAD","SKCM","COAD","BLCA","CESC","KIRP","SARC","LAML","ESCA","PAAD","PCPG","READ","TGCT","THYM","ACC","UVM","DLBC","UCS","CHOL","KICH","LIHC") #All TCGA Samples without MESO

TCGA_Samples = TCGA_Samples[order(TCGA_Samples)]

i = TCGA_Samples[2]

length(Differential_significant)
venny = function(Coxph, mutation, Differential)
{
  
  coxph_cancer = as.character(Coxph_significant$X)
  
  mutation_cancer = rownames(mutation_significant)
  
  Differential_cancer = rownames(Differential_significant)
  
  if(length(Differential_cancer != 0))
  {
    
    overlap = calculate.overlap(x = list("coxph_cancer" = coxph_cancer,"mutation_cancer" = mutation_cancer,"Differential_cancer" = Differential_cancer))
    
   venn_plot =  draw.triple.venn(area1 = length(coxph_cancer), area2 = length(mutation_cancer), area3 = length(Differential_cancer), n123 = length(overlap[[1]]), n12 = length(overlap[[2]])+length(overlap[[1]]), n13 = length(overlap[[3]])+length(overlap[[1]]),
      n23 = length(overlap[[4]])+length(overlap[[1]]), category = c("Expression_Survival", "Mutation_survival", "Differential_genes"), lty = "blank", fill = c("skyblue", "pink1", "mediumorchid"),scaled = F)  
   grid.arrange(gTree(children=venn_plot), top=colnames(Coxph_significant)[2])
  }
  
  else
  {
    
    overlap1 = calculate.overlap(x = list("coxph_cancer" = coxph_cancer,"mutation_cancer" = mutation_cancer))
    
    #draw.pairwise.venn(area1 = length(coxph_cancer), area2 = length(mutation_cancer),cross.area = length(overlap1[[3]]),category = c("Expression_Survival", "Mutation_survival"),lty = "blank", fill = c("skyblue", "pink1")) 
    venny_plot = draw.pairwise.venn(area1 = length(coxph_cancer), area2 = length(mutation_cancer),cross.area = length(overlap1[[3]]),category = c("Expression_Survival", "Mutation_survival"),lty = "blank", fill = c("skyblue", "pink1"),cat.pos = c(0,0), cat.dist = rep(0.025,2),scaled = F)
    
    grid.arrange(gTree(children = venny_plot),top = colnames(Coxph_significant)[2])

  }
  
  
}

for ( i in TCGA_Samples)
{
  Differential_significant = NULL
  
  Coxph_significant = read.csv(list.files("/media/pathway/48F918113AE39451/New_results/Significant_coxph",pattern=paste("Formatted_",i,sep=""),full.names = T))
  
  mutation_significant = read.table(list.files("/media/pathway/48F918113AE39451/New_results/Significant_genes_mutation_survival",pattern=paste("Significant_",i,sep=""),full.names = T),sep="\t", row.names = 1)
  
  if (length(list.files("/media/pathway/48F918113AE39451/New_results/Significant_genes_Differential_analysis",pattern=paste("Significant_",i,sep=""),full.names = T)))Differential_significant = read.table(list.files("/media/pathway/48F918113AE39451/New_results/Significant_genes_Differential_analysis",pattern=paste("Significant_",i,sep=""),full.names = T),header=T,sep="\t",row.names =1)
  
  venn_plot = venny(Coxph_significant, mutation_significant, Differential_significant)
  
}

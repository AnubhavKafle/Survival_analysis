files = list.files(".", pattern = "*.txt")
i = files[1]
i
par(mfrow = c(3,3))

MA_plot = function()
{
  for ( i in files)
{
  Cancer = read.table(i, sep = "\t", header = T)
  
  color = ifelse((Cancer$log2FoldChange >=1 | Cancer$log2FoldChange <= -1) & Cancer$padj <= 5, "red", "black")
  
  plot(log10(Cancer$baseMean), Cancer$log2FoldChange, pch =20, col = color, main = unlist(strsplit(i,"\\."))[[1]] )
  
  abline(0,0, col = "blue", lwd = 2)

  abline(1,0,col = "blue", lwd = 2)
  
  abline(-1,0,col = "blue", lwd = 2)
  
  pdf(file = paste(unlist(strsplit(i,"\\."))[[1]],"_MA_plots.pdf"))
  } 
} 

MA_plot()

dev.off()

rm(list=ls())

volcano_plot = function()
{
  for ( i in files)
{
  Cancer = read.table(i, sep = "\t", header = T)
  
  color = ifelse((Cancer$log2FoldChange >=1 | Cancer$log2FoldChange <= -1) & Cancer$padj <= 5, "red", "black")
  
  plot(Cancer$log2FoldChange, Cancer$padj, pch =20, col = color, main = unlist(strsplit(i,"\\."))[[1]], ylim = rev(range(c(min(Cancer$padj[!is.na(Cancer$padj)]),max(Cancer$padj[!is.na(Cancer$padj)])))))
  
  abline(v=0, col = "blue", lwd = 2)
  
  #pdf(file = paste(unlist(strsplit(i,"\\."))[[1]],"_Volcano_plots.pdf"))
} 
  
}  
  
dev.off()

rm(list=ls())



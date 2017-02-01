library(survival)

TCGA_Samples = c("BRCA","GBM","OV","LUAD","UCEC","KIRC","HNSC","LGG","THCA","LUSC","PRAD","STAD","SKCM","COAD","BLCA","CESC","KIRP","SARC","LAML","ESCA","PAAD","PCPG","READ","TGCT","THYM","ACC","MESO","UVM","DLBC","UCS","CHOL") #All TCGA Samples

TCGA_Samples = TCGA_Samples[order(TCGA_Samples)]

survival = function(expression1,clinical)
{
#expression1=t(read.table("~/Desktop/scripts/ACC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",header=T,sep="\t",row.names=1))[,-1]

#clinical=t(read.table("~/Desktop/scripts/ACC-Clinical.txt",header=T,sep="\t",row.names = 1))[,-1]

event = c() #To have all the events coded in 0 or 1

days_to_death = c() #days to death, if not, then follow-up data

expression_rna = t(apply(expression1,1,as.numeric)) #TO eliminate factor problems

total_gene_list = colnames(expression1) #Obtaining all the genes in the sample

rownames(expression_rna) <- as.vector(sapply(rownames(expression1), function(x) tolower(paste(unlist(strsplit(x,"\\."))[[1]],unlist(strsplit(x,"\\."))[[2]],unlist(strsplit(x,"\\."))[[3]],sep="."))))#To keep all the rows in small letters

colnames(expression_rna) = colnames(expression1)  ##############

rm(expression1)  ################

expression_rna_countMorethanZero = expression_rna[,which(apply(expression_rna, 2, median) >= 0.1)]  ####to filter out gene that have counts less than 0.1 in 50 or more of the sample.

expression_rna_log2 = log2(t(apply(t(apply(expression_rna_countMorethanZero,1,as.numeric)), 1, function(i) i+1 ))) ### Log2 transform of the filtered gene counts

colnames(expression_rna_log2) = colnames(expression_rna_countMorethanZero) #################

match_index = which(rownames(clinical)%in% rownames(expression_rna_log2)) #Check which patient ID are matching with the tumor data

clinical_match=as.data.frame(clinical[match_index,]) #filter the clinical match. There are Normal patients as well so less mathches


#genes_survival_pval=data.frame() #storing p values in a data frame

for ( i in 1:length(as.character(clinical_match$days_to_death))) #Take follow-day where death days are not available
{
  days_to_death[i] = ifelse(is.na(as.character(clinical_match$days_to_death)[i]),as.numeric(as.character(clinical_match$days_to_last_followup))[i],as.numeric(as.character(clinical_match$days_to_death))[i])
  
  event[i] = as.numeric(as.character(clinical_match$vital_status))[i]
}

rm(clinical)

clinical_survival = as.data.frame(cbind(days_to_death,event)) #Create a data frame for survival fit analysis

rownames(clinical_survival)= rownames(clinical_match)#Matching rownames

expression_rna_ordered_log2_median = expression_rna_log2[rownames(clinical_survival),] #Creating expression matrix ordered according to rownames of clinical_match
######## NOT NEEDED #######clinical_Survival=as.data.frame(days_todeath=ifelse(is.na(clinical$days_to_death),as.numeric(as.character(clinical$days_to_last_followup)),as.numeric(as.character(clinical$days_to_death))),event=as.numeric(as.character(clinical$vital_status)),row.names=rownames(clinical))

rm(expression_rna,expression_rna_countMorethanZero,match_index,expression_rna_log2,event, days_to_death,clinical_match)

results_coxph=list() #Creating a result list

for (genes in colnames(expression_rna_ordered_log2_median))
{
  model_analysis = summary(coxph(Surv(clinical_survival$days_to_death,clinical_survival$event) ~ as.numeric(expression_rna_ordered_log2_median[,genes]))) #P-value is 5th in the index of coefficient. Coxph model is used because of the continous event of gene expression
  
  results_coxph[[genes]] = c(model_analysis$coefficients[5], model_analysis$coefficients[2])
}

p_values = vector() #An empty vector to extract P_values from the results_coxph

hazard_ratio = vector() # To keep hazard ratio data 

p_values_and_HazardRatio_geneName = list()

hazard_ratio_geneName = list()
  
for (i in 1:length(results_coxph))
{
  p_values = c(p_values,results_coxph[[i]][1])
  
  hazard_ratio = c(hazard_ratio, results_coxph[[i]][2])
}
  
p_values.adjusted = p.adjust(p_values, method = "BH")
  
gene_names = names(results_coxph) #To extract all the gene names 

for( ii in 1:length(gene_names))
{
  p_values_and_HazardRatio_geneName[[gene_names[ii]]] = c(p_values.adjusted[ii], hazard_ratio[ii])
  
  #hazard_ratio_geneName[[gene_names[ii]]] = hazard_ratio[ii]
  
}
 
p_values_matched_genes = list() 

for (number in 1:length(total_gene_list))
{
  P_values_matched_genes[[total_gene_list[number]]] = ifelse(total_gene_list[number] %in% names(p_values_and_HazardRatio_geneName), p_values_and_HazardRatio_geneName[[which(names(p_values_and_HazardRatio_geneName) %in% total_gene_list[number])]][1], NA)
 
  hazard_ratio_geneName[[total_gene_list[number]]] = ifelse(total_gene_list[number] %in% names(p_values_and_HazardRatio_geneName), p_values_and_HazardRatio_geneName[[which(names(p_values_and_HazardRatio_geneName) %in% total_gene_list[number])]][2], NA) 
}

Sys.time()

return(list(as.numeric(as.matrix(P_values_matched_genes)), names(P_values_matched_genes), as.numeric(as.matrix(hazard_ration_geneName)))) #Return as a list to the function call assignment variable

}

##############
#Survival analysis for all the tumor samples
##############

survival_results=list()  #Function will be called using this variable

p_value_gene_matrix = c()  #This is to filter p_values from the return list below for the function

hazard_ratio_matrix = c()
                                   
for ( i in TCGA_Samples)
{
  
  Cancer_Expression = t(read.table(list.files("~/Desktop/labRotation1_AnubhavK/gene_expression/",pattern=i,full.names = T),header=T,sep="\t",row.names=1))[,-1]
  
  clinical_data = t(read.table(list.files("~/Desktop/labRotation1_AnubhavK/Clinical_Firehose/",pattern=i,full.names = T),header=T,sep="\t",row.names = 1))[,-1]
  
  survival_results = survival(Cancer_Expression,clinical_data)
  
  p_value_gene_matrix[[i]] = as.numeric(unlist(survival_results[1]))  #For each gene, P values will be coerced as list elements.
  
  hazard_ratio_matrix[[i]] = as.numeric(unlist(survival_results[3]))
}  

Gene_with_survival_adj.p_values = as.data.frame(p_value_gene_matrix) #Create a data frame from list 

rownames(Gene_with_survival_adj.p_values) = unlist(survival_results[2]) #Assign gene names to the dataframe

Genes_with_hazard_ratio =  as.data.frame(hazard_ratio_matrix)

rownames(Genes_with_hazard_ratio) = unlist(survival_results[2])

write.table(Genes_with_hazard_ratio, file = "Genes_with_hazard_ratio.txt",sep = "\t", row.names = T)

write.table(Gene_with_survival_adj.p_values, file = "Genes_with_adjPvalues.txt", sep = "\t", row.names = T)



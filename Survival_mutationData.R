###### This is to format the file suitable for Surivival analysis.Done in BASH SHELL ########
for i in $(ls *.txt) ; do awk -F'\t' '{print ($1,$9,tolower($16)) }' $i | grep -v -E "Silent|RNA" | sed 's/ /\t/g'  >> "Filtered_"$i; done 
for i in $(ls Filtered*) ; do  cut -d"-" -f1,2,3 $i >> "filtered_"$i; done
rm Filtered*
mv filtered* ./Mu
rename 's/filtered_Filtered_/Filtered_/' filtered*
for i in $(ls *.txt) ; do  sed 's/Tumor_Sample_Barcode/Patience_ID/g' $i >> "tested_"$i; done
rm Filtered*
cd Mu/
rename 's/tested_//g' tested_Filtered_*

########################################
#### R code for survival analysis  for the mutation data
ACC = read.table("ACC-Mutations-AllSamples.txt", head =T, sep="\t")
unique_genes = as.character(unique(ACC$Hugo_Symbol))
genes_acc = list()
#unique(as.character(ACC[rownames(ACC[which(ACC$Hugo_Symbol=="NOL9"),]),3]))

for(i in unique_genes)
{
  genes_acc[[i]] = unique(as.character(ACC[rownames(ACC[which(ACC$Hugo_Symbol==i),]),3]))
}

mutation_matrix = list()
for ( i in names(genes_acc))
{
  
  mutation_matrix[[i]] = cbind(rep(1,length(genes_acc[[i]])))
  rownames(mutation_matrix[[i]]) = as.character(genes_acc[[i]])
   colnames(mutation_matrix[[i]])= "Status"
}

  clinical_data = t(read.table(list.files("~/Desktop/labRotation1_AnubhavK/Clinical_Firehose/",pattern="ACC",full.names = T),header=T,sep="\t",row.names = 1))[,-1]

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

status = c()

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
  




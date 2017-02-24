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

rm(expression_rna,expression_rna_countMorethanZero,match_index,expression_rna_log2,event, days_to_death,clinical_match)
Survival_times = sapply(as.numeric(clinical_survival$days_to_death), function(x)x/(30*12))
clinical_survival = cbind(Survival_times, clinical_survival)
######## NOT NEEDED #######clinical_Survival=as.data.frame(days_todeath=ifelse(is.na(clinical$days_to_death),as.numeric(as.character(clinical$days_to_last_followup)),as.numeric(as.character(clinical$days_to_death))),event=as.numeric(as.character(clinical$vital_status)),row.names=rownames(clinical))
expression1 = t(read.table(list.files("~/Desktop/labRotation1_AnubhavK/gene_expression/",pattern=i,full.names = T),header=T,sep="\t",row.names=1))[,-1]

clinical = t(read.table(list.files("~/Desktop/labRotation1_AnubhavK/Clinical_Firehose/",pattern=i,full.names = T),header=T,sep="\t",row.names = 1))[,-1]
###########################################################################
##########################################################################
significant_genes =  read.csv("/home/pathway/Desktop/labRotation1_AnubhavK/New_results_pendrive_backup/significant_survival_ranked_Common_genes/ACC_ranked_significant_genes_expressionSurvival.csv")
Significant_expression_matrix = matrix()
index =c()

genes = sapply(as.vector(colnames(expression_rna_ordered_log2_median)),function(x)unlist(strsplit(x,"\\|"))[1])

for (i in as.vector(significant_genes$common_Expre_Mut_genes)){
  
  value = grep(paste("^",i,"\\b",sep=""), as.vector(genes))
  if(length(value) == 1)index=c(value, index)#index = c(grep(i,as.vector(genes)),index)
  #Significant_expression_matrix = cbind(as.vector(expression_rna_ordered_log2_median[,grep(paste("^",i,sep = ""),colnames(expression_rna_ordered_log2_median))]),Significant_expression_matrix)
  
  
}
index=index[order(index)]
Significant_expression_matrix = as.data.frame(expression_rna_ordered_log2_median[,index])
#####fitting parametric survival regression model on the data

#regression_fit = survreg(Surv(clinical_survival$death_years,clinical_survival$event)~.,data = Significant_expression_matrix,dist = "weibull")

#regression_fit = survreg(Surv(clinical_survival$days_to_death,clinical_survival$event)~Significant_expression_matrix$`A2ML1|144568`,Significant_expression_matrix,dist = "weibull")

#model_analysis = coxph(Surv(clinical_survival$death_years,clinical_survival$event) ~ ., data = expression_rna_ordered_log2_median) #P-value is 5th in the index of coefficient. Coxph model is used because of the continous event of gene expression
#################3
#PAMR survival modelling 
test_data_expression  = Significant_expression_matrix[-(1:ceiling(length(rownames(Significant_expression_matrix))/2)),]

training_data_expression  = Significant_expression_matrix[1:ceiling(length(rownames(Significant_expression_matrix))/2),]

clinical_survival_training = clinical_survival[(1:ceiling(length(rownames(clinical_survival))/2)),]

clinical_survival_test = clinical_survival[-(1:ceiling(length(rownames(clinical_survival))/2)),]

cancer_survival_training_pamr = list(x = as.matrix(t(training_data_expression)), survival.time = clinical_survival_training$Survival_times,censoring.status = clinical_survival_training$event, genenames  = colnames(training_data_expression) )

cancer_survival_test_pamr = list(x = as.matrix(t(test_data_expression)), survival.time = clinical_survival_test$Survival_times,censoring.status = clinical_survival_test$event, genenames  = colnames(test_data_expression) )

survival_train_model = pamr.train(cancer_survival_training_pamr, ngroup.survival = 4)

survival_class = pamr.surv.to.class2(cancer_survival_test_pamr$survival.time, cancer_survival_test_pamr$censoring.status, n.class = survival_train_model$ngroup.survival)$prob

#table(apply(survival_class, 1, function(x)which(x==max(x))))

yhat = pamr.predict(survival_train_model, cancer_survival_pamr$x, threshold = 1.0)
head(yhat)
pamr.plotsurvival(yhat, cancer_survival_pamr$survival.time, cancer_survival_pamr$censoring.status)
pamr.confusion.survival(survival_train_model, survival.time = cancer_survival_pamr$survival.time, cancer_survival_pamr$censoring.status,yhat)
##################

BLCA_RNASeq=read.table("BLCA.mRNAseq_raw_counts.txt",header = T,sep="\t")

BLCA_Cancer=BLCA_RNASeq[,grep("*\\.01$", colnames(BLCA_RNASeq))]

BLCA_Cancer=cbind(BLCA_RNASeq[,1],BLCA_Cancer)
colnames(BLCA_Cancer)[1]= "Gene_symbol"


BLCA_Normal=BLCA_RNASeq[,grep("*\\.11$", colnames(BLCA_RNASeq))]

BLCA_Normal=cbind(BLCA_RNASeq[,1],BLCA_Normal)
colnames(BLCA_Normal)[1]= "Gene_symbol"
match=c()

ss=c()
for ( i in colnames(BLCA_Cancer))
{
  a=paste(unlist(strsplit(i,"[.]"))[c(1,2,3)],sep=".",collapse = ".") 
  for ( j in colnames(BLCA_Normal))
  {
   ss=c(ss,paste(unlist(strsplit(j,"[.]"))[c(1,2,3)],sep=".",collapse = "."))
   if (length(intersect(a,ss)== 1))
  {
    match = c(match,i)
    break
   }
    ss=c()
  }
}
BLCA_Cancer_Match=BLCA_Cancer[,match]

write.table(BLCA_Cancer_Match,file="BLCA_Cancer_Matched.txt",sep="\t",row.names=F)

write.table(BLCA_Normal,file="BLCA_Normal_Matched.txt",sep="\t",row.names=F)


#############################################################

















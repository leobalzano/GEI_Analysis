##############################################
####     AMMImodelsprocessingandPlot.R     ###
##############################################
# This script is to perform and plot AMMI models including the pre-processing of the data

##############################################
# Libraries:
library(agricolae)

# FUNCTIONS:
AMMIprocessing<-function (X){
  LaTabla<-X
  LaTabla<-data.frame(Genotype=rownames(LaTabla),LaTabla)
  vectito<-c(rep(NA, dim(LaTabla)[1]))
  VectorTF<-c(rep(NA, dim(LaTabla)[1]))
  VectorOTE<-c(rep(0, dim(LaTabla)[1]))
  namevar<-LaTabla[1,1]
  #namevar<-rownames(LaTabla)[1]
  
  if(sum(namevar==LaTabla[,1])>1){
  #if(sum(namevar==rownames(LaTabla))>1){
    for (a in 1: length(unique (LaTabla[,1]))) {
      #print (a)
      namevar<-LaTabla[a,1]
      VectorTF<-namevar==LaTabla[,1]
      VectorTF[VectorTF==TRUE]<-seq(1:length(VectorTF[VectorTF==TRUE]))  
      #print(VectorTF)
      VectorOTE<-(VectorOTE +VectorTF )
    }
    
    #print(VectorOTE)
    LaTabla$Replicates<-VectorOTE
  } else {#print("Warning: Your data does not have replicates")
    LaTabla<-rbind(LaTabla,LaTabla)
    vectito<-c(rep(NA, dim(LaTabla)[1]))
    VectorTF<-c(rep(NA, dim(LaTabla)[1]))
    VectorOTE<-c(rep(0, dim(LaTabla)[1]))
    namevar<-LaTabla[1,1]
    #replicas<-length(unique(LaTabla[,1]))
    
    if(sum(namevar==LaTabla[,1])>1){
      for (a in 1: length(unique (LaTabla[,1]))) {
        #print (a)
        namevar<-LaTabla[a,1]
        VectorTF<-namevar==LaTabla[,1]
        VectorTF[VectorTF==TRUE]<-seq(1:length(VectorTF[VectorTF==TRUE]))  
        #print(VectorTF)
        VectorOTE<-(VectorOTE +VectorTF )
      }
      
      #print(VectorOTE)
      LaTabla$Replicates<-VectorOTE
    }
  }
  LaTabla<-LaTabla
  return(LaTabla)
  
}
#####

# TCH <- read.delim("~/Dropbox (UFL)/Projects2018/GEIPlatform/Analisis Ramon Rea Megaambiente/TCH.txt",sep="\t", row.names=1)
# X<-TCH
# X<-rbind(TCH,TCH)
# Tabelina<- AMMIprocessing(TCH)
# 
# # Tabelina<-AMMIprocessing(TCH2X)
#  MeltedTable<-melt(Tabelina,id = c("Genotype","Replicates"))
# # 
# ammi<-AMMI(ENV=MeltedTable$variable, GEN =MeltedTable$Genotype , Y=MeltedTable$value,
#             REP =MeltedTable$Replicates )
# plot(ammi, main= "Modelo AMMI para PolP", type =1)
# # plot(ammi,0,1)
# ammi
# #Para obtener los ASV (AMMI Stability Value)
# Idx<-index.AMMI(ammi)
# Idx

##############################################
########     PCAforWideMatrices.R     ########
##############################################
# This script is to perform and plot PCA's for wide matrices,
# meaning n<<p matrices

##############################################
# FUNCTIONS:
# PCA itself
PCAforWideMatrices<-function(X)
{
  #PCA.GENES is very useful to obtain principal components to a matrix that has more variables than individuals. 
  #R can not apply princomp is such case and when there are a lot of variables eigen(t(X)%*%X) can not be computed.
  
  #X is a matrix that has on columns the genes considered as variables in the PCA analysis.
  #First we center the matrix by columns (Xoff) and then we obtain the eigenvalues and the eigenvectors of the matrix Xoff%*%t(Xoff) and we #use the equivalences between the loadings and scores to obtain the solution
  #Llamo scores1 y loadings1 a lo que busco y scores2 y loadings2 a los scores y loadings de la transpuesta
  X <- as.matrix(X)
  n<-ncol(X)
  p<-nrow(X)
  offset<-apply(X,2,mean)
  Xoff<- X-(cbind(matrix(1,p,1))%*%rbind(offset))
  
  #eigen command sorts the eigenvalues in decreasing orden.
  
  eigen<-eigen(Xoff%*%t(Xoff)/(p-1))
  var<-cbind(eigen$values/sum(eigen$values),cumsum(eigen$values/sum(eigen$values)))
  
  loadings2<-eigen$vectors
  scores2<-t(Xoff)%*%loadings2
  
  normas2<-sqrt(apply(scores2^2,2,sum))
  
  scores1<-loadings2%*%diag(normas2)
  loadings1<-scores2%*%diag(1/normas2)
  
  output<-list(eigen,var,scores1,loadings1, Xoff)
  names(output)<-c("eigen","var.exp","scores","loadings", "Xoff")
  output
}
##############################################
# Plot of the PCA
PlotPCAWide <-function (PCAofX,A1=1,A2=2, Groups, Title="PCA plot", LegendPos="bottomleft", xyLims=NULL,ABline=TRUE) {
# Function to plot a PCA from PCAforWideMatrices
# INPUT: PCAofX: The result of the PCA from PCAforWideMatrices
#            A1: Axis X of the plot
#            A2: Axis Y of the plot
#        Groups: A vector indicating tho which group, a particular row belongs
#         Title: The plot title
#     LegendPos: The position of the legend of the plot
#        xyLims: A vector indicating the limits of x and y axes in the following order, x1,x2,y1,y2
#        ABline: Whether or not to draw the zero,zero axis
  
  
t<-PCAofX$scores
u<-PCAofX$loadings

V <- PCAofX$var.exp ; dim(V)  # varianza explicada

Names <- rownames(PCAofX$Xoff)
Groups<-as.factor(Groups)
ColorGroups <- 1:nlevels(Groups)
ColorInd <- as.integer(Groups)


if (is.null(xyLims)){
plot(t[,A1], t[,A2], pch = 21, cex = 2,
     col = "black", bg=ColorInd,
     xlab = paste("PC",A1,": ", round(V[A1,1]*100, digits = 2), "%",sep = ""),
     ylab = paste("PC",A2,": ", round(V[A2,1]*100,digits = 2), "%",sep = ""),
     main= Title
     
)
text(t[,A1], t[,A2],labels=  Names,
     cex= 1, pos=1, col=ColorInd )

legend(LegendPos,legend= unique(Groups) , cex=1,
       col= unique(ColorInd), 
       pch = 16, bty = "o", box.col = "black")
} else {
  plot(t[,A1], t[,A2], pch = 21, cex = 2,
       col = "black", bg=ColorInd,
       xlab = paste("PC",A1,": ", round(V[A1,1]*100, digits = 2), "%",sep = ""),
       ylab = paste("PC",A2,": ", round(V[A2,1]*100,digits = 2), "%",sep = ""),
       main= Title,
       xlim = c(xyLims[1],xyLims[2]), ylim = c(xyLims[3],xyLims[4])
  )
  text(t[,A1], t[,A2],labels=  Names,
       cex= 1, pos=1, col=ColorInd )
  
  legend(LegendPos,legend= unique(Groups) , cex=1,
         col= unique(ColorInd), 
         pch = 16, bty = "o", box.col = "black")
  
  
}

if (ABline) {
abline(h=0,v=0,col="gray55")
} else {}
}

GGPlotPCAWide<-function (PCAofX,A1=1,A2=2, Groups, Title="PCA plot", xyLims=NULL,
                         sizeSpot=2, IndsNames="yes") {
  #library("plotly", lib.loc="/usr/local/lib/R/site-library")
  #library("plotly")
  t<-PCAofX$scores
  u<-PCAofX$loadings
  
  V <- PCAofX$var.exp ; dim(V)  # varianza explicada
  t_v2<-data.frame(t[,c(A1,A2)])
  rownames(t_v2)<-rownames(PCAofX$Xoff)
  
  t_v2$Groups = Groups
  t_v2$Names = rownames(t_v2)
  
  if (IndsNames=="yes"){
    ggplot(as.data.frame(t_v2), aes(X1, X2, color=Groups, label=Names)) + geom_point(size=sizeSpot) +
      xlab(paste("PC",A1,": ", round(V[A1,1]*100, digits = 2), "%",sep = "")) +
      ylab(paste("PC",A2,": ", round(V[A2,1]*100,digits = 2), "%",sep = "")) + ggtitle(Title) + geom_text(size=sizeSpot*2, vjust=2)
    
  } else{
    Groups<-as.factor(Groups)
    ColorGroups <- 1:nlevels(Groups)
    ColorInd <- as.integer(Groups)
    
    
    if (is.null(xyLims)){
      ggplot(as.data.frame(t_v2), aes(X1, X2, color=Groups)) + geom_point(size=sizeSpot) +
        xlab(paste("PC",A1,": ", round(V[A1,1]*100, digits = 2), "%",sep = "")) +
        ylab(paste("PC",A2,": ", round(V[A2,1]*100,digits = 2), "%",sep = "")) + ggtitle(Title)
      
    } else {
      ggplot(as.data.frame(t_v2), aes(X1, X2, color=Groups)) + geom_point(size=sizeSpot) +
        xlab(paste("PC",A1,": ", round(V[A1,1]*100, digits = 2), "%",sep = "")) +
        ylab(paste("PC",A2,": ", round(V[A2,1]*100,digits = 2), "%",sep = "")) + ggtitle(Title) +
        xlim(xyLims[1],xyLims[2]) + ylim (xyLims[3],xyLims[4])
    }
  }
}


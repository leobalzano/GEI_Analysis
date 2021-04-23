#############################################
##########     GEIplatform.R     ############
#############################################
# Author: Leandro Balzano-Nogueira
# University of Florida

# This script is build to create a visualization of different 
# Genotype x Environment analyses under a visualization through Shiny

#############################################
# Libraries:
# Desktop
library(agricolae)
library(DT)
library(shiny)
library(shinyjs)
library(shinydashboard)
library(GGEBiplots)
library(GGEBiplotGUI)
library(reshape2)
library(ade4)
library("plotly")
#library(adegraphics)


# Laptop
# library(DT)
# library("agricolae", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
# library("shiny", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
# library("shinyjs", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
# library("shinydashboard", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
# library("GGEBiplots", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
# library("GGEBiplotGUI", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
# library("reshape2", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
#############################################
## server.R ##
source("PCAforWideMatrices.R")
source("AMMImodelsprocessingandPlot.R")


shinyServer(function(input, output, session) {
    infotableRaw <- reactiveVal(NULL);
    infotable <- reactiveVal(NULL);
  
    observeEvent({
      input$file
      input$checkbox
    }, {
      if (is.null(input$file))
        return(NULL)
                  
      infotableRaw(read.table(input$file$datapath, sep="\t", header=input$checkbox))
      infotable(read.table(input$file$datapath, sep="\t", header=input$checkbox,row.names = 1))
    })
    
    observeEvent({
      infotableRaw()
    }, {
      table1 <- infotable();
      
      if (!all(apply(table1,2,is.numeric))) {
        showModal(modalDialog("Warning: All your data must be numeric and it is not", title="Incorrect table format", size="l"))
        infotable(NULL)
        return(NULL)
      }
      
      ## Dealing with missing values
      for(i in 1:ncol(table1)){
        table1[is.na(table1[,i]), i] <- mean(table1[,i], na.rm = TRUE)
      }
      infotable(table1)
    })
    
    output$value <- DT::renderDT({
      table1 <- infotableRaw();
      table1
    })
    
    ## PCA
    output$PCA <- renderPlot({
      table1 <- infotable();
      if (is.null(table1)) {
        sregData(NULL)
        return(NULL)
      }
      
      pcaResults <- PCAforWideMatrices(table1)
      Grupos<-dim(table1)[1]
      
      GGPlotPCAWide (PCAofX=pcaResults,A1=1,A2=2, 
                     Groups=Grupos, 
                     Title="PCA plot", 
                     xyLims=NULL,IndsNames = "yes")
      
    })
    
    ## SREG model
    sregData <- reactiveVal(NULL);
    observeEvent({
      infotable()
      input$selectBiplotMetrics
      centering=input$selectCenteringType
      scaling=input$selectScalingMethod
    },{
      table1 <- infotable();
      if (is.null(table1)) {
        sregData(NULL)
        return(NULL)
      }
      ploteo <- GGEModel(table1,
                         SVP=input$selectBiplotMetrics,
                         centering=input$selectCenteringType,
                         scaling=input$selectScalingMethod)
      sregData(ploteo);
    })
    
    
    output$SREG1 <- renderPlot(GGEPlot(sregData()))
    output$SREG2 <- renderPlot(WhichWon(sregData()))
    output$SREG3 <- renderPlot(RankEnv(sregData()))
    output$SREG4 <- renderPlot(RankGen(sregData()))
      
    
    ## AMMI model 
    # output$AMMItableOut <- DT::renderDT({
    #   table1 <- infotable();
    #   if (is.null(table1)) {
    #     sregData(NULL)
    #     return(NULL)
    #   }
    #   preAMMITable<-AMMIprocessing(table1)
    #   print(preAMMITable)
    # })
      
    output$AMMIplot<- renderPlot({
      table1 <- infotable();
      if (is.null(table1)) {
        sregData(NULL)
        return(NULL)
      }
      preAMMITable<-AMMIprocessing(table1)
      print("preAMMITable calculated appropriatelly")
      print(preAMMITable)
      #colnames(preAMMITable)[1]<-"X"
      
      MeltedTable<-melt(preAMMITable,id = c("Genotype","Replicates"))
      print("MeltedTable calculated appropriatelly")
      ammi<-AMMI(ENV=MeltedTable$variable, GEN =MeltedTable$Genotype , Y=MeltedTable$value,
                  REP =MeltedTable$Replicates )
      plot(ammi, type =1, gcol="dodgerblue2",ecol="darkorange2")
    })
    
    output$AMMIresult<-DT::renderDT({
      table1 <- infotable();
      if (is.null(table1)) {
        sregData(NULL)
        return(NULL)
      }
      preAMMITable<-AMMIprocessing(table1)
      print("preAMMITable calculated appropriatelly")
      print(preAMMITable)
      #colnames(preAMMITable)[1]<-"X"
      
      MeltedTable<-melt(preAMMITable,id = c("Genotype","Replicates"))
      print("MeltedTable calculated appropriatelly")
      ammi<-AMMI(ENV=MeltedTable$variable, GEN =MeltedTable$Genotype , Y=MeltedTable$value,
                 REP =MeltedTable$Replicates )
      ANAVAR<-round(ammi$ANOVA, digits=2)
      print(ANAVAR) 
      })
      
    #Para obtener los ASV (AMMI Stability Value)
    output$AMMIasv<- DT::renderDT({
      table1 <- infotable();
      if (is.null(table1)) {
        sregData(NULL)
        return(NULL)
      }
      #print(table1)
      preAMMITable<-AMMIprocessing(table1)
      #print("preAMMITable calculated appropriatelly")
      print(preAMMITable)
      #preAMMITable<-data.frame(rownames(preAMMITable),preAMMITable)
      #colnames(preAMMITable)[1]<-"X"
      #print(preAMMITable)
      MeltedTable<-melt(preAMMITable,id = c("Genotype","Replicates"))
      #print("MeltedTable calculated appropriatelly")
      #print(MeltedTable)
      #browser()
      ammi<-AMMI(ENV=MeltedTable$variable, GEN =MeltedTable$Genotype , Y=MeltedTable$value,
                 REP =MeltedTable$Replicates )
      Idx<-index.AMMI(ammi)
      #print(Idx)
      # Crops with improved stability according AMMI.
      Idx<-round(Idx,digits=3)
      Idx[order(Idx[,3]),]
      })
    
    output$AMMItriplot<- renderPlot({
      table1 <- infotable();
      if (is.null(table1)) {
        sregData(NULL)
        return(NULL)
      }
      #print(table1)
      preAMMITable<-AMMIprocessing(table1)
      #print("preAMMITable calculated appropriatelly")
      print(preAMMITable)
      #preAMMITable<-data.frame(rownames(preAMMITable),preAMMITable)
      #colnames(preAMMITable)[1]<-"X"
      #print(preAMMITable)
      MeltedTable<-melt(preAMMITable,id = c("Genotype","Replicates"))
      print("MeltedTable calculated appropriatelly")
      #print(MeltedTable)
      #browser()
      ammi<-AMMI(ENV=MeltedTable$variable, GEN =MeltedTable$Genotype , Y=MeltedTable$value,
                 REP =MeltedTable$Replicates )
      if (requireNamespace("klaR", quietly = TRUE)) {
        plot(ammi,first=1,second=2,third=3, type=2,number=TRUE)
      }
    }) # Closing the Triplot
    #########################################################
    # Coinertia segment
    infotableRaw1 <- reactiveVal(NULL);
    infotable1 <- reactiveVal(NULL);
    infotableRaw2 <- reactiveVal(NULL);
    infotable2 <- reactiveVal(NULL);
    #browser()
    observeEvent({
      input$file1
      input$checkbox1
    }, {
      if (is.null(input$file1))
        return(NULL)
      
      infotableRaw1(read.table(input$file1$datapath, sep="\t", header=input$checkbox1))
      infotable1(read.table(input$file1$datapath, sep="\t", header=input$checkbox1,row.names = 1))
    })
    
    observeEvent({
      input$file2
      input$checkbox2
    }, {
      if (is.null(input$file2))
        return(NULL)
      
      infotableRaw2(read.table(input$file2$datapath, sep="\t", header=input$checkbox2))
      infotable2(read.table(input$file2$datapath, sep="\t", header=input$checkbox2,row.names = 1))
    })
    
    observeEvent({
      infotableRaw1()
    }, {
      table1 <- infotable1();
      if (!all(apply(table1,2,is.numeric))) {
        showModal(modalDialog("Warning: All your data must be numeric and it is not", title="Incorrect table format", size="l"))
        infotable1(NULL)
        return(NULL)
      }
      ## Dealing with missing values of Table1
      for(i in 1:ncol(table1)){
        table1[is.na(table1[,i]), i] <- mean(table1[,i], na.rm = TRUE)
      }
      infotable1(table1)
    })
      
      
    observeEvent({
        infotableRaw2()
      }, {
        table2 <- infotable2();
    
    if (!all(apply(table2,2,is.numeric))) {
      showModal(modalDialog("Warning: All your data must be numeric and it is not", title="Incorrect table format", size="l"))
      infotable2(NULL)
      return(NULL)
    }
    ## Dealing with missing values
    for(i in 1:ncol(table2)){
        table2[is.na(table2[,i]), i] <- mean(table2[,i], na.rm = TRUE)
      }
      infotable2(table2)
    })
    
      output$CoiGenotypes<- renderPlot({
      table1 <- infotable1();
      if (is.null(table1)) {
        #sregData(NULL)
        return(NULL)
      }
      table2 <- infotable2();
      if (is.null(table2)) {
        #sregData(NULL)
        return(NULL)
      }
      #browser()
      pcatable1<- dudi.pca(table1, scale = FALSE, scan = FALSE, nf = 3)
      pcatable2<- dudi.pca(table2, scale = FALSE, scan=FALSE, nf=3)
      #print(pcatable1)
      #print(pcatable2)
      coib <- coinertia( pcatable1,pcatable2, scann = FALSE)
      print (coib)
      #browser()
      plot(coib)
    }) # Closing the coinertia plot
    # plot(sepan(coib))
    # coib
    # summary(coib)
    # str(coib)
    # 
    # 
    # TCHT<-t(TCH)
    # PolPT<-t(PolP)
    # 
    # pcaPolPT<- dudi.pca(PolPT, scale = FALSE, scan = FALSE, nf = 2)
    # pcaTCHT<- dudi.pca(TCHT, scale=F, scan=F, nf=2)
    # 
    # coib <- coinertia(pcaPolPT, pcaTCHT, scann = FALSE)
    # plot(coib)
    # plot(sepan(coib))
    # coib
    # summary(coib)
  } # Closing the Function
)  # Closing shinyServer






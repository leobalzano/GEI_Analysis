#############################################
##########     GEIplatform.R     ############
#############################################
# Author: Leandro Balzano-Nogueira
# University of Florida

# This script is build to create a visualization of different 
# Genotype x Environment analyses under a visualization through Shiny
# And also it will deliver a Markdown report.

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
# library("agricolae", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
# library(DT)
# library("shiny", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
# library("shinyjs", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
# library("shinydashboard", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
# library("GGEBiplots", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
# library("GGEBiplotGUI", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
# library("reshape2", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
#############################################
# setwd("/Users/leobalzano/Dropbox (UFL)/Projects2018/GEIPlatform/GEIplatformShiny/")
# getwd()

#############################################
## ui.R ##
shinyUI(
  dashboardPage(skin = "blue", #skin color
   
    ### HEADER
    dashboardHeader(title="GEI analysis"), ## Header

    ### Side Bar
    dashboardSidebar(
      sidebarMenu(id = "sidebarmenu", ### create a menu in the sidebar
        # Input Sidebar Tab
        menuItem("Input", tabName = "input", icon = icon("code-branch")) ,
        conditionalPanel("input.sidebarmenu !== 'coinertia'",
          fluidPage(
            # File Input
            fileInput("file", label = h3("File input")),
            hr(),
  
            # Checkbox
            checkboxInput("checkbox", label = "Header", value = TRUE),
            hr()#,
          ) # Closing fluidPage of input  
        ), # Closing conditionalPanel of Input
        # PCA Sidebar Tab
        menuItem("PCA", tabName = "PCA", icon = icon("code-branch")),

        # SREG Sidebar Tab
        menuItem("SREG Model", tabName = "SREG", icon = icon("code-branch")),
        conditionalPanel("input.sidebarmenu === 'SREG'",
          fluidPage(
            # Select box of Biplot Metrics
            selectInput("selectBiplotMetrics", label = h6("Biplot Metrics"),
                        choices=list(
                          "JK Biplot (preserving row metrics)" = "row",
                          "GH Biplot (preserving column metrics)" = "column",
                          "HJ Biplot (Dual metrics)" = "dual",
                          "SQ Biplot (simmetrical representation)" = "symmetrical"),
                        selected = "dual"),
            hr(),
    
            # Select box of Centering Type
            selectInput("selectCenteringType", label = h6("Centering by Type"),
                        choices = list(
                          # "No centering" = "none",
                          "Global centering (G+E+GE)" = "global",
                          "Genotypical emphasis centering (G+GE)" = "tester",
                          "Double centering (GE)" = "double"),
                        selected = "tester"),
            hr(),
    
            # Select box of Biplot Metrics
            selectInput("selectScalingMethod", label = h6("Scaling Method"),
                        choices = list("No scaling" = "none", "Dividing by Std Deviation" = "sd"),
                        selected = "none"),
            hr()
    
          )# Closing fluidPage of second menuItem
      ),# Closing ConditionalPanel of SREG sidebar tab
      
      # AMMI Sidebar Tab
      menuItem("AMMI Model", tabName = "AMMItable", icon = icon("code-branch")),
      
      # Coinertia Analysis
      menuItem("Coinertia", tabName = "coinertia", icon = icon("code-branch")) ,
      conditionalPanel("input.sidebarmenu === 'coinertia'",
      fluidPage(
        # File Input 1
        fileInput("file1", label = h3("File input 1")),
        #hr(),
        
        # Checkbox 1
        checkboxInput("checkbox1", label = "Header", value = TRUE),
        #hr()#,
        # File Input 2
        fileInput("file2", label = h3("File input 2")),
        #hr(),
        
        # Checkbox 2
        checkboxInput("checkbox2", label = "Header", value = TRUE)
        #hr()#,
      )
      )
      )# Closing sidebarMenu
    ), #Closing dashboardSidebar
    
    ### BODY
    dashboardBody(
      tabItems(
        tabItem(tabName = "input",
          DT::DTOutput("value")
        ),
        
        tabItem(tabName = "PCA",
            fluidRow(box(plotOutput("PCA"))
            )
        ),
        
        tabItem(tabName = "SREG",
            fluidRow(box(plotOutput("SREG1")),
                     box(plotOutput("SREG2")),
                     box(plotOutput("SREG3")),
                     box(plotOutput("SREG4"))
            )
        ),
        
        tabItem(tabName = "AMMItable",
            fluidRow(#box (DT::DTOutput("AMMItableOut")),
                     box (plotOutput ("AMMIplot")),
                     box (plotOutput ("AMMItriplot")),
                     box (DT::DTOutput("AMMIresult")),
                     box (DT::DTOutput("AMMIasv"))
            ) #Closing the fluid Row of AMMI    
        ), # Closing tabItem AMMI  
        
        tabItem(tabName = "coinertia",
                fluidRow(plotOutput ("CoiGenotypes")
                ) #Closing the fluid Row of Coinertia
        ) # Clossing tabItem of Coinertia    
      ) # Closing TabItems
    ) # Closing dashboardBody
  ) # Closing dashboardPage
) # Closing shinyUI




#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(grid)
library(shinythemes)
library('shinyDirectoryInput')

# Define UI for application that analyses REMI-seq screen data
shinyUI(navbarPage(theme = shinytheme("cerulean"),"REMI-seq analyser",
                   tabPanel("Step 1: Experimental design",
                            fluidRow(
                              column(6,
                                     h4('1. Choose number of samples:'),
                                     numericInput('lgrid1', 'Number of samples (including the control)', 1,
                                                  min = 1, max = 50)),
                              
                              column(6,
                                     h4('2. Load files with experimental design'),
                                     fileInput('file1', 'Choose CSV File',accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))
                              )
                            ),
                            br(),
                            br(),
                            fluidRow(
                              column(12,
                                     h4('3. Select the reference sample:'),
                                     selectInput("sample5", "",choices= 1:3)
                              )
                            ),
                            br(),
                            br(),
                            fluidRow(
                              column(12,
                                  h4('4. Check if everything is correct'),
                                  DT::dataTableOutput("exp_design", width = "100%", height = "100%")
                              )
                            )
                   ),
                   tabPanel("Step 2: Extract tags",
                            fluidRow(
                              column(6,align="center",
                                     h4('1. Load the annotations file:'),
                                     textInput('file2', 'Choose CSV File',value="../help_files/annotated_positions_03022016")
                              ),
                              column(6,align="center",
                                     h4('2. Load the used indices:'),
                                     textInput('file3', 'Choose CSV File',value="../help_files/indices")
                              ) 
                            ),
                            br(),
                            br(),
                            fluidRow(
                              column(6,align="center",
                                      h4('3. Load the files:'),
                                      #textInput("directory", "Choose directory that contains fastq files:",value = "")
                              			 directoryInput('directory', label = 'Choose directory that contains fastq files:',value = '../raw_reads')
                              )
                            ),
                            br(),
                            hr(),
                            br(),
                            fluidRow(
                              column(5,align="center",
                                     actionButton("get_tags", label = "Start Analysis")),
                              column(2,align="center",
                                     h3('or')),
                              column(5,align="center",
                                     actionButton("load_results", label = "Load precomputed file"))
                            ),
                            fluidRow(
                              column(12,
                                     DT::dataTableOutput("stats", width = "100%", height = "100%")
                              )
                            )
                   ),
                   tabPanel("Step 3: Generate counts",
                            fluidRow(
                              column(8,
                                 h3('Choose the method to combine the counts:')
                              ),
                              column(4,align="center",
                                     h4('4. Choose cutoff:'),
                                     sliderInput("cutoff","", min = 0, max = 1000, value = 5)
                              )
                            ),
                            fluidRow(
                              column(12,offset=1,
                                     radioButtons("radio", label = "",
                                                  choices = list("Counts for each insertion point" = 1, "Counts for different vectors at each insertion point" = 2, "Counts for each gene" = 3,"Counts for each gene including 500 bp upstream" = 4,"Raw read counts per insertion point" = 5), 
                                                  selected = 1
                                      )
                              )
                            ),
                            br(),
                            hr(),
                            br(),
                            fluidRow(
                              column(12,
                                     DT::dataTableOutput("norm", width = "100%", height = "100%")
                              ),
                              column(12,
                                     downloadButton("downloadData", "Download")
                              )
                            )
                   ),
                   tabPanel("Step 4: Check replicates",
                            fluidRow(
                              column(12,
                                selectInput("sample", "Sample:",choices= 1:3),
                                checkboxInput("g0", "Remove low readcount mutants", value = FALSE, width = NULL),
                                checkboxInput("g1", "Remove mutants with at least 4x difference in the reference library", value = FALSE, width = NULL))
                            ),
                            br(),
                            hr(),
                            br(),
                            fluidRow(
                              column(12,
                              plotOutput("plot1")),
                              column(12,
                                     downloadButton("downloadPlots1", "Download Plots")
                              )
                            )
                            
                   ),
                   tabPanel("Step 5: Check for enriched and depleted mutants",
                   				 
                   				 
                   				
                            fluidRow(
                            	
                              column(12,
                                     h3('Select the sample that you want to plot:'),
                                     selectInput("samplePlot", "",choices= 1:3)
                              ),
                              column(12,
                              			 plotOutput("plot2")
                              )
                              
                              
                            )
                   ),
									 tabPanel("Step 6: Download results",
                            fluidRow(
                              column(12,
                                     DT::dataTableOutput("fold", width = "100%", height = "100%")
                              ),
                              column(12,
                                     DT::dataTableOutput("sumstats", width = "100%", height = "100%")
                              ),
                              column(12,
                              			 downloadButton("downloadData_FC", "Download table")
                              )
                            )
                   )
))
library(shiny)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
library(shinyBS)
#BiocManager::install("DESeq2")
library(DESeq2)
library(plotly)
library(gplots)
library(RColorBrewer)
library(stats)
library(grDevices)
library(eply) ##unquote
#install.packages("pheatmap")
library(pheatmap)#install.packages("flexclust") ### for dist2
library(flexclust) 
library(WGCNA)
#install.packages("genie") ## hclust2
library(genie)
library(flashClust,verbose=FALSE) ##genetree WGNA
library(DBI)
library(RSQLite)
library(Rtsne)
datapath = "/home/ubuntu/database/"   # production server
# 
sqlite  <- dbDriver("SQLite")
convert <- dbConnect( RSQLite::SQLite(), paste0(datapath, "convertIDs.db"))  #read only mode
keggSpeciesID = read.csv(paste0(datapath, "data_go/KEGG_Species_ID.csv"))
# List of GMT files in /gmt sub folder
gmtFiles = list.files(path = paste0(datapath,"pathwayDB"), pattern=".*\\.db")
gmtFiles = paste(datapath, "pathwayDB/", gmtFiles,sep="")
geneInfoFiles = list.files(path = paste0(datapath, "geneInfo"), pattern=".*GeneInfo\\.csv")
geneInfoFiles = paste(datapath, "geneInfo/", geneInfoFiles,sep="")

STRING10_species = read.csv(paste0(datapath, "data_go/STRING10_species.csv"))

###colors for heatmap
mycolors = sort(rainbow(20))[c(1,20,10,11,2,19,3,12,4,13,5,14,6,15,7,16,8,17,9,18)] # 20 colors for kNN clusters
#Each row of this matrix represents a color scheme;

hmcols <- colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF",
                                 "#E0F3F8", "#91BFDB", "#4575B4")))(75)
heatColors = rbind(  greenred(75),     bluered(75),     
                     colorpanel(75,"green", "black","magenta"),
                     colorpanel(75,"blue", "yellow","red"),hmcols )
rownames(heatColors) = c("Green-Black-Red", "Blue-White-Red", "Green-Black-Magenta",
                         "Blue-Yellow-Red", "Blue-white-brown")
colorChoices = setNames(1:dim(heatColors)[1], rownames(heatColors)) # for pull down menu

shinyUI(fluidPage(
  navbarPage(title = "RNA-SEQ WORKFLOW",
             tabPanel("Load Data", class = "active",
                      sidebarLayout(
                        sidebarPanel(
                          
                          radioButtons("action1",strong("1.Choose data type"),
                                       choices = list("Read counts data (recommended)" = 1), 
                                                      #"Normalized expression values (RNA-seq FPKM,microarray,etc.)" = 2),
                                       selected = 1),
                          radioButtons("symbols", "Select symbols", 
                                       c("File with Ensembl Gene IDs only"      = 1,            
                                         "File with Gene Symbols only"         = 2, 
                                         "File with Ensembl and Gene Symbols" = 3
                                         
                                       )),
                          fileInput("upfile",strong("2.Upload expression data(CSV or text)"),
                                    accept = c(
                                      'text/csv',
                                      'text/comma-separated-values',
                                      'text/tab-separated-values',
                                      'text/plain',
                                      '.csv',
                                      '.tsv'
                                    )##accept
                          ),##fileinput
                          
                          
                          uiOutput("choose_columns"),
                          
                          fileInput("metafile",strong("Upload an metadata file(CSV or text)"),
                                    accept = c(
                                      'text/csv',
                                      'text/comma-separated-values',
                                      'text/tab-separated-values',
                                      'text/plain',
                                      '.csv',
                                      '.tsv'          
                                    )##accept
                          ),##fileinput
                          
                          
                          uiOutput("choose_rows")
                          
                          ,strong("3. Verify guessed species. Change if neccessary.")
                          #,selectInput("selectOrg", label = NULL,"Best matching species",width='100%')
                          ,uiOutput("selectOrg")
                          
                          
                        ),##sidebarpanel
                        mainPanel(
                          tableOutput("contents"),
                          tableOutput("metadata"),
                          textOutput("texted")
                        ) ##mainpanel
                      ) ##sidebarLayout
             ), ##tabPanel
             
             tabPanel("Pre-process", calss="active",
                      sidebarLayout(
                        sidebarPanel(
                          h5("Examine the results of Scatterplot for selected comparison")
                          ,uiOutput("treatment")
                          ,uiOutput("treatment1"),
                          downloadButton('downloadNormalizeddata', 'Normalized Count data'),
                          sliderInput("range", "Range:",
                                      min = 1, max = 100,
                                      value = c(1,20))
                          
                          
                        ),
                        
                        
                        mainPanel(
                          tableOutput("transformed"),
                          plotOutput("barplot"),
                          plotOutput("boxplot"),
                          plotOutput("densityplot"),
                          plotOutput("scatterplot"),
                          plotOutput("dotplot")
                        ) # mainPanel
                      ) ##sidebarlayout
             ), #tabPanel
             
             tabPanel("Heatmap",
                      sidebarLayout(
                        sidebarPanel(
                          sliderInput("obs", "Number of observations:",
                                      min = 0, max = 12000, value = 50,step = 50),
                          downloadButton('downloadHeatmap', 'High-resolution figure'),
                          downloadButton('downloadHeatmapdata', 'Heatmap data'),
                          actionButton("showStaticHeatmap", "Interactive heatmap")
                        ),
                        
                        
                        mainPanel(
                          plotOutput("heatmap"),
                          bsModal("ModalExample","Interactive Heatmap","showStaticHeatmap",size = "large",
                                  plotlyOutput("heatmapPlotly",width = "600px", height = "800px"))
                        )#mainPanel
                      ) #sidebarlayout
             ),#tabPanel
             
             tabPanel("K-means",
                      sidebarLayout(
                        sidebarPanel(
                          sliderInput("variable_genes", "Number of observations:",
                                      min = 0, max = 12000, value = 50,step = 100),
                          sliderInput("clusters","Number of clusters",
                                      min = 2,max = 30,value=2),
                          actionButton("showElbowMethod", "Elbow Method")
                          ,HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />') 
                          ,h5("Pathway database")
                          ,htmlOutput("selectGO3")
                          ,tags$style(type='text/css', "#selectGO3 { width:100%;   margin-top:-9px}"),
                          sliderInput("redundant","Percentage ratio of redundancy",
                                      min = 0,max = 1,value=1,step = 0.01),
                          downloadButton('downloadkmeansdata', 'Enriched pathway data')
                          
                          
                        ),
                        mainPanel(
                          plotOutput("Kmeans"),
                          bsModal("ModalExample1","Elbow Method","showElbowMethod",size = "large",
                                  plotOutput("Elbowplot",width = "100%", height = "800px")),
                          tableOutput("results")
                        )#mainPanel
                      )#sidebarLayout
             ),#Tabpanel
             tabPanel("PCA",
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("PCA_MDS", "Methods", 
                                       c("Principal Component Analysis"      = 1,            
                                         "Multidimensional Scaling"         = 3, 
                                         "t-SNE"                            = 2
                                         
                                       ))
                        ),#sidebarpanel
                        
                        mainPanel(
                          plotOutput("PCA")
                        )
                        
                      )#sidebarlayout
             ),#TabPanel
             tabPanel("Differential Expression",
                      sidebarLayout(
                        sidebarPanel(
                          h5("Examine the results of DEGs for each comparison")
                          ,uiOutput("listComparisons")
                          ,uiOutput("listComparisons1")
                          ,br()
                          ,numericInput("num", 
                                        h3("Fold Change"), 
                                        value = 1)
                          ,numericInput("num1",
                                        h3("adj P-value"),
                                        value= 0.05),
                          
                          downloadButton('downloadDEG1data', 'DESeq2 output')
                        ),
                        
                        mainPanel(
                          tableOutput('text'),
                          tableOutput("table")
                        )##mainPanel
                        
                      )##sidebarLayout
                      
             ),#tabPanel
             
             tabPanel("Differential Pathways",
                      sidebarLayout(
                        sidebarPanel(
                          htmlOutput("selectGO5"),
                          sliderInput("redundant1","Percentage ratio of redundancy",
                                      min = 0,max = 1,value=1,step = 0.01),
                          downloadButton('downloadDEG2data', 'Enriched DEG pathway')
                        ),
                        
                        
                        mainPanel(
                          plotOutput("heatmap1"),
                          tableOutput("degenrichment")
                        )#mainPanel
                      ) #sidebarlayout
             ),#tabPanel
             
             tabPanel("WGCNA",
                      sidebarLayout(
                        sidebarPanel(
                          numericInput("GENES", 
                                       h3("Most Variable genes"), 
                                       value = 1000),
                          numericInput("soft", 
                                       h3("Soft Threshold"), 
                                       value = 3),
                          numericInput("Mod", 
                                       h3("Minimum Module Size"), 
                                       value = 10),
                          uiOutput("Modules"),
                          uiOutput("selectGO4"),
                          
                          numericInput("threshold", 
                                       h3("Edge Threshold"), 
                                       value = 0.3),
                          numericInput("topgenes", 
                                       h3("Top genes"), 
                                       value = 10),
                          downloadButton('downloadWGCNAdata', 'Enriched WGCNA pathway')
                        ),
                        
                        
                        mainPanel(
                          plotOutput("dendogram"),
                          plotOutput("network"),
                          tableOutput("enriched_pathways")
                        )#mainPanel
                      ) #sidebarlayout
             )#tabPanel
             
  ) ##navbarpage
)## fluidpage
)


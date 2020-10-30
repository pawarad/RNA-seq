options(shiny.maxRequestSize=30*1024^2) ##Shiny can only read files upto 5mb to exceed the limit to 30mb put this above server
##Server Code

shinyServer(function(input, output, session) {
  #options(shiny.maxRequestSize=30*1024^2)
  ##############----------------- read files --------------------------------------------
  library(DBI)
  
  datapath = "/home/ubuntu/database/"   # production server
  
  sqlite  <- dbDriver("SQLite")
  convert <- DBI::dbConnect(RSQLite::SQLite(), paste0(datapath, "convertIDs.db") )  #read only mode
  keggSpeciesID = read.csv(paste0(datapath, "data_go/KEGG_Species_ID.csv"))
  # List of GMT files in /gmt sub folder
  gmtFiles = list.files(path = paste0(datapath,"pathwayDB"), pattern=".*\\.db")
  gmtFiles = paste(datapath, "pathwayDB/", gmtFiles,sep="")
  geneInfoFiles = list.files(path = paste0(datapath, "geneInfo"), pattern=".*GeneInfo\\.csv")
  geneInfoFiles = paste(datapath, "geneInfo/", geneInfoFiles,sep="")
  
  STRING10_species = read.csv(paste0(datapath, "data_go/STRING10_species.csv"))
  
  #dl <- 61
  minFDR <- 0.05
  ##########---------------------- Functions ------------------------------------------
  ### first build ConvertID function
  cleanGeneSet <- function (x){
    # remove duplicate; upper case; remove special characters
    x <- unique( toupper( gsub("\n| ","",x) ) )
    x <- x[which( nchar(x)>1) ]  # genes should have at least two characters
    return(x)
  }
  
  add_legend <- function(...) {
    opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
                mar=c(0, 0, 0, 6), new=TRUE)
    on.exit(par(opar))
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend(...)
  }
  
  orgInfo <- dbGetQuery(convert, paste("select distinct * from orgInfo " ))
  orgInfo <- orgInfo[order(orgInfo$name),]
  speciesChoice <- setNames(as.list( orgInfo$id ), orgInfo$name2 )
  
  findSpeciesByIdName <- function (speciesID){ # find species name use id
    return( orgInfo[which(orgInfo$id == speciesID),3]  )
  }
  
  findSpeciesById <- function (speciesID){ # find species name use id
    return( orgInfo[which(orgInfo$id == speciesID),]  )
  }
  
  gmtCategory <- function(con){
    ix = grep(con$species[1,1],geneInfoFiles)
    pathway <- dbConnect(sqlite,gmtFiles[ix])
    #cat(paste("selectOrg:",selectOrg) )
    # Generate a list of geneset categories such as "GOBP", "KEGG" from file
    geneSetCategory <-  dbGetQuery(pathway, "select distinct * from categories " ) 
    geneSetCategory  <- sort( geneSetCategory[,1] )
    categoryChoices <- setNames(as.list( geneSetCategory ), geneSetCategory )
    categoryChoices <- append( setNames( "All","All available gene sets"), categoryChoices  )
    
    # move one element to the 2nd place
    move1 <- function(i) c(categoryChoices[1],categoryChoices[i],categoryChoices[-c(1,i)])
    i = which( names(categoryChoices)  == "KEGG"); categoryChoices= move1(i);	
    i = which( names(categoryChoices)  == "GOMF"); categoryChoices= move1(i);	
    i = which( names(categoryChoices)  == "GOCC"); categoryChoices= move1(i);	
    i = which( names(categoryChoices)  == "GOBP"); categoryChoices= move1(i);
    #change GOBP to the full description for display
    names(categoryChoices)[ match("GOBP",categoryChoices)  ] <- "GO Biological Process"
    names(categoryChoices)[ match("GOCC",categoryChoices)  ] <- "GO Cellular Component"
    names(categoryChoices)[ match("GOMF",categoryChoices)  ] <- "GO Molecular Function"
    
    dbDisconnect(pathway)
    return(categoryChoices )}
  
  convertID <- function(query,dl){
    querySet <- cleanGeneSet( unlist( strsplit( toupper(query),'\t| |\n|\\,')))
    # querySet is ensgene data for example, ENSG00000198888, ENSG00000198763, ENSG00000198804
    querySTMT <- paste( "select distinct id,ens,species from mapping where id IN ('", paste(querySet,collapse="', '"),"')",sep="")
    result <- dbGetQuery(convert, querySTMT)
    
    result <- result[which(result$species == dl ) ,]
    if( dim(result)[1] == 0  ) return(NULL) #stop("ID not recognized!")
    speciesMatched <- as.data.frame(paste("Using selected species ", findSpeciesByIdName(dl) )  )
    
    result <- result[which(!duplicated(result[,2]) ),] # remove duplicates in ensembl_gene_id
    result <- result[which(!duplicated(result[,1]) ),] # remove duplicates in user ID
    
    colnames(speciesMatched) = c("Matched Species (genes)" )
    conversionTable <- result[,1:2]; colnames(conversionTable) = c("User_input","ensembl_gene_id")
    conversionTable$Species = sapply(result[,3], findSpeciesByIdName )
    
    return(list(originalIDs = querySet,IDs=unique( result[,2]),
                species = findSpeciesById(result$species[1]),
                #idType = findIDtypeById(result$idType[1] ),
                speciesMatched = speciesMatched,
                conversionTable = conversionTable))}
  
  geneInfo <- function (con){
    if(is.null(con) ) return(as.data.frame("ID not recognized!") ) # no ID 
    querySet <- con$IDs
    if(length(querySet) == 0) return(as.data.frame("ID not recognized!") )
    ix = grep(con$species[1,1],geneInfoFiles)
    if(length(ix) == 1)  # if only one file           #WBGene0000001 some ensembl gene ids in lower case
    { x = read.csv(as.character(geneInfoFiles[ix]) ); x[,1]= toupper(x[,1]) } else # read in the chosen file 
    { return(as.data.frame("Multiple geneInfo file found!") )   }
    Set = match(x$ensembl_gene_id, querySet)
    Set[which(is.na(Set))]="Genome"
    Set[which(Set!="Genome")] ="List"
    #x = cbind(x,Set)  # just for debuging
    return(x)
  }
  
  detectGroups <- function (x){  # x are col names
    # Define sample groups based on column names
    # Args:
    #   x are vector of characters, column names in data file
    # Returns: 
    #   a character vector, representing sample groups.
    tem <- gsub("[0-9]*$","",x) # Remove all numbers from end
    #tem = gsub("_Rep|_rep|_REP","",tem)
    tem <- gsub("_$","",tem); # remove "_" from end
    tem <- gsub("_Rep$","",tem); # remove "_Rep" from end
    tem <- gsub("_rep$","",tem); # remove "_rep" from end
    tem <- gsub("_REP$","",tem)  # remove "_REP" from end
    return( tem )
  }
  
  FindOverlap <- function (converted,gInfo, GO) {
    Min_overlap <-2
    maxTerms =15 # max number of enriched termsx
    idNotRecognized = as.data.frame("ID not recognized!")
    
    # only coding
    gInfo <- gInfo[which( gInfo$gene_biotype == "protein_coding"),]  
    querySet <- intersect( converted$IDs, gInfo[,1]);
    
    ix = grep(converted$species[1,1],gmtFiles)
    totalGenes <- converted$species[1,7]
    
    if (length(ix) == 0 ) {return(idNotRecognized )}
    
    #pathway <- dbConnect(sqlite)
    PATH <- dbConnect(RSQLite::SQLite(),gmtFiles[ix])
    
    sqlQuery = paste("select distinct gene,pathwayID from pathway where gene IN ('", paste(querySet,collapse="', '"),"')" ,sep="")
    
    #cat(paste0("HH",GO,"HH") )
    #Go <- categoryChoices[3]
    sqlQuery = paste0(sqlQuery, " AND category ='",GO,"'")
    result <- dbGetQuery( PATH, sqlQuery  )
    if( dim(result)[1] ==0) {return(as.data.frame("No matching species or gene ID file!" )) }
    
    sharedGenesPrefered <- function(pathwayID) {
      tem <- result[which(result[,2]== pathwayID ),1]
      ix = match(tem, converted$conversionTable$ensembl_gene_id) # convert back to original
      tem2 <- unique( converted$conversionTable$User_input[ix] )
      if(length(unique(gInfo$symbol) )/dim(gInfo)[1] >.7  ) # if 70% genes has symbol in geneInfo
      { ix = match(tem, gInfo$ensembl_gene_id); 
      tem2 <- unique( gInfo$symbol[ix] )      }
      return( paste( tem2 ,collapse=" ",sep="") )}
    
    x0 = table(result$pathwayID)					
    x0 = as.data.frame( x0[which(x0>=Min_overlap)] )# remove low overlaps
    if(dim(x0)[1] <= 5 ) return(idNotRecognized) # no data
    colnames(x0)=c("pathwayID","overlap")
    pathwayInfo <- dbGetQuery( PATH, paste( " select distinct id,n,Description from pathwayInfo where id IN ('", 
                                            paste(x0$pathwayID,collapse="', '"),   "') ",sep="") )
    
    x = merge(x0,pathwayInfo, by.x='pathwayID', by.y='id')
    
    x$Pval=phyper(x$overlap-1,length(querySet),totalGenes - length(querySet),as.numeric(x$n), lower.tail=FALSE );
    x$FDR = p.adjust(x$Pval,method="fdr")
    x <- x[ order( x$FDR)  ,]  # sort according to FDR
    
    if(dim(x)[1] > maxTerms ) x = x[1:maxTerms,]	
    minFDR <- 0.05
    if(min(x$FDR) > minFDR) x=as.data.frame("No significant enrichment found!") else {
      x <- x[which(x$FDR < minFDR),] 
      
      x= cbind(x,sapply( x$pathwayID, sharedGenesPrefered ) )
      colnames(x)[7]= "Genes"
      x <- subset(x,select = c(FDR,overlap,n,description,Genes) )
      colnames(x) = c("Corrected P value (FDR)", "Genes in list", "Total genes in category","Functional Category","Genes"  )
      return(x)}
    
  }
  
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
  
  #####------------------------load data ---------------------------------------------- 
  
  output$choose_columns <- renderUI({
    # If missing input, return to avoid error later in function
    inFile <- input$upfile
    
    if(is.null(input$upfile))
      return()
    
    # Get the data set with the appropriate name
    ## condition to read csv or txt file
    df <- read.csv(inFile$datapath,row.names = NULL)
    if(dim(df)[2] <= 2 ) df <- read.table(inFile$datapath, sep="\t",header=TRUE,fill=TRUE)	# not CSV
    dat <- df
    colnames <- names(dat)
    
    # Create the checkboxes and select them all by default
    checkboxGroupInput("columns", "Choose columns", 
                       choices  = colnames,
                       selected = colnames)
  })
  
  df <- reactive({
    # If missing input, return to avoid error later in function
    inFile <- input$upfile
    
    if(is.null(input$upfile))
      return()
    
    # Get the data set
    df <- read.csv(inFile$datapath,row.names = NULL)
    if(dim(df)[2] <= 2 ) df <- read.table(inFile$datapath, sep="\t",header=TRUE,fill=TRUE,row.names =NULL)	# not CSV
    
    ##Remove duplicates and keep the highest value row
    # rs <- rowSums(df[,3:5])
    # df_ordered <- df[order(df[,3],rs,decreasing = TRUE), ]
    # df_ordered <- df_ordered[!duplicated(df_ordered[1]), ]
    
    
    dat <- df
    
    
    # Make sure columns are correct for data set (when data set changes, the
    # columns will initially be for the previous data set)
    if (is.null(input$columns) || !(input$columns %in% names(dat)))
      return()
    
    # Keep the selected columns
    dat <- dat[, input$columns, drop = FALSE]
    NSAF <<- dat
    # Return first 20 rows
    return(dat)
  })
  
  dis_data <- reactive({
    if (input$symbols == 1){
      query = df()[,1]
      converted <- convertID(query,dl())
      gInfo <- geneInfo(converted)
      
      d <- gInfo[,c(1,14)]
      x <- data.frame("ensembl_gene_id" = df()[,1])
      
      
      final_genes <- merge(x = x, y = d, by = "ensembl_gene_id", all.x = TRUE)
      df_ordered <- final_genes[!duplicated(final_genes[,1]), ]
      data <- cbind.data.frame(df_ordered,df()[,-1])
      
      ka <- grep('[A-Z]', data[,2], ignore.case = FALSE, perl = FALSE, value = FALSE,
                 fixed = FALSE, useBytes = FALSE, invert = TRUE)
      
      data[,1] <- as.character(data[,1])
      data[,2] <- as.character(data[,2])
      
      data[ka,2] <- data[ka,1]
      
      data <- data[!duplicated(data[,2]), ]
      data <- na.omit(data)
      xxx <- data ##read only first 20 as big files will take a lot of time to display
      return(xxx)
    }
    
    
    if (input$symbols == 2){
      data <- cbind.data.frame(row.names(df()),df())
      data <- data[!duplicated(data[,2]), ]
      data <- na.omit(data)
      xx <- data
      return(xx)
    }
    
    if (input$symbols == 3){
      data <- df()
      ka <- grep('[A-Z]', data[,2], ignore.case = FALSE, perl = FALSE, value = FALSE,
                 fixed = FALSE, useBytes = FALSE, invert = TRUE)
      
      data[,1] <- as.character(data[,1])
      data[,2] <- as.character(data[,2])
      
      data[ka,2] <- data[ka,1]
      data <- data[!duplicated(data[,2]), ]
      data <- na.omit(data)
      
      xxx <- data ##read only first 20 as big files will take a lot of time to display
      return(xxx)
      
    }
  })
  
  output$contents <- renderTable({
    dis <<- dis_data()
    dis_data()[1:20,]
  },
  ## to display data asthetically
  include.rownames=FALSE,striped=TRUE,bordered = TRUE, width = "auto",hover=T)
  
  
  
  output$choose_rows <- renderUI({
    
    inFile <- input$metafile
    # If missing input, return to avoid error later in function
    if(is.null(input$metafile))
      return()
    ## condition to read csv or txt file
    df_meta <- read.csv(inFile$datapath,row.names = NULL)
    
    if(dim(df_meta)[2] <= 2 ) df_meta <- read.table(inFile$datapath, sep="\t",header=TRUE,fill=TRUE,row.names = NULL)	# not CSV
    rownames(df_meta) <- df_meta[,1]
    df_meta[, 3] <- as.factor(df_meta[, 3]) ## for df_transformed
    
    # Get the data set with the appropriate name
    
    dat_meta <- df_meta
    rownames <- rownames(dat_meta)
    
    # Create the checkboxes and select them all by default
    checkboxGroupInput("rows", "Choose rows", 
                       choices  = rownames,
                       selected = rownames)
  })
  
  df_meta <- reactive({
    inFile <- input$metafile
    # If missing input, return to avoid error later in function
    if(is.null(input$metafile))
      return()
    ## condition to read csv or txt file
    df_meta <- read.csv(inFile$datapath)
    
    if(dim(df_meta)[2] <= 2 ) df_meta <- read.table(inFile$datapath, sep="\t",header=TRUE,fill=TRUE)	# not CSV
    rownames(df_meta) <- df_meta[,1]
    df_meta[, 3] <- as.factor(df_meta[, 3]) ## for df_transformed
    
    # Get the data set with the appropriate name
    
    dat_meta <- df_meta
    
    # Make sure columns are correct for data set (when data set changes, the
    # columns will initially be for the previous data set)
    if (is.null(input$rows) || !(input$rows %in% rownames(dat_meta)))
      return()
    
    # Keep the selected columns
    dat_meta <- dat_meta[input$rows, ,drop = FALSE]
    
    # Return first 20 rows
    return(dat_meta)
    
  })
  
  output$metadata <- renderTable({
    df_meta()
    
  },
  ## to display data asthetically
  include.rownames=FALSE,striped=TRUE,bordered = TRUE, width = "auto",hover=T)
  
  
  output$selectOrg <- renderUI({
    
    orgInfo <- dbGetQuery(convert, paste("select distinct * from orgInfo " ))
    orgInfo <- orgInfo[order(orgInfo$name),]
    speciesChoice <- setNames(as.list( orgInfo$id ), orgInfo$name2 )
    s <<- speciesChoice
    
    # Create the checkboxes and select them all by default
    selectInput("rowsss", "selectOrg", 
                choices  = speciesChoice,
                selected = speciesChoice )
  })
  
  dl <- reactive({
    orgInfo <- dbGetQuery(convert, paste("select distinct * from orgInfo " ))
    orgInfo <- orgInfo[order(orgInfo$name),]
    speciesChoice <- setNames(as.list( orgInfo$id ), orgInfo$name2 )
    d <- input$rowsss
    dl <- as.numeric(d)
    dl
  })
  
  
  output$texted <- renderPrint({
    dl()
    
    #speciesChoice$d
  })
  
  ########----------------------- pre-process --------------------------------------------
  df_transformed <- reactive({
    
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(dis_data()[,-c(1,2)]), colData = df_meta(), design = ~ Treatment)
    dds <- estimateSizeFactors(dds)
    df_transformed <- log2(counts(dds, normalized=TRUE)+4)
    return(df_transformed)
  })
  
  df_transformed_data <- reactive({
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(dis_data()[,-c(1,2)]), colData = df_meta(), design = ~ Treatment)
    dds <- estimateSizeFactors(dds)
    df_transformed <- log2(counts(dds, normalized=TRUE)+4)
    df_log <- cbind.data.frame(dis_data()[,c(1,2)],df_transformed)
    #colnames(df_log)[1] <- colnames(data_ui[1])
    df_log
  })
  
  output$downloadNormalizeddata <- downloadHandler(
    filename = "Normalized_counts.csv",
    
    content = function(file) {
      
      
      write.csv(df_transformed_data(),file = file,row.names = FALSE)
      
      #dev.off()
      #dev.copy2pdf(file)
    }
  )
  
  output$transformed <-  renderTable({
    if (input$symbols == 1){
      # query = df()[,1]
      # converted <- convertID(query,dl())
      # gInfo <- geneInfo(converted)
      # 
      # d <- gInfo[,c(1,14)]
      # x <- data.frame("ensembl_gene_id" = df()[,1])
      # 
      # 
      # final_genes <- merge(x = x, y = d, by = "ensembl_gene_id", all.x = TRUE)
      # df_ordered <- final_genes[!duplicated(final_genes[,1]), ]
      # 
      # data <- cbind.data.frame(df_ordered,df()[,-1])
      # data_ui <- data##read only first 20 as big files will take a lot of time to display
      # #return(xxx)
      df_log <- cbind.data.frame(dis_data()[,c(1,2)],df_transformed())
      #colnames(df_log)[1] <- colnames(data_ui[1])
      return(df_log[1:5,])
    }
    
    if (input$symbols == 2){
      df_log <- cbind.data.frame(dis_data()[,c(1,2)],df_transformed())
      return(df_log[1:5,])
      return(df_transformed()[1:5,])
      # data_ui <- df()[1:20,]
      # df_log <- cbind.data.frame(data_ui[,1],df_transformed())
      # colnames(df_log)[1] <- colnames(data_ui[1])
      # return(df_log[1:5,])
    }
    
    if (input$symbols == 3){
      df_log <- cbind.data.frame(dis_data()[,c(1,2)],df_transformed())
      return(df_log[1:5,])
    }
    
  },
  ## to display data asthetically
  include.rownames=FALSE,striped=TRUE,bordered = TRUE, width = "auto",hover=T)
  
  
  output$barplot <- renderPlot({
    #x <- as.matrix(as.data.frame(lapply(df_transformed(), as.numeric))) ##converting into numeric droping first column Ensemble
    x <- df_transformed()
    groups = as.factor(colnames(x ) )
    if(nlevels(groups)<=1 | nlevels(groups) >20 )
      col1 = 'green'  else
        col1 = rainbow(nlevels(groups))[ groups ]
    
    barplot(colSums(x)/1e6,
            col=col1,las=3, main="Total read counts (millions)")
  })
  
  output$boxplot <- renderPlot({
    # Box plot
    x = df_transformed() ##droping first coulmn Ensemble
    groups = as.factor(colnames(x ) )
    if(nlevels(groups)<=1 | nlevels(groups) >20 )
      col1 = 'green'  else
        col1 = rainbow(nlevels(groups))[ groups ]
    boxplot(x, las = 2, col=col1,
            ylab='Transformed expression levels',
            main='Distribution of transformed data')
    
  })
  
  output$treatment <- renderUI({
    
    selectInput("treatment", h3("select treatment"), 
                choices = colnames(df_transformed()))
  })
  
  output$treatment1 <- renderUI({
    
    selectInput("treatment1", h3("select second treatment"), 
                choices = colnames(df_transformed()))
  })
  
  output$scatterplot <- renderPlot({
    # Scatter plot of the first two samples
    #plot(df_transformed()[,1:2],xlab=colnames(df_transformed())[1],ylab=colnames(df_transformed())[2],
    #main='Scatter plot of first two samples')
    
    plot(df_transformed()[,c(input$treatment,input$treatment1)] ,xlab = input$treatment, ylab = input$treatment1 ,main='Scatter plot of first two samples') 
  })
  
  
  
  output$densityplot <- renderPlot({
    x <- df_transformed()
    groups = as.factor(colnames(x ) )
    if(nlevels(groups)<=1 | nlevels(groups) >20 )
      col1 = 'green'  else
        col1 = rainbow(nlevels(groups))[ groups ]
    maxDensity = max( apply(x,2, function(y) max(density(y)$y ) ) )	
    plot(density(x[,1]),col = col1[1], lwd=2,
         xlab="Expression values", main= paste("Density plot of transformed data"),
         ylim=c(0, maxDensity+0.01 )  )  #ylim=c(0,1)
    
    for( i in 2:dim(x)[2] )
      lines(density(x[,i]),col=col1[i],  lwd=1 )
    if(nlevels(groups)< 31 ) # if too many samples do not show legends
      legend("topright", levels(groups), lty=rep(1,nlevels(groups)), col=rainbow(nlevels(groups)) )	
    
  })
  
  output$dotplot <- renderPlot({
    dis <- dis_data()
    rownames(dis) <- dis[,2]
    rs <- rowSums(dis[,3:ncol(dis)])
    total <- sum(rs)
    percent <- (rs/total)*100
    rs_order <- as.data.frame(percent[order(percent,decreasing = TRUE)])
    
    dotchart(rs_order[input$range[1]:input$range[2],], labels = rownames(rs_order),
             cex = 0.6, xlab = "counts")
  })
  
  #######---------------------------------PCA------------------------------------------------------------------------
  
  output$PCA <-  renderPlot({
    if (input$PCA_MDS == 1){
      dds <- DESeqDataSetFromMatrix(countData = as.matrix(dis_data()[,-c(1,2)]), colData = df_meta(), design = ~ Treatment)
      dds <- estimateSizeFactors(dds)
      vsd <- vst(dds, blind=FALSE)
      plotPCA(vsd, intgroup=c("Treatment","Study_design"))
      #dev.off()
    }
    
    if (input$PCA_MDS == 2){
      dds <- DESeqDataSetFromMatrix(countData = as.matrix(dis_data()[,-c(1,2)]), colData = df_meta(), design = ~ Treatment)
      dds <- estimateSizeFactors(dds)
      vsd <- vst(dds, blind=FALSE)
      colData(vsd)
      sampleDists <- dist(t(assay(vsd)))
      sampleDistMatrix <- as.matrix( sampleDists )
      mds <- as.data.frame(colData(vsd))  %>%
        cbind(cmdscale(sampleDistMatrix))
      ggplot(mds, aes(x = `1`, y = `2`, color = Treatment, shape = Study_design)) +
        geom_point(size = 3) + coord_fixed()+scale_shape_manual(values=seq(0,30))
      
    }
    
    if (input$PCA_MDS == 3){
      groups = as.character (detectGroups( colnames(dis_data()[,-c(1,2)]) ) )
      g = unique(groups)# order is reversed
      
      # check for replicates, removes samples without replicates
      reps = as.matrix(table(groups)) # number of replicates per biological sample
      if ( sum( reps[,1] >= 2) <2 ) # if less than 2 samples with replicates
        return( list(results= NULL, comparisons = NULL, Exp.type="Failed to parse sample names to define groups. 
		Cannot perform DEGs and pathway analysis. Please double check column names! Use WT_Rep1, WT_Rep2 etc. ", topGenes=NULL)) 
      
      # remove samples without replicates
      g <- rownames(reps)[which(reps[,1] >1)]
      ix <- which(groups %in% g)  
      groups <- groups[ix]   
      rawCounts <- dis_data()[,-c(1,2)][,ix] 
      tsne1 <- Rtsne(t(rawCounts),dims = 2, perplexity = 1,max_iter = 300)
      plot(tsne1$Y)
    }
    
  })
  
  
  
  
  #######-------------------------------Heatmap---------------------------------------------------------------------------------------------
  heatmap <- reactive({
    
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(dis_data()[,-c(1,2)]), colData = df_meta(), design = ~ Treatment)
    # print(ncol(df()))
    # print(nrow(df_meta()))
    
    dds <- estimateSizeFactors(dds)
    df_transformed <- log2(counts(dds, normalized=TRUE)+4)
    rownames(df_transformed) <- dis_data()[,2]
    var_genes <- apply(df_transformed, 1, var)
    
    vargenes <<- var_genes
    #print(head(var_genes))
    select_var <- names(sort(var_genes, decreasing=TRUE))[1:input$obs]
    #print(head(select_var))
    highly_variable_lcpm <- df_transformed[select_var,]
    #print(dim(highly_variable_lcpm))
    as.matrix(highly_variable_lcpm)
  })
  
  plotHeatmap <- reactive({
    par(mar = c(1,1,1,1))
    lmat = rbind(c(0,4),c(0,1),c(3,2),c(5,0))
    lwid = c(2,6) # width of gene tree; width of heatmap
    lhei = c(1.5,.2,8,1.1)
    groups = detectGroups(colnames(heatmap()) )
    groups.colors = rainbow(length(unique(groups) ) )
    
    # Plot the heatmap
    heatmap.2(heatmap(),col=heatColors[as.integer(2),],trace="none", main=paste("Top",input$obs,"most variable genes across samples",sep= " "),scale="row"
              ,cexRow=0.75,margins=c(8,12)
              ,ColSideColors = groups.colors[ as.factor(groups)]
              ,srtCol=45
              ,cexCol = 1,
              lhei =c(1.5,5),lwid = c(2,4))
    
    
    if(length(unique(groups) ) <= 30 ) {  # only add legend when there is less categories
      par(lend = 1)           # square line ends for the color legend
      add_legend("bottomright",
                 legend = unique(groups), # category labels
                 col = groups.colors[ unique(as.factor(groups))],  # color key
                 lty= 1,             # line style
                 lwd = 10            # line width
      )}
    
  })
  
  output$heatmap <- renderPlot({
    
    plotHeatmap()
    
  })
  
  heatmap_download <- reactive({
    par(mar = c(1,1,1,1))
    lmat = rbind(c(0,4),c(0,1),c(3,2),c(5,0))
    lwid = c(2,6) # width of gene tree; width of heatmap
    lhei = c(1.5,.2,8,1.1)
    groups = detectGroups(colnames(heatmap()) )
    groups.colors = rainbow(length(unique(groups) ) )
    
    # Plot the heatmap
    heatmap_plot <- heatmap.2(heatmap(),col=heatColors[as.integer(2),],trace="none", main=paste("Top",input$obs,"most variable genes across samples",sep= " "),scale="row"
                              ,cexRow=0.75,margins=c(8,12)
                              ,ColSideColors = groups.colors[ as.factor(groups)]
                              ,srtCol=45
                              ,cexCol = 1,
                              lhei =c(1.5,5),lwid = c(2,4))
    heatmap_plot
  })
  output$downloadHeatmap <- downloadHandler(
    
    
    filename = "plot.pdf",
    
    content = function(filename) {
      
      pdf(file = filename,width = 10,height = 8)
      heatmap_download()
      dev.off()
      #dev.copy2pdf(file)
    }
  )
  
  heatmap_data <- reactive({
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(dis_data()[,-c(1,2)]), colData = df_meta(), design = ~ Treatment)
    # print(ncol(df()))
    # print(nrow(df_meta()))
    
    dds <- estimateSizeFactors(dds)
    df_transformed <- log2(counts(dds, normalized=TRUE)+4)
    rownames(df_transformed) <- dis_data()[,2]
    var_genes <- apply(df_transformed, 1, var)
    #vargenes <<- var_genes
    select_var <- names(sort(var_genes, decreasing=TRUE))[1:input$obs]
    highly_variable_lcpm <- df_transformed[select_var,]
    
    rownames(df_transformed) <- dis_data()[,1]
    var_genes <- apply(df_transformed, 1, var)
    select_var <- names(sort(var_genes, decreasing=TRUE))[1:input$obs]
    highly_variable_lcpm1 <- df_transformed[select_var,]
    #print(dim(highly_variable_lcpm))
    heatmap_data <- cbind.data.frame(rownames(highly_variable_lcpm1),rownames(highly_variable_lcpm),highly_variable_lcpm)
    colnames(heatmap_data[c(1,2),]) <- colnames(dis_data()[c(1,2),])
    heatmap_data
    #as.matrix(heatmap_data)
  })
  
  output$downloadHeatmapdata <- downloadHandler(
    filename = "Heatmap.csv",
    
    content = function(file) {
      
      
      write.csv(heatmap_data(),file = file,row.names = FALSE)
      
      #dev.off()
      #dev.copy2pdf(file)
    }
  )
  #p_heatmap_plotly <- eventReactive(input$showStaticHeatmap,
  #plot_ly(x=colnames(heatmap()), y=rownames(heatmap()),z =heatmap(),type = "heatmap",colors = c("red","white","blue")))
  
  output$heatmapPlotly <- renderPlotly(
    {
      groups = detectGroups(colnames(heatmap()) )
      groups.colors = rainbow(length(unique(groups) ) )
      
      plot_ly(x=colnames(heatmap()), y=rownames(heatmap()),z =heatmap(),type = "heatmap",colors = c("red","white","blue"))%>% 
        layout(margin = list(b = 150,l=200))
      
    }
  )
  
  
  
  ##############-------------------------------------K-means----------------------------------------------
  Kmeans_data <- reactive({
    df_log <- cbind.data.frame(dis_data()[,2],df_transformed())
    colnames(df_log)[1] <- colnames(dis_data()[2])
    df_log <- df_log[!duplicated(df_log[,1]), ]
    row.names(df_log) <- df_log[,1]
    df_log[,1] <- NULL
    var_genes <- apply(df_log, 1, var)
    select_var <- names(sort(var_genes, decreasing=TRUE))
    
    highly_variable_lcpm <- df_log[select_var,]
    
    x <- highly_variable_lcpm
    
    
    #x<- x1
    x = x[1:input$variable_genes,]
    x = 100* x / apply(x,1,function(y) sum(abs(y))) ###L1 Normalization
    
    x
    
  })
  
  kmeans_df <- reactive({
    x <- Kmeans_data()
    k <- input$clusters
    
    cl = kmeans(x,k,iter.max = 50)
    
    hc <- hclust2(dist(cl$centers-(apply(cl$centers,1,mean))))  #euclidean distance
    tem = match(cl$cluster,hc$order) #  new order
    x = x[order(tem),] ; 	bar = sort(tem)
    kmeans_df <- list( x = x, bar = bar)
    return(kmeans_df)
  })
  
  output$Kmeans <- renderPlot({
    x <- Kmeans_data()
    k <- input$clusters
    
    cl = kmeans(x,k,iter.max = 50)
    
    hc <- hclust2(dist(cl$centers-(apply(cl$centers,1,mean))))  #euclidean distance
    tem = match(cl$cluster,hc$order) #  new order
    x = x[order(tem),] ; 	bar = sort(tem)
    #kmeans_df <<- list( x = x, bar = bar)
    mycolor = 1
    
    heatmap.2(as.matrix(x),  Rowv =F,Colv=F, dendrogram ="none",
              col=heatColors[as.integer(mycolor),], density.info="none", trace="none", scale="none", keysize=.3
              ,key=F, labRow = F
              ,RowSideColors = mycolors[bar]
              ,margins = c(8, 24)
              ,srtCol=45)
    
    ngenes = as.character( table(bar))
    sideColors = mycolors
    
    legend.text = paste("Cluster ", toupper(letters)[unique(bar)], " (N=", ngenes,")", sep="") 
    
    par(lend = 1)           # square line ends for the color legend
    legend("topright",      # location of the legend on the heatmap plot
           legend = legend.text, # category labels
           col = sideColors,  # color key
           lty= 1,             # line style
           lwd = 10 )
  })
  
  output$Elbowplot <- renderPlot({
    x <- Kmeans_data()
    k = 30
    wss <- (nrow(x)-1)*sum(apply(x,2,var))
    for (i in 2:k) wss[i] <- sum(kmeans(x,centers=i,iter.max = 30)$withinss)
    par(mar=c(4,5,4,4))
    plot(1:k, wss, type="b", xlab="Number of Clusters (k)",
         ylab="Within groups sum of squares",
         cex=2,cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2	,xaxt="n"	 )
    axis(1, at = seq(1, 30, by = 2),cex.axis=1.5,cex=1.5)
    
    
  })
  
  output$selectGO3 <- renderUI({
    #query = rownames(kmeans_df()$x)[which(Kmeans()$bar == 1)]
    query = rownames(kmeans_df()$x)[which(kmeans_df()$bar == 1)]
    converted <- convertID(query,dl())
    gInfo <- geneInfo(converted)
    categoryChoices <- gmtCategory(converted)
    
    selectInput("selectGO3", label = NULL, # h6("Funtional Category"), 
                choices = categoryChoices[-1]) ## because we dont want all selected gene from caterogory choices  
  })
  
  output$treatment <- renderUI({
    
    selectInput("treatment", h3("select treatment"), 
                choices = colnames(df_transformed()))
  })
  
  
  
  Result <- reactive({
    results <- data.frame()
    for(i in 1:input$clusters){
      #query = rownames(kmeans_df()$x)[which(Kmeans()$bar == 1)]
      query = rownames(kmeans_df()$x)[which(kmeans_df()$bar == i)]
      converted <- convertID(query,dl())
      gInfo <- geneInfo(converted)
      categoryChoices <- gmtCategory(converted)
      result <- FindOverlap(converted,gInfo,input$selectGO3)
      result$direction = toupper(letters)[i] 
      #results <- result
      results <- rbind(results,result)
      #x <- results
      n = nrow(results)
      tem=rep(TRUE,n)
      
      geneLists = lapply(results$Genes, function(y) unlist( strsplit(as.character(y)," " )   ) )
      for( i in 2:n)
        for( j in 1:(i-1) ) { 
          if(tem[[j]]) { # skip if this one is already removed
            commonGenes = length(intersect(geneLists[[i]] ,geneLists[[j]] ) )
            if( commonGenes/ length(geneLists[[j]] ) > input$redundant )
              tem[[i]] = FALSE	
          }			
        }								
      results <- results[which(tem),]
    }
    
    minFDR  <- 0.05
    results= results[,c(6,1,2,4)]
    colnames(results)= c("Cluster","FDR","nGenes","Pathways")
    if(min(results$FDR) > minFDR ) results = as.data.frame("No signficant enrichment found.") else
      results= results[which(results$FDR < minFDR),]
    colnames(results)[2] = "adj.Pval"
    return(results)
  })
  
  output$results <- renderTable({
    Result()
  })
  
  Enriched_pathway_data <- reactive({
    results <- data.frame()
    for(i in 1:input$clusters){
      #query = rownames(kmeans_df()$x)[which(Kmeans()$bar == 1)]
      query = rownames(kmeans_df()$x)[which(kmeans_df()$bar == i)]
      converted <- convertID(query,dl())
      gInfo <- geneInfo(converted)
      categoryChoices <- gmtCategory(converted)
      result <- FindOverlap(converted,gInfo,input$selectGO3)
      result$direction = toupper(letters)[i] 
      #results <- result
      results <- rbind(results,result)
      #x <- results
      n = nrow(results)
      tem=rep(TRUE,n)
      
      geneLists = lapply(results$Genes, function(y) unlist( strsplit(as.character(y)," " )   ) )
      for( i in 2:n)
        for( j in 1:(i-1) ) { 
          if(tem[[j]]) { # skip if this one is already removed
            commonGenes = length(intersect(geneLists[[i]] ,geneLists[[j]] ) )
            if( commonGenes/ length(geneLists[[j]] ) > input$redundant )
              tem[[i]] = FALSE	
          }			
        }								
      results <- results[which(tem),]
    }
    
    minFDR  <- 0.05
    results= results[,c(6,1,2,4,5)]
    colnames(results)= c("Cluster","FDR","nGenes","Pathways","Genes")
    if(min(results$FDR) > minFDR ) results = as.data.frame("No signficant enrichment found.") else
      results= results[which(results$FDR < minFDR),]
    colnames(results)[2] = "adj.Pval"
    return(results)
  })
  
  output$downloadkmeansdata <- downloadHandler(
    filename = "Enriched pathway.csv",
    
    content = function(file) {
      
      write.csv(Enriched_pathway_data(),file = file,row.names = FALSE)
      
      #dev.off()
      #dev.copy2pdf(file)
    }
  )
  ###########-----------------------------DEG1 -----------------------------------------------------------
  
  listComparisons <- reactive({
    detectGroups <- function (x){  # x are col names
      # Define sample groups based on column names
      # Args:
      #   x are vector of characters, column names in data file
      # Returns: 
      #   a character vector, representing sample groups.
      tem <- gsub("[0-9]*$","",x) # Remove all numbers from end
      #tem = gsub("_Rep|_rep|_REP","",tem)
      tem <- gsub("_$","",tem); # remove "_" from end
      tem <- gsub("_Rep$","",tem); # remove "_Rep" from end
      tem <- gsub("_rep$","",tem); # remove "_rep" from end
      tem <- gsub("_REP$","",tem)  # remove "_REP" from end
      return( tem )
    }
    
    groups = as.character(detectGroups(colnames(dis_data()[,-c(1,2)])))
    g = unique(groups) # order is reversed
    
    groups1 = as.character ( detectGroups(df_meta()[,3]))
    g1 = unique(groups1)# order is reversed
    
    table_grp <- cbind(g,g1)
    table_grp
  })
  
  output$listComparisons <- renderUI({
    
    selectInput("select", h3("Select box"), 
                choices = listComparisons()[,1])
  })
  
  output$listComparisons1 <- renderUI({
    
    selectInput("select1", h3("Select box"), 
                choices = listComparisons()[,1],
                selected = listComparisons()[2,1])
  })
  
  
  dds <- reactive({ 
    dds <- DESeqDataSetFromMatrix(countData = dis_data()[,-c(1,2)], colData = df_meta(), design = ~ Treatment)
    
    dds <- DESeq(dds)
    dds
  })
  
  deseq2_results <- reactive({
    trt <- input$select
    #which(table_grp[,1] == trt)
    ind = which(listComparisons()[,1] == trt)
    t1 <- listComparisons()[ind,2]
    print(t1)
    
    trt1 <- input$select1
    ind1 = which(listComparisons()[,1] == trt1)
    t2 <- listComparisons()[ind1,2]
    print(t2)
    
    pair <- paste0(trt,'-',trt1)
    print(pair)
    
    deseq2_results <- results(dds(),contrast = c("Treatment", t1 ,t2))
    deseq2_results <- cbind.data.frame(dis_data()[,2],deseq2_results)
    colnames(deseq2_results)[1] <- colnames(df()[1])
    deseq2_results
  })
  
  output$text <- renderTable({
    
    trt <- input$select
    #which(table_grp[,1] == trt)
    ind = which(listComparisons()[,1] == trt)
    t1 <- listComparisons()[ind,2]
    print(t1)
    
    trt1 <- input$select1
    ind1 = which(listComparisons()[,1] == trt1)
    t2 <- listComparisons()[ind1,2]
    print(t2)
    
    pair <- paste0(trt,'-',trt1)
    print(pair)
    
    
    resFC <- subset(deseq2_results() ,log2FoldChange < input$num)
    #resFC <- subset(deseq2_results ,log2FoldChange < 2)
    resSig <- subset(resFC, padj < input$num1)
    
    
    Up <- sum(resSig$log2FoldChange > 0, na.rm=TRUE)
    Uregulated <<- rownames(subset(resSig,log2FoldChange > 0))
    
    Down <- sum(resSig$log2FoldChange < 0, na.rm=TRUE)
    Doregulated <<- rownames(subset(resSig,log2FoldChange < 0))
    
    tabel <- cbind(pair,Up,Down)
    tabel
  })
  
  output$table <- renderTable({
    #dese_results <- cbind.data.frame(rownames(deseq2_results()),deseq2_results())
    #print(colnames(deseq2_results()[1]))
    
    trt <- input$select
    #which(table_grp[,1] == trt)
    ind = which(listComparisons()[,1] == trt)
    t1 <- listComparisons()[ind,2]
    print(t1)
    
    trt1 <- input$select1
    ind1 = which(listComparisons()[,1] == trt1)
    t2 <- listComparisons()[ind1,2]
    print(t2)
    
    pair <- paste0(trt,'-',trt1)
    print(pair)
    
    #colnames(deseq2_results()[1]) <- colnames(df()[1])
    resFC <- subset(deseq2_results() ,log2FoldChange < input$num)
    #resFC <- subset(deseq2_results ,log2FoldChange < 2)
    resSig <- subset(resFC, padj < input$num1)
    resSig <- resSig[order(resSig$log2FoldChange, decreasing = TRUE), ]
    
    
    #colnames(resSig[1]) <- colnames(df()[1])
    resSig[1:15,]
    #rownames(deseq2_results())
  })
  
  DESeq2_output <- reactive({
    #dese_results <- cbind.data.frame(rownames(deseq2_results()),deseq2_results())
    #print(colnames(deseq2_results()[1]))
    
    trt <- input$select
    #which(table_grp[,1] == trt)
    ind = which(listComparisons()[,1] == trt)
    t1 <- listComparisons()[ind,2]
    print(t1)
    
    trt1 <- input$select1
    ind1 = which(listComparisons()[,1] == trt1)
    t2 <- listComparisons()[ind1,2]
    print(t2)
    
    pair <- paste0(trt,'-',trt1)
    print(pair)
    
    #colnames(deseq2_results()[1]) <- colnames(df()[1])
    resFC <- subset(deseq2_results() ,log2FoldChange < input$num)
    #resFC <- subset(deseq2_results ,log2FoldChange < 2)
    resSig <- subset(resFC, padj < input$num1)
    resSig <- resSig[order(resSig$log2FoldChange, decreasing = TRUE), ]
    
    
    #colnames(resSig[1]) <- colnames(df()[1])
    return(resSig)
    #rownames(deseq2_results())
  })
  
  output$downloadDEG1data <- downloadHandler(
    filename = "DESeq2_output.csv",
    
    content = function(file) {
      
      write.csv(DESeq2_output(),file = file,row.names = FALSE)
      
      #dev.off()
      #dev.copy2pdf(file)
    }
  )
  
  ###############--------------------DEG2----------------------------------------------------------------------------
  output$heatmap1 <- renderPlot({
    
    df1 <- as.data.frame(deseq2_results())
    colnames(df1)[1] <- "Genes"
    deq <<- df1
    
    resFC <- subset(df1 ,log2FoldChange < input$num)
    resSig <- subset(resFC, padj < input$num1)
    resOrdered <- resSig[order(resSig$padj),]
    
    #resOrdered <- deseq2_results()[order(deseq2_results()$padj),]
    topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 0, ], 
                         resOrdered[ resOrdered[,'log2FoldChange'] < 0, ] )
    
    
    topResults <- topResults[!duplicated(topResults[,1]), ]
    row.names(topResults) <- topResults[,1]
    topResults[,1] <- NULL
    
    hmcol <- brewer.pal(11,'RdBu')
    nCounts <- counts(dds(), normalized=TRUE)
    
    X <- cbind.data.frame(dis_data()[,2],nCounts)
    colnames(X)[1] <- "Genes"
    X <- X[!duplicated(X[,1]), ]
    row.names(X) <- X[,1]
    X[,1] <- NULL
    
    
    tem1 <- as.matrix(X[row.names(topResults), ])
    
    
    fir <- grep(input$select,colnames(tem1))
    up <<- tem1[,fir]
    
    sec <- grep(input$select1,colnames(tem1))
    down <<- tem1[,sec]
    
    sect <- tem1[,c(fir,sec)]
    AA <<- sect
    stats::heatmap(sect, Rowv = NA, col = hmcol,Colv = NA, mar = c(8,2))
    #dev.off()
  })
  
  output$selectGO5 <- renderUI({
    resFC <- subset(deseq2_results() ,log2FoldChange < input$num)
    resSig <- subset(resFC, padj < input$num1)
    
    Upregulated <<- rownames(subset(resSig,log2FoldChange > 0))#query = rownames(kmeans_df()$x)[which(Kmeans()$bar == 1)]
    query = Upregulated
    converted <- convertID(query,dl())
    gInfo <- geneInfo(converted)
    categoryChoices <- gmtCategory(converted)
    
    selectInput("selectGO5", label = h3("Select the Pathway database"), 
                choices = categoryChoices[-1]) ## because we dont want all selected gene from caterogory choices  
  })
  
  differential_pathway <- reactive({
    resFC <- subset(deseq2_results() ,log2FoldChange < input$num)
    resSig <- subset(resFC, padj < input$num1)
    
    Upregulated <<- rownames(subset(resSig,log2FoldChange > 0))
    Doregulated <<- rownames(subset(resSig,log2FoldChange < 0))
    
    query = Uregulated
    converted <- convertID(query,dl())
    gInfo <- geneInfo(converted)
    categoryChoices <- gmtCategory(converted)
    result <- FindOverlap(converted,gInfo,input$selectGO5)
    result$direction = "Up Regulated"
    
    query = Doregulated
    converted <- convertID(query,dl())
    gInfo <- geneInfo(converted)
    categoryChoices <- gmtCategory(converted)
    result1 <- FindOverlap(converted,gInfo,input$selectGO5)
    result1$direction = "Down Regulated"
    
    results <- rbind(result,result1)
    #x <- results
    n = nrow(results)
    tem=rep(TRUE,n)
    
    geneLists = lapply(results$Genes, function(y) unlist( strsplit(as.character(y)," " )   ) )
    for( i in 2:n)
      for( j in 1:(i-1) ) { 
        if(tem[[j]]) { # skip if this one is already removed
          commonGenes = length(intersect(geneLists[[i]] ,geneLists[[j]] ) )
          if( commonGenes/ length(geneLists[[j]] ) > input$redundant1 )
            tem[[i]] = FALSE	
        }			
      }								
    results <- results[which(tem),]
    
    minFDR  <- 0.05
    results= results[,c(6,1,2,4)]
    colnames(results)= c("Cluster","FDR","nGenes","Pathways")
    if(min(results$FDR) > minFDR ) results = as.data.frame("No signficant enrichment found.") else
      results= results[which(results$FDR < minFDR),]
    colnames(results)[2] = "adj.Pval"
    return(results)
  })
  
  output$degenrichment <- renderTable({
    differential_pathway()
  })
  
  differential_pathway_dw <- reactive({
    resFC <- subset(deseq2_results() ,log2FoldChange < input$num)
    resSig <- subset(resFC, padj < input$num1)
    
    Upregulated <<- rownames(subset(resSig,log2FoldChange > 0))
    Doregulated <<- rownames(subset(resSig,log2FoldChange < 0))
    
    query = Uregulated
    converted <- convertID(query,dl())
    gInfo <- geneInfo(converted)
    categoryChoices <- gmtCategory(converted)
    result <- FindOverlap(converted,gInfo,input$selectGO5)
    result$direction = "Up Regulated"
    
    query = Doregulated
    converted <- convertID(query,dl())
    gInfo <- geneInfo(converted)
    categoryChoices <- gmtCategory(converted)
    result1 <- FindOverlap(converted,gInfo,input$selectGO5)
    result1$direction = "Down Regulated"
    
    results <- rbind(result,result1)
    #x <- results
    n = nrow(results)
    tem=rep(TRUE,n)
    
    geneLists = lapply(results$Genes, function(y) unlist( strsplit(as.character(y)," " )   ) )
    for( i in 2:n)
      for( j in 1:(i-1) ) { 
        if(tem[[j]]) { # skip if this one is already removed
          commonGenes = length(intersect(geneLists[[i]] ,geneLists[[j]] ) )
          if( commonGenes/ length(geneLists[[j]] ) > input$redundant1 )
            tem[[i]] = FALSE	
        }			
      }								
    results <- results[which(tem),]
    
    minFDR  <- 0.05
    results= results[,c(6,1,2,4,5)]
    colnames(results)= c("Cluster","FDR","nGenes","Pathways","Genes")
    if(min(results$FDR) > minFDR ) results = as.data.frame("No signficant enrichment found.") else
      results= results[which(results$FDR < minFDR),]
    colnames(results)[2] = "adj.Pval"
    return(results)
  })
  output$downloadDEG2data <- downloadHandler(
    filename = "Enriched Differential pathway.csv",
    
    content = function(file) {
      
      write.csv(differential_pathway_dw(),file = file,row.names = FALSE)
      
      #dev.off()
      #dev.copy2pdf(file)
    }
  )
  
  #######################-----------------------------WGCNA----------------------------------------
  
  wgcna <- reactive({
    df_log <- cbind.data.frame(dis_data()[,2],df_transformed())
    colnames(df_log)[1] <- colnames(dis_data()[2])
    df_log <- df_log[!duplicated(df_log[,1]), ]
    row.names(df_log) <- df_log[,1]
    df_log[,1] <- NULL 
    aaaaa <<- df_log
    
    
    
    
    
    var_genes <- apply(df_log, 1, var)
    select_var <- names(sort(var_genes, decreasing=TRUE))[1:input$GENES]
    suuu <<- select_var
    
    highly_variable_lcpm <- df_log[select_var,]
    supp <<- highly_variable_lcpm
    
    WGCNA_matrix <- t(highly_variable_lcpm)
    
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    
    
    sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 5,corFnc = cor,corOptions = list(use = 'p'),networkType = "unsigned")
    f <<- sft
    
    adj= adjacency(WGCNA_matrix,type = "unsigned", power = input$soft);
    
    
    #turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
    TOM=TOMsimilarityFromExpr(WGCNA_matrix,networkType = "unsigned", TOMType = "unsigned", power = input$soft);
    subGeneNames <- colnames(WGCNA_matrix)
    colnames(TOM) =subGeneNames
    rownames(TOM) =subGeneNames
    
    
    
    geneTree1 = flashClust(as.dist(1-TOM), method = 'average')
    
    #module identification using dynamic tree cut algorithm
    modules = cutreeDynamic(dendro = geneTree1,method = "tree",
                            minClusterSize = input$Mod)
    
    modules <<- modules
    #assign module colours
    module.colours = labels2colors(modules)
    
    moduleInfo = cbind( subGeneNames, module.colours, modules)
    moduleInfo = moduleInfo[which(moduleInfo[,2] != "grey") ,] # remove genes not in any modules
    moduleInfo = moduleInfo[order(moduleInfo[,3]),] # sort
    
    n.modules = length(unique(module.colours) ) -1 ; nGenes = dim(moduleInfo)[1]	
    
    yess <- list(x = highly_variable_lcpm,powers=powers,sft=sft, TOM = TOM, dynamicColors = module.colours, moduleInfo = moduleInfo,n.modules=n.modules, nGenes =nGenes) 
    yess
    
    yo <<- yess
  })
  
  output$dendogram <- renderPlot({
    diss1 = 1-wgcna()$TOM;
    dynamicColors = wgcna()$dynamicColors
    
    hier1=flashClust(as.dist(diss1), method="average" )
    #set the diagonal of the dissimilarity to NA 
    diag(diss1) = NA;
    plotDendroAndColors(hier1, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
  })
  
  moduleinfo <- reactive({
    df_log <- cbind.data.frame(dis_data()[,2],df_transformed())
    colnames(df_log)[1] <- colnames(dis_data()[2])
    df_log <- df_log[!duplicated(df_log[,1]), ]
    row.names(df_log) <- df_log[,1]
    df_log[,1] <- NULL
    
    var_genes <- apply(df_log, 1, var)
    select_var <- names(sort(var_genes, decreasing=TRUE))[1:input$GENES]
    bbbb <<-  select_var
    
    
    highly_variable_lcpm <- df_log[select_var,]
    
    WGCNA_matrix <- t(highly_variable_lcpm)
    
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    
    sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 5,corFnc = cor,corOptions = list(use = 'p'),networkType = "unsigned")
    
    adj= adjacency(WGCNA_matrix,type = "unsigned", power = input$soft);
    
    
    #turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
    TOM=TOMsimilarityFromExpr(WGCNA_matrix,networkType = "unsigned", TOMType = "unsigned", power = input$soft);
    subGeneNames <- colnames(WGCNA_matrix)
    colnames(TOM) =subGeneNames
    rownames(TOM) =subGeneNames
    
    geneTree1 = flashClust(as.dist(1-TOM), method = 'average')
    
    #module identification using dynamic tree cut algorithm
    modules = cutreeDynamic(dendro = geneTree1,method = "tree",
                            minClusterSize = input$Mod)
    
    #assign module colours
    module.colours = labels2colors(modules)
    
    moduleInfo = cbind( subGeneNames, module.colours, modules)
    moduleInfo = moduleInfo[which(moduleInfo[,2] != "grey") ,] # remove genes not in any modules
    moduleInfo = moduleInfo[order(moduleInfo[,3]),] # sort
    #return(moduleInfo)
    
    moduleInfo
  })
  
  output$Modules <- renderUI({
    modules = unique(moduleinfo()[, c("modules","module.colours")] )
    moduleList = apply(modules,1,paste,collapse=". ")
    moduleList  = paste0( moduleList, " (", table(moduleinfo()[,"modules"] )," genes)"  )
    moduleList = c(moduleList,"Entire network")
    
    selectInput("selectmodule", h3("Select Module"), 
                choices = moduleList)
  })
  
  output$network <- renderPlot({
    modules = unique(moduleinfo()[, c("modules","module.colours")] )
    moduleList = apply(modules,1,paste,collapse=". ")
    moduleList  = paste0( moduleList, " (", table(moduleinfo()[,"modules"] )," genes)"  )
    moduleList = c(moduleList,"Entire network")
    
    
    module = unlist(strsplit(input$selectmodule ," " ) )[2]
    
    moduleColors = moduleinfo()[,"module.colours"]
    inModule = (moduleColors==module);
    
    datExpr = t(wgcna()$x )
    probes = colnames(datExpr)
    modProbes = probes[inModule];
    
    modTOM = wgcna()$TOM[inModule, inModule];
    dimnames(modTOM) = list(modProbes, modProbes)
    
    nTop = input$topgenes;
    
    if( nTop > 1000) nTop = 1000; 
    
    IMConn = softConnectivity(datExpr[, modProbes]);
    top = (rank(-IMConn) <= nTop)
    
    query <- probes
    
    
    
    converted <- convertID(query,dl())
    gInfo <- geneInfo(converted)
    
    probeToGene = gInfo[,c("ensembl_gene_id","symbol")]
    probeToGene$symbol = gsub(" ","",probeToGene$symbol)
    
    ix = which( is.na(probeToGene$symbol) |
                  nchar(probeToGene$symbol)<2 | 
                  toupper(probeToGene$symbol)=="NA" |  
                  toupper(probeToGene$symbol)=="0"  ) 			
    probeToGene[ix,2] = probeToGene[ix,1]  # use gene ID
    
    net <- modTOM[top,top] > input$threshold
    
    ix = match( colnames(net), probeToGene[,2])		
    colnames(net) = probeToGene[ix,2]
    ix = match( rownames(net), probeToGene[,2])		
    rownames(net) = probeToGene[ix,2]	
    
    library(igraph,verbose=FALSE)
    #plot(graph_from_data_frame(d=data.frame(1:10,ncol=2)  ,directed=F) )
    # http://www.kateto.net/wp-content/uploads/2016/01/NetSciX_2016_Workshop.pdf
    plot( graph_from_adjacency_matrix( net, mod ="undirected" ), 
          vertex.label.color="black", vertex.label.dist=3,vertex.size=7)
    
    
    
  })
  
  Enriched_pathway <- reactive ({
    
    modules = unique(moduleinfo()[, c("modules","module.colours")] )
    moduleList = apply(modules,1,paste,collapse=". ")
    moduleList  = paste0( moduleList, " (", table(moduleinfo()[,"modules"] )," genes)"  )
    moduleList = c(moduleList,"Entire network")
    
    
    module = unlist(strsplit(input$selectmodule ," " ) )[2]
    
    moduleColors = moduleinfo()[,"module.colours"]
    inModule = (moduleColors==module);
    
    datExpr = t(wgcna()$x)
    probes = colnames(datExpr)
    modProbes = probes[inModule];
    
    query <- modProbes
    converted <- convertID(query,dl())
    gInfo <- geneInfo(converted)
    categoryChoices <- gmtCategory(converted)
    result <- FindOverlap(converted,gInfo,input$selectGO4)
    #result$direction = toupper(letters) 
    results <- result
    #results <- rbind(results,result)
    
    
    minFDR  <- 0.05
    results= results[,c(1,2,4)]
    colnames(results)= c("FDR","nGenes","Pathways")
    if(min(results$FDR) > minFDR ) results = as.data.frame("No signficant enrichment found.") else
      results= results[which(results$FDR < minFDR),]
    colnames(results)[1] = "adj.Pval"
    return(results)
  })
  
  output$selectGO4 <- renderUI({
    modules = unique(moduleinfo()[, c("modules","module.colours")] )
    moduleList = apply(modules,1,paste,collapse=". ")
    moduleList  = paste0( moduleList, " (", table(moduleinfo()[,"modules"] )," genes)"  )
    moduleList = c(moduleList,"Entire network")
    
    
    module = unlist(strsplit(input$selectmodule ," " ) )[2]
    
    moduleColors = moduleinfo()[,"module.colours"]
    inModule = (moduleColors==module);
    
    datExpr = t(wgcna()$x)
    probes = colnames(datExpr)
    modProbes = probes[inModule];
    
    query <- modProbes
    converted <- convertID(query,dl())
    gInfo <- geneInfo(converted)
    categoryChoices <- gmtCategory(converted)
    
    selectInput("selectGO4",h3("Select the Pathway"), 
                choices = categoryChoices[-1]) ## because we dont want all selected gene from caterogory choices  
  })
  
  output$enriched_pathways <- renderTable({
    as.data.frame(Enriched_pathway())
  })
  
  Enriched_pathway_dw <- reactive({
    modules = unique(moduleinfo()[, c("modules","module.colours")] )
    moduleList = apply(modules,1,paste,collapse=". ")
    moduleList  = paste0( moduleList, " (", table(moduleinfo()[,"modules"] )," genes)"  )
    moduleList = c(moduleList,"Entire network")
    
    
    module = unlist(strsplit(input$selectmodule ," " ) )[2]
    
    moduleColors = moduleinfo()[,"module.colours"]
    inModule = (moduleColors==module);
    
    datExpr = t(wgcna()$x)
    probes = colnames(datExpr)
    modProbes = probes[inModule];
    
    query <- modProbes
    converted <- convertID(query,dl())
    gInfo <- geneInfo(converted)
    categoryChoices <- gmtCategory(converted)
    result <- FindOverlap(converted,gInfo,input$selectGO4)
    #result$direction = toupper(letters) 
    results <- result
    #results <- rbind(results,result))
    minFDR  <- 0.05
    results= results[,c(6,1,2,4,5)]
    colnames(results)= c("Cluster","FDR","nGenes","Pathways","Genes")
    if(min(results$FDR) > minFDR ) results = as.data.frame("No signficant enrichment found.") else
      results= results[which(results$FDR < minFDR),]
    colnames(results)[2] = "adj.Pval"
    return(results)
  })
  
  output$downloadWGCNAdata <- downloadHandler(
    filename = "Enriched WGCNA pathway.csv",
    
    content = function(file) {
      
      write.csv(Enriched_pathway_dw(),file = file,row.names = FALSE)
      
      #dev.off()
      #dev.copy2pdf(file)
    }
  )
})


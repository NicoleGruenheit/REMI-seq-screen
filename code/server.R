suppressMessages(library(tidyverse))
suppressMessages(library(GGally))
suppressMessages(library(data.table))
suppressMessages(library(gtools))
source('normalise_reads.R', encoding = 'UTF-8')

# Define server logic 
server <- function(input, output,session) {
  # step 1
  samples = reactive({
    inFile <- input$file1
    f = read.table(inFile$datapath,header=T)
    f = as.data.frame(f)
    updateSelectInput(session,"sample",label="",choices=unique(f[,1]))
    updateSelectInput(session,"sample5",label="",choices=unique(f[,1]))
    updateSelectInput(session,"samplePlot",label="",choices=unique(f[,1]))
    f
  })
  
  output_expdesign_data <- reactive({
    samples()
  })
  
  output$exp_design <- DT::renderDataTable({
    inFile <- input$file1
    if(is.null(inFile)){
      f = data.frame("Sample","Replicates","Filename","Name")
      DT::datatable(f, extensions = c('ColReorder'), rownames = FALSE, options = list(scrollX = T, dom = 'rltip', colReorder = TRUE,fixedHeader = TRUE),escape = F)
    }else{
      DT::datatable(output_expdesign_data(), extensions = c('ColReorder'), rownames = FALSE, options = list(scrollX = T, dom = 'rltip', colReorder = TRUE,fixedHeader = TRUE),escape = F)
    }  
  })
  
  # step 2
  
  stats <- reactiveValues(df_data = NULL)
  pos <- reactiveValues(df_data = NULL)
  s <- reactiveValues(df_data = NULL)
  directory <- reactiveValues(df_data = NULL)
  
  observeEvent(ignoreNULL = TRUE,
  	eventExpr = {
  		input$directory
  	},
  	handlerExpr = {
  		if (input$directory > 0) {
  			# condition prevents handler execution on initial app launch
  			directory = choose.dir(default = readDirectoryInput(session, 'directory'),
  												caption="Choose a directory...")
  			updateDirectoryInput(session, 'directory', value = directory)
  		}
  	}
  )
  
  output$directory = renderText({
  	readDirectoryInput(session, 'directory')
  })
  
  
  observeEvent(input$get_tags, {
    directory = readDirectoryInput(session, 'directory')
    annots = input$file2
    indices = input$file3
    cmd <- paste("perl", "get_tags.pl", "-files", directory, "-annot",annots,"-indices",indices,"-cutoff 0")
    system(cmd)
    stats$df_data = read.table(paste0(dirname(directory),"/",basename(directory),"_stats"),header=T)
    s$df_data = length(levels(stats$df_data$filename))
    pos$df_data = read.table(paste0(dirname(directory),"/",basename(directory),"_positions_separate"),header=F,sep="\t",stringsAsFactors = FALSE)
  })
  
  observeEvent(input$load_results, {
    directory = readDirectoryInput(session, 'directory')
    stats$df_data = read.table(paste0(dirname(directory),"/",basename(directory),"_stats"),header=T)
    #print(length(levels(stats$df_data$filename)))
    s$df_data = length(levels(stats$df_data$filename))
    pos$df_data = read.table(paste0(dirname(directory),"/",basename(directory),"_positions_separate"),header=F,sep="\t",stringsAsFactors = FALSE)
  })
  
  output$stats <- DT::renderDataTable({
    
    if(is.null(stats)){
      f = data.frame("filename","reads","reads_with_tags","percentage_reads_with_tags","unique_tags","tags_in_inverted_repeat","not_unique_tags")
      DT::datatable(f, extensions = c('ColReorder'), rownames = FALSE, options = list(scrollX = T, dom = 'rltip', colReorder = TRUE,fixedHeader = TRUE),escape = F)
    }else{
      directory = input$directory
      DT::datatable(stats$df_data, extensions = c('ColReorder'), rownames = FALSE, options = list(scrollX = T, dom = 'rltip', colReorder = TRUE,fixedHeader = TRUE),escape = F)
    }  
  })
  
  # step 3
  norm <- reactiveValues(df_data = NULL)
  summary <- reactiveValues(df_data = NULL) 
  
  positions = reactive({
    comb =  input$radio
    directory = input$directory
    cutoff = input$cutoff
    norm$df_data = normalise_reads(pos$df_data,comb,cutoff)
    test = as.data.table(norm$df_data)
    if(!(comb == 3 | comb == 4)){
      print(comb)
      colnames(test)[colnames(test)=="V6"] <- "Chromosome"
      colnames(test)[colnames(test)=="V7"] <- "Position"
    }
    as.data.frame(test)
  })
  
  output$norm <- DT::renderDataTable({
    if(is.null(positions)){
    
    }else{
      DT::datatable(positions(), extensions = c('ColReorder'), rownames = FALSE, options = list(scrollX = T, dom = 'rltip', colReorder = TRUE,fixedHeader = TRUE),escape = F)
    }  
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$directory, ".csv", sep = "")
    },
    content = function(file) {
      write.table(positions(), file, row.names = FALSE,sep="\t")
    }
  )
  
  output$summary <- DT::renderDataTable({
    if(is.null(positions)){
    }else{
      DT::datatable(summary_stats(), extensions = c('ColReorder'), rownames = FALSE, options = list(scrollX = T, dom = 'rltip', colReorder = TRUE,fixedHeader = TRUE),escape = F)
    }  
  })
  
  #step 4
  Plot1 <- reactive({
    names = filter(samples(),Sample == input$sample)$Name
    number = length(unique(names))
    df.subset <- positions()[, names]
    if(input$g0){
      collist = names(df.subset)
      sel <- apply(df.subset[,collist],1,function(row) !1 %in% row)
      df.subset = df.subset[sel,]
    }
    if(input$g1){
      
    }
    ggpairs(log(df.subset,base=10), aes(alpha = 0.4))
  })
  
  output$plot1 <- renderPlot({
    validate(
      need(input$file1 != "", "Please load the experimental design on tab 1")
    )
    p=Plot1()
    print(p)
  })
  
  output$downloadPlots1 <- downloadHandler(
    filename = function() { paste(input$directory, '_correlation_plots.pdf', sep='') },
    content = function(file) {
      pdf(file)
      print(Plot1())
      dev.off()
    }
  )
  
  #step 5
  folds = reactive({
    validate(
      need(input$file1 != "", "Please load the experimental design on tab 1")
    )
    start = input$sample5
    names_start = filter(samples(),Sample == input$sample5)$Name
    names_other = filter(samples(),!(Sample == input$sample5))$Name
    data = positions() %>% select(contains(start)) %>% mutate(ref = rowSums(.)) %>% select(-contains(start))
    data1 = positions() %>% select(-contains(start))
    data1 = cbind(data1,data)
    data1 = gather(data1,sample,value,one_of(names_other))
    data1 = mutate(data1,FC = foldchange(value,ref),lFC = foldchange2logratio(FC,base=10))
    data1 = mutate(data1,bin = ifelse(ref < 100,100,ifelse(ref<1000,1000,10000)))
    data2 = spread(select(data1,Chromosome,Position,sample,value),sample,value)
    data1 = select(data1,-value,-FC)
    data1 = spread(data1,sample,lFC)
    data1 = left_join(data1,data2,by=c("Chromosome","Position"))
    data1
  })
  
  output$fold <- DT::renderDataTable({
    if(is.null(positions)){
    }else{
      DT::datatable(folds(), extensions = c('ColReorder'), rownames = FALSE, options = list(scrollX = T, dom = 'rltip', colReorder = TRUE,fixedHeader = TRUE),escape = F)
    }  
  })
  
  output$downloadData_FC <- downloadHandler(
    filename = function() {
      paste(input$directory, ".csv", sep = "")
    },
    content = function(file) {
      write.table(folds(), file, row.names = FALSE,sep="\t")
    }
  )
  
  Plot2 <- reactive({
    choice1 = input$samplePlot
    print(choice1) 
    validate(
      need(input$samplePlot != input$sample5, "Please choose a different sample")
    )
    help = folds() %>% select(contains("Chromosome","Position",choice1)) %>% select("Chromosome","Position",contains(".y"))
    help = help %>% gather(sample,value,3:dim(help)[2])
    help = filter(help,value==1)
    data = folds() %>% select(contains(choice1)) %>% select(contains(".x"))
    reps = dim(data)[2]
    data1 <- data.frame(x=data[,1], y=data[,2], cor=paste("1","2",sep="_"))
    if(reps > 2){
      for(i in 1:reps) {
        for(j in 1:reps){
          if(!i == j & !(i == 1 & j == 2)){
            data1 <- rbind(data1,data.frame(x=data[,i], y=data[,j], cor=paste(i,j,sep="_")))
          }
        }
      }
    }
    print(dim(data1))
    ggplot(data1, aes(x, y))+geom_point(size=0.5,alpha=0.5)  + geom_hline(yintercept=0, color="red",linetype="dashed")+ geom_vline(xintercept=0, color="red",linetype="dashed") + facet_wrap(~cor, scales="free")
  })
  
  output$plot2 <- renderPlot({
    validate(
      need(input$file1 != "", "Please load the experimental design on tab 1")
    )
    p=Plot2()
    print(p)
  })
  
  summary_stats = reactive({
    data = positions()
    data1 = data.frame()
    all = dim(data)[1]
    print(head(data))
    print(head(folds()))
    data2 = gather(data,sample,count,(dim(data)[2]-s$df_data + 1):dim(data)[2])
    data2 = filter(data2,count > 1) %>% group_by(sample) %>% summarise(g1 = n())
    data2 = spread(data2,sample,g1)
    print(head(data2))
    data1 = cbind(all,data2)
    data1
  })
  
  output$sumstats <- DT::renderDataTable({
    if(is.null(positions)){
    }else{
      DT::datatable(summary_stats(), extensions = c('ColReorder'), rownames = FALSE, options = list(scrollX = T, dom = 'rltip', colReorder = TRUE,fixedHeader = TRUE),escape = F)
    }  
  })
}
#   Setting     #################

library(knitr)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(writexl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(fgsea)
library(enrichplot)
library(DT)
library(stringr)
library(shiny)
library(bslib)
library(readxl)
library(RColorBrewer)
library(purrr)
library(ggvenn)
library(forcats)

MyPalette <- c("#9933aa","#ffdd22","#aa4400","#ff0000","#337722","#00ff66","#005566","#002277",
               "#441144","#aa0077","#00bbff","#003333","#4422cc","#116611","#330077","#111111",
               "#667700","#ddaa00","#33ffff","#ff22ff","#ffff33","#00ff00","#0000ff","#444444")


#   Functions     ####################

Display_Venn <- function(Markers, colpalette = NULL, set.names = NULL) {
  
  # package install checking
  pkgs <- c("ggvenn","purrr")
  for(pkg in pkgs){
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
  
  library(ggvenn)
  library(purrr)
  
  # Arguments validation
  if (is.null(Markers)) {
    stop("Argument 'Markers' is NULL. Please provide a list of sets.")
  }
  
  if (!is.list(Markers)) {
    stop("Argument 'Markers' must be a list of vectors.")
  }
  
  if (!all(sapply(Markers, is.vector))) {
    stop("All elements in the 'Markers' list must be vectors.")
  }
  
  
  n <- ifelse(length(Markers) >= 4, 4, length(Markers))
  
  if (n < 2) {
    stop("You need at least 2 sets to make a comparison.")
  }
  
  if (length(Markers) > 4) {
    warning("The function supports up to 4 sets only. Taking the first 4 vectors as arguments")
  }
  
  # Handling colours
  default_palette <- c("navy", "red", "darkgreen", "violet")
  if (is.null(colpalette) || length(colpalette) < n) {
    warning("The provided colour palette is too short or NULL. Using default colours.")
    colpalette <- default_palette[1:n]
  } else if(length(colpalette) >= n){
    colpalette <- colpalette[1:n]
  }
  
  # Handling sets
  sets <- Markers[1:n]
  if (!(is.null(set.names))) {
    if (length(set.names) != n) {
      stop("The length of 'set.names' must match the number of sets in 'Markers'")
    }
    names(sets) <- set.names
  } else if (is.null(names(sets))) {
    names(sets) <- paste0("Set", 1:n)
  }
  
  # Calculate intersections
  common_genes <- list()
  for(i in 2:n){
    x <- combn(sets, i, simplify = F)
    for(s in 1:length(x)){
      if(i == 2){
        common_genes[[paste(names(x[[s]]), collapse = " ∩ ")]] <- intersect(x[[s]][[1]], x[[s]][[2]])
      }
      
      if(i == 3){
        common_genes[[paste(names(x[[s]]), collapse = " ∩ ")]] <- purrr::reduce(list(x[[s]][[1]],
                                                                                     x[[s]][[2]],
                                                                                     x[[s]][[3]]), intersect)
      }
      
      if(i == 4){
        common_genes[["Common all"]] <- purrr::reduce(x[[s]], intersect)
      }
    }
  }
  
  # Calculate group specific genes
  specific_genes <- list()
  for(i in 2:n){
    x <- combn(sets, i, simplify = F)
    for(s in 1:length(x)){
      if(i == 2){
        specific_genes[[paste(names(x[[s]]), collapse = " or ")]] <- setNames(list(setdiff(x[[s]][[1]], 
                                                                                           intersect(x[[s]][[1]], x[[s]][[2]])),
                                                                                   setdiff(x[[s]][[2]], 
                                                                                           intersect(x[[s]][[1]], x[[s]][[2]]))),
                                                                              c(names(x[[s]][1]), names(x[[s]][2])))
      }
      
      if(i == 3){
        specific_genes[[paste(names(x[[s]]), collapse = " or ")]] <- setNames(list(setdiff(x[[s]][[1]], 
                                                                                           union(x[[s]][[2]], x[[s]][[3]])),
                                                                                   setdiff(x[[s]][[2]], 
                                                                                           union(x[[s]][[1]], x[[s]][[3]])),
                                                                                   setdiff(x[[s]][[3]], 
                                                                                           union(x[[s]][[2]], x[[s]][[1]]))),
                                                                              c(names(x[[s]][1]), names(x[[s]][2]), names(x[[s]][3])))
      }
      
      if(i == 4){
        specific_genes[["unique each"]] <- setNames(list(setdiff(x[[s]][[1]],
                                                                 purrr::reduce(list(x[[s]][[2]], x[[s]][[3]], x[[s]][[4]]), 
                                                                               union)),
                                                         setdiff(x[[s]][[2]],
                                                                 purrr::reduce(list(x[[s]][[1]], x[[s]][[3]], x[[s]][[4]]), 
                                                                               union)),
                                                         setdiff(x[[s]][[3]],
                                                                 purrr::reduce(list(x[[s]][[2]], x[[s]][[1]], x[[s]][[4]]), 
                                                                               union)),
                                                         setdiff(x[[s]][[4]],
                                                                 purrr::reduce(list(x[[s]][[2]], x[[s]][[3]], x[[s]][[1]]), 
                                                                               union))),
                                                    c(names(x[[s]][1]),names(x[[s]][2]),names(x[[s]][3]),names(x[[s]][4])))
      }
    }
  }
  
  # Plot
  p <- ggvenn(
    sets,
    fill_color = colpalette,
    stroke_size = .8,
    show_elements = F,
    stroke_linetype = "solid",
    set_name_color = colpalette,
    set_name_size = 5,
    text_color = "black",
    text_size = 5,
    padding = 0.03, 
    show_stats = "c", 
    fill_alpha = 0.4
  )
  
  return(list(plot = p, intersections = common_genes, group_specific = specific_genes))
}


#   Shiny Application   ##################

options(shiny.maxRequestSize = 900*1024^2)




#  User interface    ####################

ui <- fluidPage(
  theme = bs_theme(bootswatch = "flatly"),
  h2("1- Importing the data", 
     style = "color:gold ; font-weight:700 ; background-color:black ; text-align:left ; margin-top:20px ; margin-bottom:30px"),
  
  ##  I-  Import the data     ---------------------------------------
  
  fluidRow(
    column(
      width = 4,
      fileInput("mk1", "List markers Lib1", accept = ".rds"),
      verbatimTextOutput("str1")
    ),
    column(
      width = 4,
      fileInput("mk2", "List markers Lib2", accept = ".rds"),
      verbatimTextOutput("str2")
    ),
    column(
      width = 4,
      fileInput("mk3", "List markers Lib3", accept = ".rds"),
      verbatimTextOutput("str3")
    )
  ),
  
  h2("2- Setting the contrasts", 
     style = "color:gold ; font-weight:700 ; background-color:black ; text-align:left ; margin-top:50px ; margin-bottom:20px"),
  
  ##  II- Set the comparisons     --------------------------------------------------
  
  fluidRow(
    column(
      width = 4,
      selectInput("sel1","Select from Lib1", choices = c(), multiple = T, selectize = F, width = "400px", size = 10),
      p("", style = "margin-top:20px ; margin-bottom:30px"),
      radioButtons("radio1","Select genes", choices = c("Up","Down","All"), selected = ""),
      p("", style = "margin-top:20px ; margin-bottom:30px"),
      actionButton("act1", "Set up selection", width = "150px", class = "btn-sm btn-primary")
    ),
    column(
      width = 4,
      selectInput("sel2","Select from Lib2", choices = c(), multiple = T, selectize = F, width = "400px", size = 10),
      p("Selected genes", style = "margin-top:35px ; color:purple ; font-weight:600"),
      verbatimTextOutput("selected")
    ),
    column(
      width = 4,
      selectInput("sel3","Select from Lib3", choices = c(), multiple = T, selectize = F, width = "400px", size = 10)
    )
  ),
  
  ##  III-  Comparison      ----------------------------------------------------------
  
  h2("3- Comparison", 
     style = "color:gold ; font-weight:700 ; background-color:black ; text-align:left ; margin-top:50px ; margin-bottom:20px"),
  
  fluidRow(
    column(
      width = 4, offset = 4,
      actionButton("act2","Display ggvenn", width = "400px", class = "btn-lg btn-success"),
      p("", style = "margin-top:30px ; margin-bottom:20px")
    )
  ),
  
  fluidRow(
    column(
      width = 8, offset = 2,
      plotOutput("plt1"),
      p("", style = "margin-top:30px ; margin-bottom:20px")
    )
  ),
  
  ##  IV- Displaying genes      -------------------------------------------------
  
  fluidRow(
    column(
      width = 6,
      p("See the list of the common genes",
        style = "color:purple ; font-weight:600 ; font-size:120% ; margin-top:20px"),
      selectInput("sel.set","Select your intersection", choices = c(), width = "600px"),
      verbatimTextOutput("nb_common_genes"),
      verbatimTextOutput("common_genes"),
    ),
    column(
      width = 6,
      p("Genes specific to a condition", 
        style = "color:purple ; font-weight:600 ; font-size:120% ; margin-top:20px"),
      selectInput("specific.genes", "Select the contrast", choices = c(), width = "600px"),
      selectInput("specific.genes2", "Select the condition", choices = c(), width = "600px"),
      verbatimTextOutput("nb_genes"),
      verbatimTextOutput("specific_genes")
    )
  )
  
)











#   Server function     ########################

server <- function(input, output, session){
  
  #  I-  Import the data     -------------------------------
  
  Lib1 <- reactive({
    req(input$mk1)
    Lib1 <- readRDS(input$mk1$datapath)
  })
  Lib2 <- reactive({
    req(input$mk2)
    Lib2 <- readRDS(input$mk2$datapath)
  })
  Lib3 <- reactive({
    req(input$mk3)
    Lib3 <- readRDS(input$mk3$datapath)
  })
  
  output$str1 <- renderPrint({
    req(Lib1())
    if(length(names(Lib1())) != 0){print("Lib1 markers loaded")}
  })
  output$str2 <- renderPrint({
    req(Lib2())
    if(length(names(Lib2())) != 0){print("Lib2 markers loaded")}
  })
  output$str3 <- renderPrint({
    req(Lib3())
    if(length(names(Lib3())) != 0){print("Lib3 markers loaded")}
  })
  
  
  #  II- Setting the contrast    ---------------------------------------------
  
  observeEvent(Lib1(),{
    updateSelectInput(session, "sel1", choices = names(Lib1()))
  })
  observeEvent(Lib2(),{
    updateSelectInput(session, "sel2", choices = names(Lib2()))
  })
  observeEvent(Lib3(),{
    updateSelectInput(session, "sel3", choices = names(Lib3()))
  })
  
  
  sets <- reactiveVal(list())
  
  Sets <- eventReactive(input$act1, {
    req(Lib1(),Lib2(),Lib3(), input$radio1)
    
    AllLists <- c(Lib1(),Lib2(),Lib3())
    tmp <- c(input$sel1, input$sel2, input$sel3)
    
    tmp_list <- list()
    for(i in tmp){
      tmp_list[[i]] <- AllLists[[i]][[input$radio1]]
    }
    
    sets(tmp_list)
    sets()
  })
  
  output$selected <- renderPrint({
    req(Sets())
    print(str(Sets()))
  })
  
  
  # III- Display ggvenn     -----------------------
  
  vennObj <- eventReactive(input$act2, {
    req(sets())
    sets <- sets()
    
    Display_Venn(sets)
  })
  
  VennPlot <- eventReactive(vennObj(), {
    req(vennObj())
    vennObj <- vennObj()
    vennObj[["plot"]]
  })
  
  output$plt1 <- renderPlot({
    req(VennPlot())
    VennPlot()
  })
  
  
  # IV- Display genes
  
  observeEvent(vennObj(), {
    req(vennObj())
    vennObj <- vennObj()
    
    updateSelectInput(session, "sel.set",
                      choices = names(vennObj[["intersections"]]), selected = "")
    updateSelectInput(session, "specific.genes",
                      choices = names(vennObj[["group_specific"]]), selected = "")
  })
  
  observeEvent(input$specific.genes, {
    req(vennObj())
    vennObj <- vennObj()
    
    updateSelectInput(session, "specific.genes2", 
                      choices = names(vennObj[["group_specific"]][[input$specific.genes]]), 
                      selected = "")
  })
  
  output$nb_common_genes <- renderPrint({
    req(vennObj())
    vennObj <- vennObj()
    print(length(vennObj[["intersections"]][[input$sel.set]]))
  })
  
  output$common_genes <- renderPrint({
    req(vennObj())
    vennObj <- vennObj()
    print(vennObj[["intersections"]][[input$sel.set]])
  })
  
  output$nb_genes <- renderPrint({
    req(vennObj())
    vennObj <- vennObj()
    print(length(vennObj[["group_specific"]][[input$specific.genes]][[input$specific.genes2]]))
  })
  
  output$specific_genes <- renderPrint({
    req(vennObj())
    vennObj <- vennObj()
    print(vennObj[["group_specific"]][[input$specific.genes]][[input$specific.genes2]])
  })
  
  
}





#   RUN shiny   #####################

shinyApp(ui, server)


























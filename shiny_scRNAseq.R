#   Setting     #################
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(writexl)
# library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ComplexHeatmap)
library(DT)
library(stringr)
library(shiny)
library(shinythemes)
library(shinyBS)
library(bslib)
library(circlize)
library(purrr)
library(ggvenn)
library(forcats)

MyPalette <- c("#9933aa","#ffdd22","#aa4400","#ff0000","#337722","#00ff66","#005566","#002277",
               "#441144","#aa0077","#00bbff","#003333","#4422cc","#116611","#330077","#111111",
               "#667700","#ddaa00","#33ffff","#ff22ff","#ffff33","#00ff00","#0000ff","#444444")


# Functions : ------------------------------------------------------------------
Assign_colors <- function(object, pal){
  n = length(unique(object))
  
  if(is.numeric(object)){
    message("This function can't handle numeric argument !!!")
    stop("Please provide a character vector as an object")
  }
  
  if(length(pal) < n){
    warning("Not enough colours within the selected palette to match your variables !! ",
            "Selecting a default palette from base R")
    pal = sample(colours(distinct = F), size = n)
  }
  
  res = sample(pal, size = n) %>% setNames(unique(object))
  
  return(res)
}


extract_term_genes <- function(object, goID){
  
  tryCatch({
    if(nrow(object) == 0){
      stop("The object you provided is empty (no go terms could be found) !")
    }
    
    if(!(goID %in% object[["ID"]])){
      stop("The goID you provided doesn't exist within the results, please check again !")
    }
    
    term_genes <- str_split(object[["geneID"]][which(object[["ID"]] == goID)], "/", simplify = F)
    res <- term_genes[[1]]
    
    
    return(res)
  }, error = function(e){
    message("Error: ", e$message)
  })
}


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


my_umap_theme <- function(){
  theme_minimal()+
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(size = 20, face = "bold", colour = "darkred", hjust = 0.5, margin = margin(b= 0.2, unit = "cm")),
          plot.subtitle = element_text(size = 17, face = "bold", colour = "#771144", hjust = 0.5, margin = margin(b= 0.5, unit = "cm")),
          legend.title = element_text(size = 17, face = "bold", colour = "darkred", hjust = 0.5, margin = margin(b= 0.5, unit = "cm")),
          legend.box.margin = margin(l= 0.5, unit = "cm"),
          legend.text = element_text(size = 14, colour = "black", face = "bold"))
}




#   Shiny App options   ##################
options(shiny.maxRequestSize = 900*1024^2)








#  Shiny ui    ####################
ui <- fluidPage(
  # theme = bs_theme(bootswatch = "flatly"),
  theme = shinytheme("flatly"),
  h1("I- Importing the data",
     style = "color:gold ; font-weight:700 ; background-color:black ; margin-top:10px ; margin-bottom:20px"),
  
  ## I- Importing data   --------------------------------------------------------
  fluidRow(
    column(width = 6,
           fileInput(inputId = "LoadSeurat",
                     label = p("Load your seurat object",
                               style = "color:darkred ; font-weight:600 ; font-size:130% ; margin-bottom:10px"),
                     accept = ".rds",
                     width = "400px")),
    column(width = 5,
           verbatimTextOutput("InfoSeurat"))
  ),
  
  ## II- Metadata selection   ----------------------------------------------------
  h1("II- Metadata selection",
     style = "color:gold ; font-weight:700 ; background-color:black ; margin-top:200px ; margin-bottom:20px"),
  fluidRow(
    column(width = 7,
           selectInput(inputId = "Sel.metadata",
                       label = "Select you metadata",
                       choices = c()),
           p("Select the column from your metadata from which you want to apply 
             a differential expressed analysis on the next steps",
             style = "font-style:italic ; color:darkgray ; font-size:90% ; margin-bottom:30px"),
           plotOutput("plt1", height = "500px")),
    column(width = 3,
           offset = 1,
           actionButton("act1", "View the summary", class = "btn-primary", 
                        width = "100%", style = "font-size:24px ; color:gold ; font-weight:600 ; border-radius:15px"),
           p("", style = "margin-top:30px"),
           tableOutput("tab1"))
  ),
  
  ## III- DEA   ------------------------------------------------------------------
  h1("III- Differential expression analysis",
     style = "color:gold ; font-weight:700 ; background-color:black ; margin-top:200px ; margin-bottom:20px"),
  
  ### 1- Making the contrasts   --------------------------------------------------
  p("1- Select the contrasts",
    style = "color:darkred ; font-weight:600 ; font-size:160% ; background-color:gold ; margin-bottom:20px"),
  fluidRow(
    column(width = 4,
           p("From the selected metadata above, choose the cluster you want to test and against which cluster ! 
             (e.g Relapse vs Diagnosis)",
             style = "font-size:110% ; font-weight:600 ; color:darkgreen ; margin-bottom:20px"),
           selectInput("ident.1", "Selection 1", choices = c(), width = "300px"),
           selectInput("ident.2", "Selection 2", choices = c(), width = "300px")),
    column(width = 6,
           offset = 1,
           plotOutput("plt2", height = "500px"))
  ),
  
  ### 2- Setting the parameters   ------------------------------------------------
  p("2- Setting the parameters & Run the DEA",
    style = "color:darkred ; font-weight:600 ; font-size:160% ; background-color:gold ; margin-top: 50px ; margin-bottom:20px"),
  fluidRow(
    column(width = 3,
           p("Initial parameters", 
             style = "color:darkgreen ; font-weight:700 ; font-size:130% ; margin-bottom:10px ; 
             margin-top:10px ; border-bottom: 3px solid darkred; padding-bottom: 3px ;
             border-top: 3px solid darkred; padding-top: 3px"),
           div(
             tags$span(
               style = "font-weight:600 ; font-size:110%",
               "Log2FC min.threshold ",
               icon("info-circle", id= "info-numeric", style = "color:midnightblue ; cursor:pointer")
             ),
             numericInput("log2fc",NULL,value = 0.5, min = 0, max = 3, step = 0.1, width = "100px")
           ),
           bsTooltip(id = "info-numeric",
                     title = "Define the minimum Log2FC threshold to be applied before running the DEA (in absolute value).",
                     placement = "top", trigger = "hover", options = list(container = "body")),
           div(
             tags$span(
               style = "font-weight:600 ; font-size:110%",
               "Min.pct threshold ",
               icon("info-circle", id= "info-numeric2", style = "color:midnightblue ; cursor:pointer")
             ),
             numericInput("min.pct",NULL,value = 0.1, min = 0, max = 1, step = 0.05, width = "100px")
           ),
           bsTooltip(id = "info-numeric2",
                     title = "Define the minimum percentage of cells within any cluster that have to express the gene",
                     placement = "top", trigger = "hover", options = list(container = "body")),
           div(
             tags$span(
               style = "font-weight:600 ; font-size:110%",
               "Min.diff.pct",
               icon("info-circle", id= "info-numeric3", style = "color:midnightblue ; cursor:pointer")
             ),
             numericInput("min.diff.pct",NULL,value = 0.1, min = 0, max = 1, step = 0.025, width = "100px")
           ),
           bsTooltip(id = "info-numeric3",
                     title = "Define the minimum difference in the fraction of detection between the 2 groups for each gene",
                     placement = "top", trigger = "hover", options = list(container = "body")),
           
           
           p("Additional filtering", 
             style = "color:darkgreen ; font-weight:700 ; font-size:130% ; margin-bottom:10px ; 
             margin-top:30px ; border-bottom: 3px solid darkred; padding-bottom: 3px ;
             border-top: 3px solid darkred; padding-top: 3px"),
           div(
             tags$span(
               style = "font-weight:600 ; font-size:110%",
               "Filter Up log2FC ",
               icon("info-circle", id= "info-slider1", style = "color:midnightblue ; cursor:pointer")
             ),
             sliderInput("log2fc_up",NULL,value = c(0.5,10), min = 0, max = 10, step = 0.25, width = "250px")
           ),
           bsTooltip(id = "info-slider1",
                     title = "Filter the up regulated genes by the selected range of log2FC",
                     placement = "top", trigger = "hover", options = list(container = "body")),
           div(
             tags$span(
               style = "font-weight:600 ; font-size:110%",
               "Filter Down log2FC ",
               icon("info-circle", id= "info-slider2", style = "color:midnightblue ; cursor:pointer")
             ),
             sliderInput("log2fc_down",NULL,value = c(-10,-0.5), min = -10, max = 0, step = 0.25, width = "250px")
           ),
           bsTooltip(id = "info-slider2",
                     title = "Filter the down regulated genes by the selected range of log2FC",
                     placement = "top", trigger = "hover", options = list(container = "body")),
           
           div(
             style = "margin-top:30px ; margin-bottom:20px",
             actionButton("act2","Run DEA", class= "btn-primary", width = "70%", 
                          style = "font-size:24px ; font-weight:600 ; border-radius:15px")
           )),
    
    ### 3- Results   ------------------------------------------------
    column(width = 9,
           navbarPage(
             title = tags$span(
               style = "color:gold ; font-weight:600",
               "DEA Results"
             ),
             
             #### p1: Datatable     ----------------------------------------------
             tabPanel(
               "Table",
               fluidRow(
                 dataTableOutput("tab2"),
                 div(
                   style = "margin-top:30px ; margin-bottom:30px",
                   downloadButton("download.tab","Download as excel", class = "btn-sm",
                                  style = "font-size:20px ; background-color:darkgreen ; padding:5px 150px ; border-radius:10px")
                 )
               )
             ),
             
             #### p2: Bar plot    ----------------------------------------------------
             tabPanel(
               "BarPlot",
               fluidRow(
                 column(width = 3,
                        div(
                          style = "margin-bottom:40px",
                          tags$span(
                            style = "font-weight:600 ; font-size:110%",
                            "Up genes to display",
                            icon("info-circle", id= "info-numeric4", style = "color:midnightblue ; cursor:pointer")
                          ),
                          numericInput("nUp",NULL,value = 10, width = "100px")
                        ),
                        bsTooltip(id = "info-numeric4",
                                  title = "From the most up-regulated to the least",
                                  placement = "top", trigger = "hover", options = list(container = "body"))),
                 column(width = 6,
                        div(
                          style = "margin-bottom:40px",
                          tags$span(
                            style = "font-weight:600 ; font-size:110%",
                            "Down genes to display",
                            icon("info-circle", id= "info-numeric5", style = "color:midnightblue ; cursor:pointer")
                          ),
                          numericInput("nDown",NULL,value = 10, width = "100px")
                        ),
                        bsTooltip(id = "info-numeric5",
                                  title = "From the most down-regulated to the least",
                                  placement = "top", trigger = "hover", options = list(container = "body"))),
                 
                 plotOutput("plt3", width = "100%", height = "auto"),
                 div(
                   style = "margin-top:30px ; margin-bottom:30px",
                   actionButton("download.bplot","Download as pdf", 
                                class = "btn-sm", icon = icon("download"),
                                style = "font-size:20px ; background-color:darkgreen ; padding:5px 150px ; border-radius:10px")
                 )
               )
             ),
             
             #### p3: Dot plot   ----------------------------------------------------
             tabPanel(
               "DotPlot",
               fluidRow(
                 column(width = 4,
                        div(
                          style = "margin-bottom:40px",
                          tags$span(
                            style = "font-weight:600 ; font-size:110%",
                            "Up genes",
                            icon("info-circle", id= "info-sel1", style = "color:midnightblue ; cursor:pointer")
                          ),
                          selectInput("Sel.up.dplot",NULL, choices = c(), multiple = T, selectize = F, width = "90%", size = 10)
                        ),
                        bsTooltip(id = "info-sel1",
                                  title = "From the most up-regulated to the least",
                                  placement = "top", trigger = "hover", options = list(container = "body"))),
                 column(width = 4,
                        div(
                          style = "margin-bottom:40px",
                          tags$span(
                            style = "font-weight:600 ; font-size:110%",
                            "Down genes",
                            icon("info-circle", id= "info-sel2", style = "color:midnightblue ; cursor:pointer")
                          ),
                          selectInput("Sel.down.dplot",NULL, choices = c(), multiple = T, selectize = F, width = "90%", size = 10)
                        ),
                        bsTooltip(id = "info-sel2",
                                  title = "From the most down-regulated to the least",
                                  placement = "top", trigger = "hover", options = list(container = "body"))),
                 column(width = 4,
                        div(
                          style = "margin-bottom:40px",
                          tags$span(
                            style = "font-weight:600 ; font-size:110%",
                            "Select metadata",
                            icon("info-circle", id= "info-sel3", style = "color:midnightblue ; cursor:pointer")
                          ),
                          selectInput("Sel.metadata2",NULL, choices = c(), width = "90%"),
                          tags$span(
                            style = "font-weight:600 ; font-size:110%",
                            "Specify xlab angle"
                          ),
                          sliderInput("angle", NULL, value = 0, min = 0, max = 90, step = 5, width = "90%"),
                          actionButton("display.dplot", "Display", class = "btn-sm", width = "90%",
                                       style = "font-size:16px ; background-color:midnightblue ; 
                                       padding:5px 40px ; border-radius:10px ; margin-top:20px ; border-color:gold")
                        ),
                        bsTooltip(id = "info-sel3",
                                  title = "Select the group of clusters where you want to see the expression of the selected genes",
                                  placement = "top", trigger = "hover", options = list(container = "body"))),
                 
                 plotOutput("plt4", width = "100%", height = "auto"),
                 div(
                   style = "margin-top:30px ; margin-bottom:30px",
                   actionButton("download.dplot","Download as pdf", 
                                class = "btn-sm", icon = icon("download"),
                                style = "font-size:20px ; background-color:darkgreen ; padding:5px 150px ; border-radius:10px")
                 )
               )
             ),
             
             #### p4: Heatmap   -----------------------------------------------------
             tabPanel(
               "Heatmap",
               fluidRow(
                 div(
                   style = "margin-bottom:10px",
                   tags$span(
                     style = "font-weight:600; font-size:130%; color:darkred; 
                            border-bottom:4px solid darkgreen; padding-bottom: 4px",
                     "Heatmap annotation :"
                   )
                 )
               ),
               
               fluidRow(
                 column(width = 4,
                        div(
                          style = "margin-bottom:20px",
                          tags$span(
                            style = "font-weight:600 ; font-size:110%",
                            "Up genes",
                            icon("info-circle", id= "info-selHmUp", style = "color:midnightblue ; cursor:pointer")
                          ),
                          selectInput("Sel.up.hmplot",NULL, choices = c(), multiple = T, selectize = F, width = "90%", size = 15)
                        ),
                        bsTooltip(id = "info-selHmUp",
                                  title = "From the most up-regulated to the least",
                                  placement = "top", trigger = "hover", options = list(container = "body"))),
                 column(width = 4,
                        div(
                          style = "margin-bottom:20px",
                          tags$span(
                            style = "font-weight:600 ; font-size:110%",
                            "Down genes",
                            icon("info-circle", id= "info-selHmDown", style = "color:midnightblue ; cursor:pointer")
                          ),
                          selectInput("Sel.down.hmplot",NULL, choices = c(), multiple = T, selectize = F, width = "90%", size = 15)
                        ),
                        bsTooltip(id = "info-selHmDown",
                                  title = "From the most down-regulated to the least",
                                  placement = "top", trigger = "hover", options = list(container = "body"))),
                 column(width = 4,
                        div(
                          style= "margin-bottom:20px",
                          tags$span(
                            style= "font-weight:600; font-size:110%",
                            "Select annotation",
                            icon("info-circle", id= "info-selHmAnnot", style= "color:midnightblue; cursor:pointer")
                          ),
                          selectInput("Sel.HmAnnot",NULL, choices = c(), multiple = T, selectize = F, width = "90%", size = 6),
                          actionButton("GenerateData.hmAnnot", "Upload the data", icon = icon("hand-point-right"),
                                       style = "font-size:16px ; background-color:purple ;  
                                       padding:5px 40px ; border-radius:10px ; margin-top:10px ; border-color:black"),
                          p("Then", style = "font-weight:600; font-size:14px; margin-top:10px"),
                          actionButton("Setup.hm", "Generate the plot", icon = icon("hand-point-right"),
                                       style = "font-size:16px ; background-color:purple ;  
                                       padding:5px 40px ; border-radius:10px ; border-color:black"),
                          p("Then", style = "font-weight:600; font-size:14px; margin-top:15px"),
                          actionButton("plot.hm", "Display the plot", class = "btn-sm", width = "90%",
                                       style = "font-size:21px ; background-color:midnightblue ; 
                                       padding:5px 40px ; border-radius:10px ; margin-top:10px ; border-color:red")
                        ),
                        bsTooltip(id = "info-selHmAnnot",
                                  title = "Select the metadata you want to annotate your heatmap with",
                                  placement = "top", trigger = "hover", options = list(container = "body")))
               ),
               
               fluidRow(
                 div(
                   style = "margin-bottom:30px; margin-top:20px",
                   tags$span(
                     style = "font-weight:600; font-size:130%; color:darkred; 
                            border-bottom:4px solid darkgreen; padding-bottom: 4px",
                     "Heatmap plot :"
                   )
                 ),
                 div(
                  style = "overflow-x:auto; height:auto; width:100%",
                  plotOutput("plt5", width = "auto", height = "auto") 
                 ),
                 div(
                   style = "margin-top:30px; margin-bottom:30px",
                   actionButton("download.hmplot", "Download as pdf",
                                class = "btn-sm", icon = icon("download"),
                                style = "font-size:20px ; background-color:darkgreen ; padding:5px 150px ; border-radius:10px")
                 )
               )
             )
             
           ))
  )
)















#  Shiny server     ###################
server <- function(input, output, session){
  # I- Importing data   ----------------------------------------------------
  seurat <- reactive({
    req(input$LoadSeurat)
    seurat <- readRDS(input$LoadSeurat$datapath)
  })
  
  output$InfoSeurat <- renderPrint({
    req(seurat())
    print(seurat())
  })
  
  
  # II- Metadata selection    --------------------------------------------------
  observeEvent(seurat(), {
    req(seurat())
    seurat <- seurat()
    clust_cols <- names(seurat@meta.data)[sapply(seurat@meta.data, function(x) !is.numeric(x))]
    updateSelectInput(session, "Sel.metadata", choices = clust_cols, selected = "")
    updateSelectInput(session, "Sel.metadata2", choices = clust_cols, selected = "")
  })
  
  summary_table <- eventReactive(input$act1, {
    req(seurat(), input$Sel.metadata)
    seurat <- seurat()
    seurat@meta.data %>% 
      dplyr::group_by(!!sym(input$Sel.metadata)) %>% 
      dplyr::summarise(`Nb of cells` = n())
  })
  
  output$tab1 <- renderTable({
    req(summary_table())
    summary_table()
  })
  
  output$plt1 <- renderPlot({
    req(seurat(), input$Sel.metadata)
    seurat <- seurat()
    DimPlot(object = seurat,
            group.by = input$Sel.metadata,
            cols = MyPalette,
            pt.size = 1.1)+
      my_umap_theme()+
      guides(colour = guide_legend(override.aes = list(size=3.5)))
  }, height = 500)
  
  
  # III- DEA    -----------------------------------------------------------------
  ## 1- Making the contrasts    --------------------------------------------------
  observeEvent(input$Sel.metadata, {
    req(seurat(), input$Sel.metadata)
    seurat <- seurat()
    cluster_levels <- unique(seurat@meta.data[[input$Sel.metadata]])
    updateSelectInput(session, "ident.1", choices = cluster_levels, selected = "")
    updateSelectInput(session, "ident.2", choices = cluster_levels, selected = "")
  })
  
  output$plt2 <- renderPlot({
    req(seurat(), input$ident.1, input$ident.2)
    seurat <- seurat()
    seurat@meta.data %>% 
      mutate(cond = ifelse(.data[[input$Sel.metadata]] == input$ident.1 | 
                             .data[[input$Sel.metadata]] == input$ident.2, T, F)) %>% 
      ggplot()+
      geom_point(aes(x= umap_1,
                     y= umap_2,
                     colour= .data[[input$Sel.metadata]],
                     alpha= cond),
                 size= 1.1)+
      labs(title = paste0(input$ident.1, "   vs   ", input$ident.2),
           colour= "Clusters",
           alpha= "")+
      scale_colour_manual(values = MyPalette)+
      scale_alpha_manual(values = c("TRUE"=1, "FALSE"=0.02), guide = "none")+
      my_umap_theme()+
      guides(colour = guide_legend(override.aes = list(size = 3.5)))
  }, height = 500)
  
  
  ## 2- Setting the parameters      ------------------------------------------------
  observeEvent(input$act2, {
    showNotification("DEA is running >>> ", duration = 5, type = "message")
    message("==============================================")
    message("DEA is running with the following parameters: ")
    message("==============================================")
    message(sprintf("%-25s: %s", "Selected metadata", input$Sel.metadata))
    message(sprintf("%-25s: %s", "Ident.1", input$ident.1))
    message(sprintf("%-25s: %s", "Ident.2", input$ident.2))
    message(sprintf("%-25s: %s", "Min.log2FC.threshold", input$log2fc))
    message(sprintf("%-25s: %s", "Min.pct.threshold", input$min.pct))
    message(sprintf("%-25s: %s", "Min.diff.pct", input$min.diff.pct))
    message("==============================================")
  })
  
  observeEvent(input$log2fc, {
    updateSliderInput(session, "log2fc_up", value = c(input$log2fc,10))
    updateSliderInput(session, "log2fc_down", value = c(-10,-(input$log2fc)))
  })
  
  
  ## 3- Running DEA     ----------------------------------------------------------
  DEA_results <- eventReactive(input$act2, {
    req(seurat())
    seurat <- seurat()
    Idents(seurat) <- as.factor(seurat@meta.data[[input$Sel.metadata]])
    seurat %>% 
      FindMarkers(ident.1 = input$ident.1,
                  ident.2 = input$ident.2,
                  only.pos = FALSE,
                  min.pct = input$min.pct,
                  min.diff.pct = input$min.diff.pct,
                  logfc.threshold = input$log2fc,
                  verbose = TRUE) %>% 
      dplyr::mutate(Diff_pct = round((pct.1 - pct.2), 3),
                    avg_log2FC = round(avg_log2FC, 3)) %>% 
      dplyr::filter(p_val_adj < 0.05,
                    avg_log2FC >= input$log2fc_up[1] & avg_log2FC <= input$log2fc_up[2] |
                      avg_log2FC >= input$log2fc_down[1] & avg_log2FC <= input$log2fc_down[2]) %>% 
      dplyr::arrange(desc(avg_log2FC)) %>%
      dplyr::select(avg_log2FC, pct.1, pct.2, Diff_pct, p_val_adj) %>%
      rownames_to_column(var = "Genes") 
  })
  
  ## 4- Initiate reactive values      ----------------------------------------------
  Up_Genes <- reactiveVal(character(0))
  Down_Genes <- reactiveVal(character(0))
  Dotplot_Genes <- reactiveVal(character(0))
  Hmplot_Genes <- reactiveVal(character(0))
  
  
  observeEvent(DEA_results(), {
    res <- DEA_results()
    max_fc <- max(res[["avg_log2FC"]])
    min_fc <- min(res[["avg_log2FC"]])
    message("______________________________________________________")
    message("Analysis completed successfully !!!")
    message("==============================================")
    message(sprintf("%-25s: %s", "Total Nb of DEGs", nrow(res)))
    message(sprintf("%-25s: %s", "Nb of Up genes", nrow(res[which(res$avg_log2FC > 0),])))
    message(sprintf("%-25s: %s", "Nb of Down genes", nrow(res[which(res$avg_log2FC < 0),])))
    message("==============================================")
    updateSliderInput(session, "log2fc_up", value = c(input$log2fc, max_fc))
    updateSliderInput(session, "log2fc_down", value = c(min_fc, -(input$log2fc)))
    # Fill the initiated vectors:
    Up_Genes(res %>% 
               dplyr::filter(avg_log2FC > 0) %>% 
               dplyr::pull(Genes))
    Down_Genes(res %>% 
                 dplyr::filter(avg_log2FC < 0) %>% 
                 dplyr::arrange(avg_log2FC) %>% 
                 dplyr::pull(Genes))
  })
  
  ## 5- Displaying results      ---------------------------------------------------
  
  ### p1: Datatable     -----------------------------------------------------------
  output$tab2 <- renderDataTable({
    req(DEA_results())
    DEA_results()
  }, options = list(pageLength = 5, scrollX = T))
  
  ### p2: Bar plot      -----------------------------------------------------------
  output$plt3 <- renderPlot({
    req(DEA_results(), Up_Genes(), Down_Genes())
    results <- DEA_results()
    up <- Up_Genes()
    down <- Down_Genes()
    
    rbind(results %>% 
            dplyr::filter(Genes %in% up) %>% 
            head(input$nUp),
          results %>% 
            dplyr::filter(Genes %in% down) %>% 
            tail(input$nDown)) %>% 
      ggplot()+
      geom_bar(aes(x= avg_log2FC,
                   y= reorder(Genes, avg_log2FC),
                   fill= p_val_adj),
               stat = "identity",
               width = 0.7)+
      theme_bw()+
      scale_fill_gradient(low = "navy", high = "darkred")+
      labs(title = "Top Up & Down genes",
           y= "")+
      theme(
        axis.title = element_text(face= "bold", size= 15),
        axis.text.y = element_text(size= 15, face = "bold"),
        axis.text.x = element_text(size= 13, face = "bold"),
        legend.title = element_text(face= "bold", size= 16, colour= "darkred", margin = margin(b=0.5, unit = "cm")),
        legend.text = element_text(size= 13, colour= "black"),
        plot.title = element_text(size = 18, colour = "darkred", face = "bold")
      )
  }, width = 750, height = 700)
  
  # Saving parameters:
  observeEvent(input$download.bplot, {
    showModal(
      modalDialog(
        title = "Save Barplot as PDF",
        numericInput("width_bplot", "Width (in inches):", value = 6, min = 4, step = 0.5),
        numericInput("height_bplot", "Height (in inches):", value = 8, min = 4, step = 0.5),
        textInput("bplot_name", "Filename:", value = ""),
        footer = tagList(
          modalButton("Cancel"),
          downloadButton("download.bplot.modal", "Download", class = "btn-success")
        )
      )
    )
  })
  
  ### p3: Dot plot      -----------------------------------------------------------
  
  # Update gene list selection:
  observeEvent(Up_Genes(),{
    updateSelectInput(session, "Sel.up.dplot", choices = Up_Genes())
  })
  observeEvent(Down_Genes(),{
    updateSelectInput(session, "Sel.down.dplot", choices = Down_Genes())
  })
  
  # Initiate the plot:
  Dplot <- eventReactive(input$display.dplot, {
    req(DEA_results(), seurat())
    Dotplot_Genes(c(input$Sel.up.dplot, input$Sel.down.dplot))
    seurat <- seurat()
    DotPlot(seurat,
            features = Dotplot_Genes(),
            group.by = input$Sel.metadata2)+
      RotatedAxis()+
      coord_flip()+
      theme_bw()+
      labs(title = paste0("Genes expression within :  ", input$Sel.metadata2),
           x="",
           y="")+
      scale_color_gradient(low = "orange", high = "navy")+
      theme(axis.title = element_text(face= "bold", size= 15),
            axis.text.y = element_text(size= 16, face = "bold"),
            axis.text.x = element_text(size= 16, face = "bold", angle = input$angle),
            legend.title = element_text(face= "bold", size= 16, colour= "darkred"),
            legend.text = element_text(size= 13, colour= "black"),
            plot.title = element_text(size = 18, colour = "darkred", face = "bold"),
            panel.grid = element_blank())
  })
  
  # Plot:
  output$plt4 <- renderPlot({
    req(Dplot())
    Dplot()
  }, width = 800, height = 700)
  
  # Saving parameters:
  observeEvent(input$download.dplot, {
    showModal(
      modalDialog(
        title = "Save Dotplot as PDF",
        numericInput("width_dplot", "Width (in inches):", value = 6, min = 4, step = 0.5),
        numericInput("height_dplot", "Height (in inches):", value = 8, min = 4, step = 0.5),
        textInput("dplot_name", "Filename:", value = ""),
        footer = tagList(
          modalButton("Cancel"),
          downloadButton("download.dplot.modal", "Download", class = "btn-success")
        )
      )
    )
  })
  
  ### p4: Heatmap      -----------------------------------------------------------
  
  # update gene list selection:
  observeEvent(Up_Genes(),{
    updateSelectInput(session, "Sel.up.hmplot", choices = Up_Genes())
  })
  observeEvent(Down_Genes(),{
    updateSelectInput(session, "Sel.down.hmplot", choices = Down_Genes())
  })
  
  # update heatmap annotation list:
  observeEvent(DEA_results(), {
    req(seurat())
    seurat <- seurat()
    meta <- names(seurat@meta.data)[sapply(seurat@meta.data, function(x) !is.numeric(x))]
    updateSelectInput(session, "Sel.HmAnnot", choices = meta)
  })
  
  #### Expression assay       -----------------------------------------------------
  ExpressionAssay <- eventReactive(input$GenerateData.hmAnnot, {
    req(seurat())
    Hmplot_Genes(c(input$Sel.up.hmplot, input$Sel.down.hmplot))
    seurat <- seurat()
    FetchData(seurat, vars = Hmplot_Genes()) %>% t() %>% as.matrix() %>% log1p()
  })
  
  observeEvent(ExpressionAssay(), {
    req(ExpressionAssay())
    showNotification("The expression assay has been created", duration = 5, type = "message")
    message(paste0("> Heatmap expression assay successfully created ! (", nrow(ExpressionAssay()), " genes selected)"))
  })
  
  #### Annotation data      -----------------------------------------------------------
  hm.metadata <- eventReactive(input$GenerateData.hmAnnot, {
    req(seurat())
    seurat <- seurat()
    selected_vars <- input$Sel.HmAnnot
    seurat@meta.data %>% 
      dplyr::select(all_of(selected_vars))
  })
  
  observeEvent(hm.metadata(), {
    req(hm.metadata())
    showNotification("The annotation data has been created", duration = 5, type = "message")
    message(paste0("> Heatmap annotation data successfully created ! (", ncol(hm.metadata()), " variable(s) used)"))
    print(head(hm.metadata(),4))
  })
  
  ##### Top annot       -----------------------------------------------------------
  TopAnnotation <- eventReactive(input$Setup.hm, {
    req(hm.metadata())
    hm.metadata <- hm.metadata()
    
    # Three cases: (1, 2 or 3 variables):
    if(ncol(hm.metadata) == 1){
      # color palette:
      col1 <- Assign_colors(hm.metadata[[input$Sel.HmAnnot]],
                            MyPalette)
      # Set annotation:
      HeatmapAnnotation(df = hm.metadata,
                        col = setNames(list(col1),
                                       c(input$Sel.HmAnnot)),
                        annotation_name_side = "left",
                        annotation_name_gp = list(fontsize = 13,
                                                  col = "navy",
                                                  fontface = "bold"),
                        annotation_legend_param = list(grid_height = unit(0.9,"cm"),
                                                       grid_width = unit(0.7,"cm"),
                                                       labels_gp = gpar(col = "black",
                                                                        fontsize = 12),
                                                       title_gp = gpar(col = "darkred",
                                                                       fontsize = 13,
                                                                       fontface = "bold")))
    } else if(ncol(hm.metadata) == 2){
      # color palette:
      col1 <- Assign_colors(hm.metadata[[input$Sel.HmAnnot[1]]],
                            MyPalette)
      col2 <- Assign_colors(hm.metadata[[input$Sel.HmAnnot[2]]],
                            MyPalette[which(!(MyPalette %in% col1))])
      # Set annotation:
      HeatmapAnnotation(df = hm.metadata,
                        col = setNames(list(col1, col2),
                                       c(input$Sel.HmAnnot[1],
                                         input$Sel.HmAnnot[2])),
                        annotation_name_side = "left",
                        annotation_name_gp = list(fontsize = 13,
                                                  col = "navy",
                                                  fontface = "bold"),
                        annotation_legend_param = list(grid_height = unit(0.9,"cm"),
                                                       grid_width = unit(0.7,"cm"),
                                                       labels_gp = gpar(col = "black",
                                                                        fontsize = 12),
                                                       title_gp = gpar(col = "darkred",
                                                                       fontsize = 13,
                                                                       fontface = "bold")))
    } else if(ncol(hm.metadata) >= 3){
      # color palette:
      col1 <- Assign_colors(hm.metadata[[input$Sel.HmAnnot[1]]],
                            MyPalette)
      col2 <- Assign_colors(hm.metadata[[input$Sel.HmAnnot[2]]],
                            MyPalette[which(!(MyPalette %in% col1))])
      col3 <- Assign_colors(hm.metadata[[input$Sel.HmAnnot[3]]],
                            MyPalette[which(!(MyPalette %in% c(col1, col2)))])
      # Set annotation:
      HeatmapAnnotation(df = hm.metadata,
                        col = setNames(list(col1, col2, col3),
                                       c(input$Sel.HmAnnot[1],
                                         input$Sel.HmAnnot[2],
                                         input$Sel.HmAnnot[3])),
                        annotation_name_side = "left",
                        annotation_name_gp = list(fontsize = 13,
                                                  col = "navy",
                                                  fontface = "bold"),
                        annotation_legend_param = list(grid_height = unit(0.9,"cm"),
                                                       grid_width = unit(0.7,"cm"),
                                                       labels_gp = gpar(col = "black",
                                                                        fontsize = 12),
                                                       title_gp = gpar(col = "darkred",
                                                                       fontsize = 13,
                                                                       fontface = "bold")))
      
      warning("NOTE: 3 variables is a maximum possible to annotate with, any more would be ignored !")
    }
  })
  
  observeEvent(TopAnnotation(), {
    req(TopAnnotation())
    message("> Heatmap top-annotation is ready !")
  })
  
  ##### Right annot       ---------------------------------------------------------
  RightAnnotation <- eventReactive(input$Setup.hm, {
    req(ExpressionAssay())
    exp <- ExpressionAssay()
    dea.res <- DEA_results()
    rowAnnotation(
      "log2fc" = anno_barplot(dea.res$avg_log2FC[which(dea.res$Genes %in% rownames(exp))],
                              axis = T, border = T, cex = 1, bar_width = 1)
    )
  })
  
  observeEvent(RightAnnotation(), {
    req(RightAnnotation())
    message("> Heatmap right-annotation is ready !")
  })
  
  #### Plot heatmap       ---------------------------------------------------------
  hm <- eventReactive(input$plot.hm, {
    req(ExpressionAssay(), RightAnnotation(), TopAnnotation(), hm.metadata())
    exp <- ExpressionAssay()
    meta <- hm.metadata()
    RightAnnotation <- RightAnnotation()
    TopAnnotation <- TopAnnotation()
    
    # color gradient setting:
    ColExpr <- colorRamp2(breaks = c(0, max(exp)/4, max(exp)/2, max(exp)*4/5, max(exp)),
                          colors = c("#FFFFFF","#FFCC33","#DD8833","#660000","#220000"))
    
    showNotification("Heatmap on its way !", duration = 5, type = "message")
    
    # plot setting:
    ComplexHeatmap::Heatmap(exp,
                            top_annotation = TopAnnotation,
                            right_annotation = RightAnnotation,
                            cluster_rows = TRUE,
                            cluster_columns = TRUE,
                            show_row_names = T,
                            show_column_names = F,
                            col = ColExpr,
                            show_row_dend = F,
                            row_dend_side = "left",
                            row_dend_width = unit(2,"cm"),
                            show_column_dend = F,
                            row_names_side = "left",
                            row_names_gp = gpar(fontface = "bold",
                                                fontsize = 12,
                                                col = "darkred"),
                            column_title = paste0("Expression of ",nrow(exp), " genes across ", ncol(meta), " condition(s)"),
                            heatmap_legend_param = list(legend_height = unit(3,"cm"),
                                                        direction = "horizontal",
                                                        legend_width = unit(7,"cm"),
                                                        title = "log_expression",
                                                        title_position = "topcenter",
                                                        title_gp = gpar(fontsize=13, col="darkred", fontface="bold"),
                                                        label_gp = gpar(fontsize=10, col="black"),
                                                        legend_gp = gpar(fontsize=10)),
                            column_title_gp = gpar(col = "darkred",
                                                   fontsize = 15,
                                                   fontface = "bold"))
    
  })
  
  observeEvent(hm(), {
    req(hm())
    message("> Heatmap plot successfully created !")
  })
  
  output$plt5 <- renderPlot({
    req(hm())
    hm <- hm()
    ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom")
  }, width = 900, height = 700)
  
  
  
  
  
  # X- Downloading      -----------------------------------------------------------
  
  # Data table:
  output$download.tab <- downloadHandler(
    filename = function(){paste0("DEA_res_", input$ident.1,"_vs_",input$ident.2,".xlsx")},
    content = function(file){write_xlsx(DEA_results(), path = file)}
  )
  
  # Bar plot:
  output$download.bplot.modal <- downloadHandler(
    filename = function(){paste0(input$bplot_name,".pdf")},
    content = function(file){
      pdf(file, width = input$width_bplot, height = input$height_bplot)
      results <- DEA_results()
      up <- Up_Genes()
      down <- Down_Genes()
      rbind(results %>% 
              dplyr::filter(Genes %in% up) %>% 
              head(input$nUp),
            results %>% 
              dplyr::filter(Genes %in% down) %>% 
              tail(input$nDown)) %>% 
        ggplot()+
        geom_bar(aes(x= avg_log2FC,
                     y= reorder(Genes, avg_log2FC),
                     fill= p_val_adj),
                 stat = "identity",
                 width = 0.7)+
        theme_bw()+
        scale_fill_gradient(low = "navy", high = "darkred")+
        labs(title = "Top Up & Down genes",
             y= "")+
        theme(axis.title = element_text(face= "bold", size= 15),
              axis.text.y = element_text(size= 15, face = "bold"),
              axis.text.x = element_text(size= 13, face = "bold"),
              legend.title = element_text(face= "bold", size= 15, colour= "navy"),
              legend.text = element_text(size= 13, colour= "black"),
              plot.title = element_text(size = 17, colour = "darkred", face = "bold")) -> plot
      
      print(plot)
      dev.off()
    }
  )
  
  # Dot plot:
  output$download.dplot.modal <- downloadHandler(
    filename = function(){paste0(input$dplot_name,".pdf")},
    content = function(file){
      pdf(file, width = input$width_dplot, height = input$height_dplot)
      Dotplot_Genes(c(input$Sel.up.dplot, input$Sel.down.dplot))
      seurat <- seurat()
      DotPlot(seurat,
              features = Dotplot_Genes(),
              group.by = input$Sel.metadata2)+
        RotatedAxis()+
        coord_flip()+
        theme_bw()+
        labs(title = paste0("Genes expression within :  ", input$Sel.metadata2),
             x="",
             y="")+
        scale_color_gradient(low = "orange", high = "navy")+
        theme(axis.title = element_text(face= "bold", size= 15),
              axis.text.y = element_text(size= 16, face = "bold"),
              axis.text.x = element_text(size= 16, face = "bold", angle = input$angle),
              legend.title = element_text(face= "bold", size= 16, colour= "darkred"),
              legend.text = element_text(size= 13, colour= "black"),
              plot.title = element_text(size = 18, colour = "darkred", face = "bold"),
              panel.grid = element_blank()) -> plot
      
      print(plot)
      dev.off()
    }
  )
  
  
}




#  ShinyApp run     ########################
shinyApp(ui, server)














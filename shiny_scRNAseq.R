#   Setting     #################
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(writexl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ComplexHeatmap)
library(DT)
library(stringr)
library(shiny)
library(shinythemes)
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
  theme = bs_theme(bootswatch = "flatly"),
  h2("I- Importing the data",
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
  h2("II- Metadata selection",
     style = "color:gold ; font-weight:700 ; background-color:black ; margin-top:150px ; margin-bottom:20px"),
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
           actionButton("act1", "View the summary", class = "btn-lg btn-primary", width = "100%"),
           p("", style = "margin-top:30px"),
           tableOutput("tab1"))
  ),
  
  ## III- DEA   ------------------------------------------------------------------
  h2("III- Differential expression analysis",
     style = "color:gold ; font-weight:700 ; background-color:black ; margin-top:150px ; margin-bottom:20px"),
  
  ### 1- Making the contrasts   --------------------------------------------------
  p("1- Select the contrasts",
    style = "color:darkred ; font-weight:600 ; font-size:130% ; background-color:gold ; margin-bottom:20px"),
  fluidRow(
    column(width = 4,
           p("From the selected metadata above, choose the cluster you want to test and againt which cluster ! 
             (e.g Relapse vs Diagnosis)",
             style = "font-size:110% ; font-weight:600 ; color:darkgreen"),
           selectInput("ident.1", "Selection 1", choices = c(), width = "300px"),
           selectInput("ident.2", "Selection 2", choices = c(), width = "300px")),
    column(width = 6,
           offset = 1,
           plotOutput("plt2", height = "500px"))
  ),
  
  ### 2- Setting the parameters   ------------------------------------------------
  p("2- Setting the parameters (prior DEA)",
    style = "color:darkred ; font-weight:600 ; font-size:130% ; background-color:gold ; margin-top: 20px ; margin-bottom:20px"),
  fluidRow(
    column(width = 2,
           sliderInput("log2fc", "Logfc min.threshold", min = 0, max = 3, step = 0.1, value = 0.5, width = "200px"),
           p("Define the minimum Log2FC threshold to be kept in an absolute value",
             style = "color:darkgray ; font-style:italic ; font-size:90%"),
           sliderInput("min.pct", "Min.pct threshold", min = 0, max = 1, step = 0.05, value = 0.1, width = "200px"),
           p("Define the minimum percentage of cells within any cluster that have to express the gene",
             style = "color:darkgray ; font-style:italic ; font-size:90%")),
    column(width = 2.5,
           sliderInput("logfc_filt_up", "Filter Up log2FC", min = 0, max = 12, value = c(0,12), step = 0.1, width = "250px"),
           sliderInput("logfc_filt_down", "Filter down log2FC", min = -12, max = 0, value = c(-12,0), step = 0.1, width = "250px"),
           sliderInput("pct.1_filt", "Filter pct.1", min = 0, max = 1, value = c(0,1), width = "250px"))
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
      labs(title = paste0(input$ident.1, "  _vs_  ", input$ident.2),
           colour= "Clusters",
           alpha= "")+
      scale_colour_manual(values = MyPalette)+
      scale_alpha_manual(values = c("TRUE"=1, "FALSE"=0.02), guide = "none")+
      my_umap_theme()+
      guides(colour = guide_legend(override.aes = list(size = 3.5)))
  }, height = 500)
  
}




#  ShinyApp run     ########################
shinyApp(ui, server)














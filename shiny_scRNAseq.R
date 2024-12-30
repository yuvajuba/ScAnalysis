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
           actionButton("act1", "View the summary", class = "btn-lg btn-primary", 
                        width = "100%", style = "font-size:28px ; color:gold ; font-weight:600 ; border-radius:15px"),
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
           p("Initial parameters", style = "color:darkgreen ; font-weight:600 ; font-size:120% ; margin-bottom:30px"),
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
             style = "color:darkgreen ; font-weight:600 ; font-size:120% ; margin-bottom:30px ; margin-top:20px"),
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
             style = "margin-top:20px ; margin-left:20px ; margin-bottom:20px",
             actionButton("act2","Run DEA", class= "btn-primary", width = "70%", 
                          style = "font-size:26px ; font-weight:600 ; border-radius:10px")
           )),
    
    ### 3- Results   ------------------------------------------------
    column(width = 9,
           navbarPage(
             "DEA Results",
             
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
  
  
  ## 2- Setting the parameters      ------------------------------------------------
  observeEvent(input$act2, {
    showNotification("DEA is running >>> ", duration = 10, type = "message")
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
    updateSliderInput(session, "log2fc_up", value = c(input$log2fc, max_fc))
    updateSliderInput(session, "log2fc_down", value = c(min_fc, -(input$log2fc)))
  })
  
  ## 4- Displaying results      ---------------------------------------------------
  ### p1: datatable     -----------------------------------------------------------
  output$tab2 <- renderDataTable({
    req(DEA_results())
    DEA_results()
  }, options = list(pageLength = 5, scrollX = T))
  
  
  
  
  
  
  
  
  
  
  # X- Downloading      -----------------------------------------------------------
  output$download.tab <- downloadHandler(
    filename = function(){paste0("DEA_res_", input$ident.1,"_vs_",input$ident.2,".xlsx")},
    content = function(file){write_xlsx(DEA_results(), path = file)}
  )
  
  
  
}




#  ShinyApp run     ########################
shinyApp(ui, server)














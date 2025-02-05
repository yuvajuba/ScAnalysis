#   Setting     #################
library(Seurat)
library(yulab.utils)
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
library(ggupset)

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


Display_Venn <- function(Markers, 
                         colpalette = NULL, 
                         set.names = NULL, 
                         set.name.size = 5,
                         text.size = 5,
                         Padding = 0.03) {
  
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
                                                                                     x[[s]][[3]]), 
                                                                                intersect)
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
        specific_genes[[paste(names(x[[s]]), collapse = " | ")]] <- setNames(list(setdiff(x[[s]][[1]], 
                                                                                          intersect(x[[s]][[1]], 
                                                                                                    x[[s]][[2]])),
                                                                                  setdiff(x[[s]][[2]], 
                                                                                          intersect(x[[s]][[1]], 
                                                                                                    x[[s]][[2]]))),
                                                                             c(names(x[[s]][1]), 
                                                                               names(x[[s]][2])))
      }
      
      if(i == 3){
        specific_genes[[paste(names(x[[s]]), collapse = " | ")]] <- setNames(list(setdiff(x[[s]][[1]], 
                                                                                          union(x[[s]][[2]], 
                                                                                                x[[s]][[3]])),
                                                                                  setdiff(x[[s]][[2]], 
                                                                                          union(x[[s]][[1]], 
                                                                                                x[[s]][[3]])),
                                                                                  setdiff(x[[s]][[3]], 
                                                                                          union(x[[s]][[2]], 
                                                                                                x[[s]][[1]]))),
                                                                             c(names(x[[s]][1]), 
                                                                               names(x[[s]][2]), 
                                                                               names(x[[s]][3])))
      }
      
      if(i == 4){
        specific_genes[["unique each"]] <- setNames(list(setdiff(x[[s]][[1]],
                                                                 purrr::reduce(list(x[[s]][[2]], 
                                                                                    x[[s]][[3]], 
                                                                                    x[[s]][[4]]), 
                                                                               union)),
                                                         setdiff(x[[s]][[2]],
                                                                 purrr::reduce(list(x[[s]][[1]], 
                                                                                    x[[s]][[3]], 
                                                                                    x[[s]][[4]]), 
                                                                               union)),
                                                         setdiff(x[[s]][[3]],
                                                                 purrr::reduce(list(x[[s]][[2]], 
                                                                                    x[[s]][[1]], 
                                                                                    x[[s]][[4]]), 
                                                                               union)),
                                                         setdiff(x[[s]][[4]],
                                                                 purrr::reduce(list(x[[s]][[2]], 
                                                                                    x[[s]][[3]], 
                                                                                    x[[s]][[1]]), 
                                                                               union))),
                                                    c(names(x[[s]][1]),
                                                      names(x[[s]][2]),
                                                      names(x[[s]][3]),
                                                      names(x[[s]][4])))
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
    set_name_size = set.name.size,
    text_color = "black",
    text_size = text.size,
    padding = Padding, 
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
          plot.title = element_text(size = 20, face = "bold", colour = "darkred", 
                                    hjust = 0.5, margin = margin(b= 0.2, unit = "cm")),
          plot.subtitle = element_text(size = 17, face = "bold", colour = "#771144", 
                                       hjust = 0.5, margin = margin(b= 0.5, unit = "cm")),
          legend.title = element_text(size = 17, face = "bold", colour = "darkred", 
                                      hjust = 0.5, margin = margin(b= 0.5, unit = "cm")),
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
    column(width = 4,
           offset = 1,
           actionLink("act1", "View the summary", icon = icon("hand-point-right"),
                      style = "font-size:24px ; color:purple ; font-weight:600"),
           p("", style = "margin-bottom:30px"),
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
    style = "color:darkred ; font-weight:600 ; font-size:160% ; 
    background-color:gold ; margin-top: 50px ; margin-bottom:20px"),
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
               style = "color:gold ; font-weight:600 ; font-size:120%",
               "DEA Results"
             ), 
             collapsible = T,
             
             #### p1: Datatable     ----------------------------------------------
             tabPanel(
               tags$span(
                 style = "color:white ; font-weight:600 ; font-size:120%",
                 "Table"
               ),
               fluidRow(
                 dataTableOutput("tab2"),
                 div(
                   style = "margin-top:30px ; margin-bottom:30px",
                   downloadButton("download.tab","Download as excel", class = "btn-sm",
                                  style = "font-size:20px ; background-color:darkgreen ; 
                                  padding:5px 150px ; border-radius:10px")
                 )
               )
             ),
             
             #### p2: Bar plot    ----------------------------------------------------
             tabPanel(
               tags$span(
                 style = "color:white ; font-weight:600 ; font-size:120%",
                 "BarPlot"
               ),
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
                                style = "font-size:20px ; background-color:darkgreen ; 
                                padding:5px 150px ; border-radius:10px")
                 )
               )
             ),
             
             #### p3: Dot plot   ----------------------------------------------------
             tabPanel(
               tags$span(
                 style = "color:white ; font-weight:600 ; font-size:120%",
                 "DotPlot"
               ),
               fluidRow(
                 column(width = 4,
                        div(
                          style = "margin-bottom:40px",
                          tags$span(
                            style = "font-weight:600 ; font-size:110%",
                            "Up genes",
                            icon("info-circle", id= "info-sel1", style = "color:midnightblue ; cursor:pointer")
                          ),
                          selectInput("Sel.up.dplot",NULL, choices = c(), 
                                      multiple = T, selectize = F, width = "90%", size = 10)
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
                          selectInput("Sel.down.dplot",NULL, choices = c(), 
                                      multiple = T, selectize = F, width = "90%", size = 10)
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
                                style = "font-size:20px ; background-color:darkgreen ; 
                                padding:5px 150px ; border-radius:10px")
                 )
               )
             ),
             
             #### p4: Heatmap   -----------------------------------------------------
             tabPanel(
               tags$span(
                 style = "color:white ; font-weight:600 ; font-size:120%",
                 "Heatmap"
               ),
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
                          selectInput("Sel.up.hmplot",NULL, choices = c(), 
                                      multiple = T, selectize = F, width = "90%", size = 13)
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
                          selectInput("Sel.down.hmplot",NULL, choices = c(), 
                                      multiple = T, selectize = F, width = "90%", size = 13)
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
                          selectInput("Sel.HmAnnot",NULL, choices = c(), 
                                      multiple = T, selectize = F, width = "90%", size = 7),
                          actionLink("GenerateData.hmAnnot", "Upload the data", icon = icon("hand-point-right"),
                                     style = "font-size:16px ; margin-top:10px ; color:purple ; font-weight:600"),
                          p("Then", style = "font-weight:600; font-size:14px; margin-top:10px"),
                          actionLink("Setup.hm", "Generate the plot", icon = icon("hand-point-right"),
                                     style = "font-size:16px; color:purple; margin-top:10px; 
                                     font-weight:600; margin-bottom:20px")
                        ),
                        bsTooltip(id = "info-selHmAnnot",
                                  title = "Select the metadata you want to annotate your heatmap with. You can select up to 3 !",
                                  placement = "top", trigger = "hover", options = list(container = "body")),
                        textOutput("hm.ready"),
                        actionButton("plot.hm", "Display the plot", class = "btn-sm", width = "90%",
                                     style = "font-size:20px; background-color:midnightblue; font-weight:600
                                     padding:5px 40px; border-radius:10px; margin-top:20px; border-color:cadetblue"))
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
                                style = "font-size:20px ; background-color:darkgreen ; 
                                padding:5px 150px ; border-radius:10px")
                 )
               )
             )
             
           ))
  ),
  
  ## IV- Markers inspection   ---------------------------------------------------------
  h1("IV- Markers inspection",
     style = "color:gold ; font-weight:700 ; background-color:black ; margin-top:200px ; margin-bottom:20px"),
  
  ### 1- Displaying the list of genes     ------------------------------------------
  p("1- Displaying the list of genes",
    style = "color:darkred ; font-weight:600 ; font-size:160% ; background-color:gold ; margin-bottom:20px"),
  fluidRow(
    column(width = 6,
           div(
             style = "margin-bottom:20px",
             tags$span(
               style = "font-weight:600; font-size:130%; color:darkgreen; 
                            border-bottom:4px solid darkred; padding-bottom: 4px",
               "Up regulated genes: "
             )
           ),
           radioButtons("show_upgenes", NULL, choices = c("Top10", "Top50", "All"), 
                        inline = T, selected = ""),
           textOutput("genes1")),
    column(width = 6,
           div(
             style = "margin-bottom:20px",
             tags$span(
               style = "font-weight:600; font-size:130%; color:darkgreen; 
                            border-bottom:4px solid darkred; padding-bottom: 4px",
               "Down regulated genes: "
             )
           ),
           radioButtons("show_downgenes", NULL, choices = c("Top10", "Top50", "All"), 
                        inline = T, selected = ""),
           textOutput("genes2"))
  ),
  
  fluidRow(
    column(width = 6,
           div(
             style = "margin-bottom:20px; margin-top:20px",
             tags$span(
               style = "font-weight:600; font-size:130%; color:darkgreen; 
                            border-bottom:4px solid darkred; padding-bottom: 4px",
               "Save the list"
             )
           ),
           p("Add the list of genes from this DEA into a cumulative list of markers to be 
             compared later",
             style = "font-size:110%; font-weight:600; margin-bottom:10px"),
           textInput("name.list", "Name your list"),
           actionLink("add.markers", "Click to add these markers", icon = icon("hand-point-right"),
                      style = "font-size:16px; color:purple; font-weight:600; margin-bottom:20px")),
    column(width = 6,
           p("", style = "margin-top:20px"),
           verbatimTextOutput("list_markers"),
           downloadButton("download.markers", "Download the list",
                          icon = icon("download"),
                          style = "font-size:20px ; background-color:darkgreen ; 
                                  padding:5px 150px ; border-radius:10px ; margin-top:10px"),
           actionButton("ac","Clear you list", class = "btn-sm btn-danger",
                        style = "margin-top:20px; margin-bottom:30px; padding:4px 80px ; 
                        border-radius:10px; font-size:18px"))
  ),
  
  ### 2- Compare markers     ----------------------------------------------------
  p("2- Compare the markers",
    style = "color:darkred ; font-weight:600 ; font-size:160% ; background-color:gold ; 
    margin-bottom:20px ; margin-top:30px"),
  p("Let's find the common and specific genes among the contrasts !",
    style = "color:darkgreen ; font-weight:600 ; font-size:100%"),
  
  fluidRow(
    column(width = 3,
           div(
             style = "margin-top:20px; margin-bottom:20px",
             tags$span(
               style = "font-size:110%; font-weight:600",
               "Select from the current list",
               icon("info-circle", id= "info-selMarkers", style= "color:midnightblue; cursor:pointer")
             ),
             selectInput("sel.markers.list", NULL, choices = c(), multiple = T,
                         selectize = F, size = 7, width = "300px")
           ),
           bsTooltip(id = "info-selMarkers",
                     title = "You can make up to 4 comparison (any selection over the 4th will be ignored)",
                     placement = "top", trigger = "hover", options = list(container = "body")),
           radioButtons("up_down", "Which genes to compare ?", 
                        choices = c("Up genes","Down genes","All genes"), inline = F, selected = ""),
           actionButton("display_venn", "Display Venn diagram", class = "btn-sm", width = "90%",
                        style = "font-size:18px; background-color:midnightblue; font-weight:600; margin-bottom:30px;
                        border-radius:10px; margin-top:10px; border-color:cadetblue")),
    column(width = 8,
           offset = 1,
           div(
             style = "margin-top:20px; margin-bottom:30px",
             tags$span(
               style = "font-size:120%; font-weight:600; color:darkgreen; border-bottom:4px solid darkred; margin-bottom:20px",
               "Venn plot:"
             )
           ),
           column(width = 4,
                  numericInput("vplot.setname.size", "Set name size", value = 7, width = "150px")),
           column(width = 4,
                  numericInput("vplot.text.size", "Text size", value = 7, width = "150px")),
           column(width = 4,
                  numericInput("vplot.padding", "Padding", value = 0.03, step = 0.01, width = "150px")),
           div(
             style = "overflow-x:auto; width:100%; margin-top:20px",
             plotOutput("vennplot", height = "auto")
           ),
           div(
             style = "margin-top:30px; margin-bottom:30px",
             actionButton("download.vplot", "Download as pdf",
                          class = "btn-sm", icon = icon("download"),
                          style = "font-size:20px ; background-color:darkgreen ; 
                                padding:5px 150px ; border-radius:10px")
           ),
           div(
             style = "margin-top:20px; margin-bottom:20px",
             tags$span(
               style = "font-size:120%; font-weight:600; color:darkgreen; 
             border-bottom:4px solid darkred; margin-top:20px; margin-bottom:20px",
             "Common genes:"
             )
           ),
           selectInput("sel.intersection", "Select you intersection", choices = c(), width = "650px"),
           verbatimTextOutput("nb_common_genes"),
           verbatimTextOutput("common_genes"),
           div(
             style = "margin-top:30px; margin-bottom:20px",
             tags$span(
               style = "font-size:120%; font-weight:600; color:darkgreen; 
             border-bottom:4px solid darkred; margin-top:20px; margin-bottom:20px",
             "Group-specific genes:"
             )
           ),
           selectInput("specific.genes", "Select the contrast", choices = c(), width = "650px"),
           selectInput("specific.genes2", "Select the condition", choices = c(), width = "650px"),
           verbatimTextOutput("nb_genes"),
           verbatimTextOutput("specific_genes"))
  ),
  
  ## V- PEA   ---------------------------------------------------------
  h1("V- Enrichment Analysis",
     style = "color:gold ; font-weight:700 ; background-color:black ; margin-top:200px ; margin-bottom:20px"),
  p("Now we will apply an enrichment analysis on an individual set of markers and on a group of sets to 
    compare between them, from the list of genes filled above.",
    style = "font-size:130% ; font-weight:600 ; color:navy ; margin-bottom:20px"),
  
  ### 1- Individual Enrichment Analysis     -----------------------------------
  p("1- Individual analysis",
    style = "color:darkred ; font-weight:600 ; font-size:160% ; background-color:gold ; 
    margin-bottom:10px ; margin-top:30px"),
  
  #### 1-1- GO     -----------------------------------------------------------
  p("1-1 Pathway Enrichment Analysis",
    style = "color:white ; font-weight:600 ; font-size:150% ; background-color:red ; 
    margin-bottom:10px ; margin-top:10px"),
  p("For PEA, we'll use the GO database to map the terms with the genes",
    style = "font-size:120% ; font-weight:600 ; color:black ; margin-bottom:20px"),
  
  navbarPage(
    title = tags$span(
      style = "color:gold ; font-weight:600 ; font-size:120%",
      "GO results"
    ), 
    collapsible = T,
    
    ##### p1: Datatable       ------------------------------------------------
    tabPanel(
      tags$span(
        style = "color:white ; font-weight:600 ; font-size:120%",
        "Run GO"
      ),
      
      fluidRow(
        column(width = 3,
               selectInput("go.list.genes", "Select gene list", choices = c(), width = "200px"),
               selectInput("go.list.genes.reg", "Select the regulation", choices = c("Up","Down","All"), width = "200px"),
               selectInput("go.keytype", "Key Type", choices = c("SYMBOL","ENTREZID","ENSEMBL"), width = "200px"),
               selectInput("go.ont", "Ontology", choices = c("BP","MF", "CC"), selected = "BP", width = "200px"),
               sliderInput("go.maxgeneset", "Limit the geneset size", value = 500,
                           min = 500, max = 10000, step = 200, width = "300px"),
               actionButton("run.go", "Run GO analysis", class = "btn-sm", width = "90%",
                            style = "font-size:18px; background-color:midnightblue; font-weight:600; margin-bottom:30px;
                            border-radius:10px; margin-top:10px; border-color:cadetblue")),
        column(width = 9,
               dataTableOutput("tab3"),
               div(
                 style = "margin-top:20px; margin-bottom:40px",
                 downloadButton("download.tab.go","Download as excel", class = "btn-sm",
                                style = "font-size:20px ; background-color:darkgreen ; 
                                padding:5px 150px ; border-radius:10px")
               ),
               column(width = 4,
                      actionLink("addgo", "click to save the result into a list for a comparative analysis",
                                 icon = icon("hand-point-right"), 
                                 style = "font-size:120%; font-weight:600"),
                      p("", style = "margin-bottom:30px"),
                      actionLink("clear.go.list", "Clear the list", icon = icon("eraser"),
                                 style = "font-size:120%; font-weight:600; color:darkred")),
               column(width = 6,
                      offset = 1,
                      verbatimTextOutput("showGOobj", placeholder = T),
                      downloadButton("dwnload.GoObjList", "Download the list", class = "btn-sm btn-success", 
                                     style = "font-size:14px; border-radius:10px; padding:5px 50px; font-weight:600")))
      )
    ),
    
    ##### p2: Dot plot        --------------------------------------------------
    tabPanel(
      tags$span(
        style = "color:white ; font-weight:600 ; font-size:120%",
        "DotPlot"
      ),
      
      fluidRow(
        column(width = 4,
               div(
                 style = "margin-bottom:40px",
                 selectInput("dotplot.show.cat", "Select terms", choices = c(), 
                             size = 20, multiple = T, selectize = F, width = "400px"),
                 verbatimTextOutput("dotp_selected_terms", placeholder = T)
               ),
               actionButton("go.dotplot", "Plot Dotplot", class = "btn-sm", width = "90%",
                            style = "font-size:18px; background-color:midnightblue; font-weight:600; margin-bottom:40px;
                            border-radius:10px; margin-top:20px; border-color:cadetblue")),
        column(width = 8,
               div(
                 style = "overflow-x:auto; width:100%; margin-top:10px",
                 plotOutput("plt6", height = "auto")
               ),
               div(
                 style = "margin-top:30px ; margin-bottom:30px",
                 actionButton("download.go.dplot","Download as pdf", 
                              class = "btn-sm", icon = icon("download"),
                              style = "font-size:20px ; background-color:darkgreen ; 
                                padding:5px 150px ; border-radius:10px")
               ))
      )
    ),
    
    ##### p3: Correlation plot        --------------------------------------------------
    tabPanel(
      tags$span(
        style = "color:white ; font-weight:600 ; font-size:120%",
        "CorrPlot"
      ),
      p("A correlation score between the terms is calculated using Jaccard's similarity 
        index (JC). A hierarchical clustering will be then made using that JC scores, and 
        also for the clustering, several method can be applied !",
        style = "font-weight:600 ; color:darkgreen ; margin-top:20px ; margin-bottom:30px"),
      
      fluidRow(
        column(width = 4,
               div(
                 style = "margin-bottom:20px",
                 selectInput("corrplot.show.cat", "Select terms", choices = c(), 
                             size = 20, multiple = T, selectize = F, width = "400px"),
                 verbatimTextOutput("corrp_selected_terms", placeholder = T)
               ),
               div(
                 style = "margin-bottom:40px",
                 selectInput("corr.methode", "Correlation method",
                             choices = c("circle","square","ellipse","number","shade","color","pie"),
                             selected = "circle"),
                 selectInput("corr.hclust","Select clustering method",
                             choices = c("ward", "ward.D", "ward.D2", "single", "average", "median", "centroid"),
                             selected = "ward.D"),
                 numericInput("tl.cex", "Labels size", value = 0.8, step = 0.1, width = "150px"),
                 numericInput("tl.srt", "Labels rotation", value = 45, step = 5, min = 0, max = 90, width = "150px")
               ),
               actionButton("go.corrplot", "Plot Corrplot", class = "btn-sm", width = "90%",
                            style = "font-size:18px; background-color:midnightblue; font-weight:600; margin-bottom:40px;
                            border-radius:10px; margin-top:20px; border-color:cadetblue")),
        column(width = 8,
               div(
                 style = "overflow-x:auto; margin-top:10px; width:100%",
                 plotOutput("plt.corrplot", height = "auto")
               ))
      )
    ),
    
    ##### p4: Cnet Plot        --------------------------------------------------
    tabPanel(
      tags$span(
        style = "color:white ; font-weight:600 ; font-size:120%",
        "CnetPlot"
      ),
      p("Now, to concider the potentially biological complexities in which a gene may belong 
        to multiple annotation categories and provide information of numeric changes if 
        available, the cnetplot() function will be used to extract complex associations. 
        It depicts the linkages of genes and biological concepts as a network.",
        style = "font-weight:600 ; color:darkgreen ; margin-top:20px ; margin-bottom:30px"),
      
      fluidRow(
        column(width = 4,
               div(
                 style = "margin-bottom:20px",
                 selectInput("cnetplot.show.cat", "Select terms", choices = c(), 
                             size = 15, multiple = T, selectize = F, width = "400px")
               )),
        column(width = 4,
               div(
                 style = "margin-bottom:20px",
                 p("Selected terms", style = "margin-bottom:5px ; font-weight:600"),
                 verbatimTextOutput("cnetp_selected_terms", placeholder = T)
               )),
        column(width = 4,
               div(
                 style = "margin-bottom:40px",
                 selectInput("cnet.layout", "Cnet layout",
                             choices = c("circle","kk","grid","fr"),
                             selected = "kk"),
                 numericInput("cnet.genesize", "Gene labels size", value = 0.8, step = 0.1, width = "150px"),
                 numericInput("cnet.catsize", "Term labels size", value = 1.4, step = 0.1, width = "150px")
               ),
               actionButton("go.cnetplot", "Plot Correlation network", class = "btn-sm", width = "90%",
                            style = "font-size:18px; background-color:midnightblue; font-weight:600; margin-bottom:40px;
                            border-radius:10px; margin-top:20px; border-color:cadetblue")),
        column(width = 10,
               offset = 1,
               div(
                 style = "overflow-x:auto; margin-top:20px; width:100%",
                 plotOutput("plt.cnetplot", height = "auto")
               ),
               div(
                 style = "margin-top:30px ; margin-bottom:30px ; margin-left:100px",
                 actionButton("download.go.cnetplot","Download as pdf", 
                              class = "btn-sm", icon = icon("download"),
                              style = "font-size:20px ; background-color:darkgreen ; 
                                padding:5px 150px ; border-radius:10px")
               ))
      )
    ),
    
    ##### p5: Tree Plot        --------------------------------------------------
    tabPanel(
      tags$span(
        style = "color:white ; font-weight:600 ; font-size:120%",
        "TreePlot"
      ),
      p("A hierarchical clustering between the selected terms below to see the connections !",
        style = "font-weight:600 ; color:darkgreen ; margin-top:20px ; margin-bottom:30px"),
      
      fluidRow(
        column(width = 4,
               div(
                 style = "margin-bottom:20px",
                 selectInput("treeplot.show.cat", "Select terms", choices = c(), 
                             size = 25, multiple = T, selectize = F, width = "400px")
               )),
        column(width = 4,
               div(
                 style = "margin-bottom:20px",
                 p("Selected terms", style = "margin-bottom:5px ; font-weight:600"),
                 verbatimTextOutput("treep_selected_terms", placeholder = T)
               )),
        column(width = 4,
               div(
                 style = "margin-bottom:40px",
                 selectInput("tree.clust.method", "Clustering method",
                             choices = c("ward.D","complete","ward.D2","median","single"),
                             selected = "ward.D"),
                 radioButtons("treeplot_lab_format", "Labels format", choices = c("1","2","3","4"), inline = T),
                 numericInput("treeplot_fontsize", "Font size", value = 5, step = 0.5, width = "150px"),
                 numericInput("treeplot_hexpand", "Hexpand plot", value = 0.15, step = 0.01, width = "150px"),
                 numericInput("treeplot_nCluster", "Split into n clusters", value = 6, step = 1, width = "150px"),
                 numericInput("treeplot_bar_tree", "Bar tree extend", value = 8, step = 0.5, width = "150px"),
                 numericInput("treeplot_tiplab", "Tiplabs", value = 0.6, step = 0.1, width = "150px"),
                 checkboxInput("treeplot_Clust_highlight", "Highlight clusters", value = T, width = "150px")
               ),
               actionButton("go.treeplot", "Plot Treeplot", class = "btn-sm", width = "90%",
                            style = "font-size:18px; background-color:midnightblue; font-weight:600; margin-bottom:40px;
                            border-radius:10px; margin-top:20px; border-color:cadetblue")),
        column(width = 10,
               offset = 1,
               div(
                 style = "overflow-x:auto; margin-top:20px; width:100%",
                 plotOutput("plt.treeplot", height = "auto")
               ),
               div(
                 style = "margin-top:30px ; margin-bottom:30px ; margin-left:100px",
                 actionButton("download.go.treeplot","Download as pdf", 
                              class = "btn-sm", icon = icon("download"),
                              style = "font-size:20px ; background-color:darkgreen ; 
                                padding:5px 150px ; border-radius:10px")
               ))
      )
    ),
    
    ##### p6: Enrichment map Plot        --------------------------------------------------
    tabPanel(
      tags$span(
        style = "color:white ; font-weight:600 ; font-size:120%",
        "EmapPlot"
      ),
      p("Enrichment map organizes enriched terms into a network with edges connecting 
        overlapping gene sets. In this way, mutually overlapping gene sets are tend to 
        cluster together, making it easy to identify functional module.",
        style = "font-weight:600 ; color:darkgreen ; margin-top:20px ; margin-bottom:30px"),
      
      fluidRow(
        column(width = 4,
               div(
                 style = "margin-bottom:20px",
                 selectInput("emapplot.show.cat", "Select terms", choices = c(), 
                             size = 30, multiple = T, selectize = F, width = "400px")
               )),
        column(width = 4,
               div(
                 style = "margin-bottom:20px",
                 p("Selected terms", style = "margin-bottom:5px ; font-weight:600"),
                 verbatimTextOutput("emapp_selected_terms", placeholder = T)
               )),
        column(width = 4,
               div(
                 style = "margin-bottom:40px",
                 selectInput("emap.layout", "Select layout", 
                             choices = c("circle","kk","grid","fr"),
                             selected = "kk"),
                 numericInput("emap.cex.node", "Node size", value = 1.1, step = 0.1),
                 numericInput("emap.cex.label", "Label size", value = 1, step = 0.1),
                 numericInput("emap.cex.line", "edge width", value = .4, step = 0.1),
                 checkboxGroupInput("emap.displayparams", "Parameters to display",
                                    choices = c("Show edges", "Group into clusters", "Show cluster legend")),
                 numericInput("emap.edge.min", "Edge min.similarity", value = 0.4, min = 0, max = 1, step = 0.05),
                 numericInput("emap.clust.n", "Nb of clusters", value = 3, step = 1),
                 numericInput("emap.clust.labs", "Cluster n words", value = 3, step = 1)),
               actionButton("go.emapplot", "Plot Enrichment map", class = "btn-sm", width = "90%",
                            style = "font-size:18px; background-color:midnightblue; font-weight:600; margin-bottom:40px;
                            border-radius:10px; margin-top:20px; border-color:cadetblue"),
               ),
        column(width = 10,
               offset = 1,
               div(
                 style = "overflow-x:auto; margin-top:20px; width:100%",
                 plotOutput("plt.emapplot", height = "auto")
               ),
               div(
                 style = "margin-top:30px ; margin-bottom:30px ; margin-left:100px",
                 actionButton("download.go.emapplot","Download as pdf", 
                              class = "btn-sm", icon = icon("download"),
                              style = "font-size:20px ; background-color:darkgreen ; 
                                padding:5px 150px ; border-radius:10px")
               ))
      )
    )
  ),
  
  #### 1-2- GSEA     -----------------------------------------------------------
  p("1-2 Gene Set Enrichment Analysis",
    style = "color:white ; font-weight:600 ; font-size:150% ; background-color:red ; 
    margin-bottom:10px ; margin-top:60px"),
  
  
  
  ### 2- Grouped Enrichment analysis      ---------------------------------------
  p("2- Grouped enrichment analysis",
    style = "color:darkred ; font-weight:600 ; font-size:160% ; background-color:gold ; 
    margin-bottom:10px ; margin-top:120px"),
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
    message("==========================================")
    message("> DEA is running with the following parameters: ")
    message("==========================================")
    message(sprintf("%-25s: %s", "Selected metadata", input$Sel.metadata))
    message(sprintf("%-25s: %s", "Ident.1", input$ident.1))
    message(sprintf("%-25s: %s", "Ident.2", input$ident.2))
    message(sprintf("%-25s: %s", "Min.log2FC.threshold", input$log2fc))
    message(sprintf("%-25s: %s", "Min.pct.threshold", input$min.pct))
    message(sprintf("%-25s: %s", "Min.diff.pct", input$min.diff.pct))
    message(sprintf("%-25s: %s", "p_val_adj cutoff", "0.05"))
    message("==========================================")
    message("Raw DEA is running on the background with no filtering")
    message("==========================================")
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
  
  DEA_raw <- eventReactive(input$act2, {
    req(seurat())
    seurat <- seurat()
    Idents(seurat) <- as.factor(seurat@meta.data[[input$Sel.metadata]])
    seurat %>% 
      FindMarkers(ident.1 = input$ident.1,
                  ident.2 = input$ident.2,
                  only.pos = FALSE,
                  min.pct = 0,
                  min.diff.pct = 0,
                  logfc.threshold = 0,
                  verbose = F) %>% 
      dplyr::arrange(desc(avg_log2FC)) %>%
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
    message("__________________________________________")
    message("> DE Analysis completed successfully !!!")
    message("==========================================")
    message(sprintf("%-25s: %s", "Total Nb of DEGs", nrow(res)))
    message(sprintf("%-25s: %s", "Nb of Up genes", nrow(res[which(res$avg_log2FC > 0),])))
    message(sprintf("%-25s: %s", "Nb of Down genes", nrow(res[which(res$avg_log2FC < 0),])))
    message("==========================================")
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
  
  observeEvent(DEA_raw(), {
    req(DEA_raw())
    showNotification("Raw DEA results completed from the background !!")
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
  
  # initiate text-info for the heatmap as reactive value:
  hm.text_info <- reactiveVal("Heatmap is not ready yet to be displayed ! 
                              make sure to click on the links above first")
  
  # Print the text:
  output$hm.ready <- renderText({
    hm.text_info()
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
    message(paste0("> Heatmap expression assay successfully created ! (", 
                   nrow(ExpressionAssay()), " genes selected)"))
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
                        annotation_name_gp = list(fontsize = 14,
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
                        annotation_name_gp = list(fontsize = 14,
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
                        annotation_name_gp = list(fontsize = 14,
                                                  col = "navy",
                                                  fontface = "bold"),
                        annotation_legend_param = list(grid_height = unit(0.9,"cm"),
                                                       grid_width = unit(0.7,"cm"),
                                                       labels_gp = gpar(col = "black",
                                                                        fontsize = 12),
                                                       title_gp = gpar(col = "darkred",
                                                                       fontsize = 13,
                                                                       fontface = "bold")))
    }
  })
  
  observeEvent(TopAnnotation(), {
    req(TopAnnotation(), hm.metadata())
    hm.metadata <- hm.metadata()
    message("> Heatmap top-annotation is ready !")
    showNotification("Top-annotation DONE !", type = "message")
    if(ncol(hm.metadata) >= 3){
      warning("NOTE: 3 variables is a maximum possible to annotate with, any more would be ignored !")
    }
    
  })
  
  ##### Right annot       ---------------------------------------------------------
  RightAnnotation <- eventReactive(input$Setup.hm, {
    req(ExpressionAssay())
    exp <- ExpressionAssay()
    dea.res <- DEA_results()
    rowAnnotation(
      "log2fc" = anno_barplot(dea.res$avg_log2FC[which(dea.res$Genes %in% rownames(exp))],
                              axis = T, border = T, cex = 1, bar_width = 1),
      width = unit(1.5, units = "cm")
    )
  })
  
  observeEvent(RightAnnotation(), {
    req(RightAnnotation())
    message("> Heatmap right-annotation is ready !")
    showNotification("Right-annotation DONE !", type = "message")
  })
  
  # Update text-info:
  observeEvent(input$Setup.hm, {
    req(RightAnnotation(), TopAnnotation())
    hm.text_info("Heatmap is now ready to be plotted ! Click the bottom below")
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
                            column_title = paste0("Expression of ",
                                                  nrow(exp), " genes across ", 
                                                  ncol(meta), " condition(s)"),
                            heatmap_legend_param = list(legend_height = unit(3,"cm"),
                                                        direction = "horizontal",
                                                        legend_width = unit(7,"cm"),
                                                        title = "log_expression",
                                                        title_position = "topcenter",
                                                        title_gp = gpar(fontsize=13, col="darkred", fontface="bold"),
                                                        label_gp = gpar(fontsize=10, col="black"),
                                                        legend_gp = gpar(fontsize=10)),
                            column_title_gp = gpar(col = "darkred",
                                                   fontsize = 17,
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
  }, width = 900, height = 800)
  
  # Saving parameters:
  observeEvent(input$download.hmplot, {
    showModal(
      modalDialog(
        title = "Save Heatmap as PDF",
        numericInput("width_hmplot", "Width (in inches):", value = 6, min = 4, step = 0.5),
        numericInput("height_hmplot", "Height (in inches):", value = 8, min = 4, step = 0.5),
        textInput("hmplot_name", "Filename:", value = ""),
        footer = tagList(
          modalButton("Cancel"),
          downloadButton("download.hmplot.modal", "Download", class = "btn-success")
        )
      )
    )
  })
  
  
  # IV- Markers inspections      ---------------------------------------------------
  
  ## 1- Displaying the genes    ----------------------------------------------------
  output$genes1 <- renderPrint({
    req(input$show_upgenes)
    Up_Genes <- Up_Genes()
    if(input$show_upgenes == "All"){
      cat(Up_Genes, sep = " | ")
      
    } else if(input$show_upgenes == "Top10"){
      cat(head(Up_Genes, 10), sep = " | ")
      
    } else if(input$show_upgenes == "Top50"){
      cat(head(Up_Genes, 50), sep = " | ")
    }  
  })
  
  output$genes2 <- renderPrint({
    req(input$show_downgenes)
    Down_Genes <- Down_Genes()
    if(input$show_downgenes == "All"){
      cat(Down_Genes, sep = " | ")
      
    } else if(input$show_downgenes == "Top10"){
      cat(head(Down_Genes, 10), sep = " | ")
      
    } else if(input$show_downgenes == "Top50"){
      cat(head(Down_Genes, 50), sep = " | ")
    }  
  })
  
  ## 2- Fill the cumulative list of markers     -----------------------------------
  
  # Initiate the cumulative list:
  List_markers <- reactiveVal(list())
  
  # Fill the list:
  observeEvent(input$add.markers, {
    req(DEA_results(), DEA_raw())
    
    res <- DEA_results()
    fc <- res[["avg_log2FC"]] %>% setNames(res[["Genes"]])
    
    tmp_list <- list(
      Up = Up_Genes(),
      Down = Down_Genes(),
      All = c(Up_Genes(), Down_Genes()),
      All_FC = fc,
      DEA_res = res,
      all_DEA = DEA_raw()
    )
    
    current_list <- List_markers()
    current_list[[input$name.list]] <- tmp_list
    List_markers(current_list)
    
    showNotification(paste0(input$name.list, " added to the list"), type = "message", duration = 5)
  })
  
  output$list_markers <- renderPrint({
    str(List_markers())
  })
  
  # Clear the list:
  observeEvent(input$ac, {
    if(length(List_markers()) == 0){
      showNotification("The list is already empty.", type = "error")
    } else {
      showModal(
        modalDialog(
          title = "Confirm Action",
          "Are you sure you want to clear the list of markers?",
          footer = tagList(
            modalButton("Cancel"),
            actionButton("confirm_clear", "Yes, clear it!", class = "btn-danger")
          )
        )
      )
    }
  })
  
  observeEvent(input$confirm_clear, {
    removeModal()
    List_markers(list())
    showNotification("The list of markers has been cleared.", type = "warning")
  })
  
  ## 3- VennDiagram plot     -----------------------------------------------
  
  # Update selectinput:
  observeEvent(List_markers(), {
    req(List_markers())
    List_markers <- List_markers()
    updateSelectInput(session, "sel.markers.list", choices = names(List_markers))
  })
  
  # Venn Object:
  Venn_obj <- eventReactive(input$display_venn, {
    req(List_markers())
    List_markers <- List_markers()
    genes <- switch(input$up_down,
                    "Up genes" = "Up",
                    "Down genes" = "Down",
                    "All genes" = "All")
    sets <- list()
    for(set in input$sel.markers.list){sets[[set]] <- List_markers[[set]][[genes]]}
    if(length(sets) <= 4){
      Display_Venn(sets, 
                   set.names = input$sel.markers.list, 
                   colpalette = MyPalette,
                   Padding = input$vplot.padding,
                   text.size = input$vplot.text.size,
                   set.name.size = input$vplot.setname.size)
    } else {
      stop("Please select only up to 4 set of genes to compare !")
    }
    
  })
  
  # Venn Plot:
  Venn_plot <- eventReactive(Venn_obj(), {
    req(Venn_obj())
    Venn_obj <- Venn_obj()
    if(length(Venn_obj) <= 4 & length(Venn_obj) >= 2){
      Venn_obj[["plot"]]
    } else {
      stop("Please check out your gene set selection ! can't handle over 4 selections")
    }
  })
  
  output$vennplot <- renderPlot({
    req(Venn_plot())
    Venn_plot()
  }, height = 600, width = 700)
  
  # Saving parameters:
  observeEvent(input$download.vplot, {
    showModal(
      modalDialog(
        title = "Save VennPlot as PDF",
        numericInput("width_vplot", "Width (in inches):", value = 6, min = 4, step = 0.5),
        numericInput("height_vplot", "Height (in inches):", value = 8, min = 4, step = 0.5),
        textInput("vplot_name", "Filename:", value = ""),
        footer = tagList(
          modalButton("Cancel"),
          downloadButton("download.vplot.modal", "Download", class = "btn-success")
        )
      )
    )
  })
  
  ## 4- Common & group-specific genes     ----------------------------------------
  
  # Update selection inputs:
  observeEvent(Venn_obj(), {
    req(Venn_obj())
    Venn_obj <- Venn_obj()
    
    updateSelectInput(session, "sel.intersection", 
                      choices = names(Venn_obj[["intersections"]]), selected = "")
    updateSelectInput(session, "specific.genes", 
                      choices = names(Venn_obj[["group_specific"]]), selected = "")
  })
  
  observeEvent(input$specific.genes, {
    req(Venn_obj())
    Venn_obj <- Venn_obj()
    
    updateSelectInput(session, "specific.genes2", 
                      choices = names(Venn_obj[["group_specific"]][[input$specific.genes]]), 
                      selected = "")
  })
  
  # Display common genes:
  output$nb_common_genes <- renderPrint({
    req(Venn_obj())
    Venn_obj <- Venn_obj()
    print(length(Venn_obj[["intersections"]][[input$sel.intersection]]))
  })
  
  output$common_genes <- renderPrint({
    req(Venn_obj())
    Venn_obj <- Venn_obj()
    print(Venn_obj[["intersections"]][[input$sel.intersection]])
  })
  
  # Display group-specific genes:
  output$nb_genes <- renderPrint({
    req(Venn_obj())
    Venn_obj <- Venn_obj()
    print(length(Venn_obj[["group_specific"]][[input$specific.genes]][[input$specific.genes2]]))
  })
  
  output$specific_genes <- renderPrint({
    req(Venn_obj())
    Venn_obj <- Venn_obj()
    print(Venn_obj[["group_specific"]][[input$specific.genes]][[input$specific.genes2]])
  })
  
  
  # V- PEA      ---------------------------------------------------------------
  
  # Initiate a list to store GO objects:
  List_GO <- reactiveVal(list())
  
  # Upload list of saved markers:
  observeEvent(List_markers(), {
    req(List_markers())
    List_markers <- List_markers()
    updateSelectInput(session, "go.list.genes", choices = names(List_markers))
  })
  
  ## 1- Individual enrichment     ---------------------------------------------
  
  ### 1-1- GO analysis          --------------------------------------------
  
  # Run GO analysis:
  GOobject <- eventReactive(input$run.go, {
    req(List_markers())
    List_markers <- List_markers()
    MyGenes <- List_markers[[input$go.list.genes]][[input$go.list.genes.reg]]
    enrichGO(gene = MyGenes,
             OrgDb = org.Hs.eg.db,
             keyType = input$go.keytype,
             ont = input$go.ont,
             qvalueCutoff = 0.1,
             pvalueCutoff = 0.05,
             readable = T,
             maxGSSize = input$go.maxgeneset,
             minGSSize = 10)
  })
  
  observeEvent(input$run.go, {
    showNotification("GO is running ...", duration = 10, type = "message")
    message(paste0("> Running GO analysis on ... ", input$go.list.genes, "_", input$go.list.genes.reg))
  })
  
  observeEvent(GOobject(), {
    req(GOobject())
    GOobject <- GOobject()
    showNotification("GO analysis completed !", duration = 5, type = "message")
    message("> GO analysis completed successfully !!!")
    message("==========================================")
    message(sprintf("%-25s: %s", "Total Nb of terms found", nrow(GOobject@result)))
    
  })
  
  # Add the GO object to the list:
  observeEvent(input$addgo, {
    req(GOobject())
    current_list <- List_GO()
    contrast_name <- paste0("GOobj_",input$go.list.genes,"_",input$go.list.genes.reg)
    if(contrast_name %in% names(current_list)){
      current_list[[contrast_name]] <- GOobject()
      List_GO(current_list)
      showNotification(paste0("The contrast ", contrast_name, " already exists in the list!"), 
                       type = "warning")
    } else {
      current_list[[contrast_name]] <- GOobject()
      List_GO(current_list)
      showNotification(paste0("GO analysis for ", contrast_name, " added to the list!"), type = "message")
    }
  })
  
  output$showGOobj <- renderPrint({
    req(List_GO())
    List_GO <- names(List_GO())
    cat(List_GO, sep = "\n")
  })
  
  # Clear the list:
  observeEvent(input$clear.go.list, {
    if(length(List_GO()) == 0){
      showNotification("The list is already empty.", type = "error")
    } else {
      showModal(
        modalDialog(
          title = "Confirm Action",
          "Are you sure you want to clear this list?",
          footer = tagList(
            modalButton("Cancel"),
            actionButton("confirm_clear.go", "Yes, clear it!", class = "btn-danger")
          )
        )
      )
    }
  })
  
  observeEvent(input$confirm_clear.go, {
    removeModal()
    List_GO(list())
    showNotification("The list of GO objects is now empty.", type = "warning")
  })
  
  #### p1: Result table     ---------------------------------------------------------
  go.tab <- eventReactive(GOobject(), {
    req(GOobject())
    GOobject <- GOobject()
    GOobject@result %>% 
      dplyr::filter(p.adjust < 0.05) %>% 
      dplyr::mutate(RichFactor = round(Count / as.numeric(sub("/\\d+","",BgRatio)),5)) %>% 
      dplyr::select(ID, Description, RichFactor, p.adjust, GeneRatio, BgRatio, Count, geneID)
  })
  
  observeEvent(go.tab(), {
    req(go.tab())
    go.tab <- go.tab()
    message(sprintf("%-25s: %s", "Nb of significant term", nrow(go.tab)))
    
    # Updating terms lists:
    updateSelectInput(session, "dotplot.show.cat", choices = go.tab[["Description"]])
    updateSelectInput(session, "cnetplot.show.cat", choices = go.tab[["Description"]])
    updateSelectInput(session, "corrplot.show.cat", choices = go.tab[["Description"]])
    updateSelectInput(session, "treeplot.show.cat", choices = go.tab[["Description"]])
    updateSelectInput(session, "emapplot.show.cat", choices = go.tab[["Description"]])
  })
  
  output$tab3 <- renderDataTable({
    req(go.tab())
    go.tab() %>% 
      dplyr::select(ID, Description, GeneRatio, BgRatio, RichFactor, p.adjust) %>% 
      datatable(options = list(pageLength = 5, scrollX = TRUE))
  })
  
  #### p2: Dot plot   ----------------------------------------------------------
  
  # Selected terms:
  output$dotp_selected_terms <- renderPrint({
    selected_terms <- input$dotplot.show.cat
    cat(selected_terms, sep = "\n")
  })
  
  # Plot Dot plot:
  go.dplot <- eventReactive(input$go.dotplot, {
    req(go.tab())
    go.tab() %>% 
      dplyr::filter(Description %in% input$dotplot.show.cat) %>% 
      dplyr::mutate(Description = ifelse(nchar(Description) <= 50,
                                         Description,
                                         paste0(substr(Description,1,46), "...."))) %>% 
      ggplot(aes(x= RichFactor, y= fct_reorder(Description, RichFactor)))+
      geom_segment(aes(xend= 0, yend= Description))+
      geom_point(aes(color= p.adjust, size= Count))+
      scale_color_viridis_c(guide = guide_colorbar(reverse = T))+
      scale_size_continuous(range = c(3,10))+
      theme_linedraw()+
      theme(panel.grid = element_blank(),
            panel.border = element_blank(),
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5, colour = "darkred", 
                                      margin = margin(b=0.2, unit = "in")),
            plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"),
            axis.title.x = element_text(size = 16, face = "bold", colour = "darkred", 
                                        margin = margin(t=0.2, unit = "in")),
            axis.title.y = element_blank(),
            axis.text = element_text(size = 14, face = "bold"),
            legend.title = element_text(size = 15, face = "bold", colour = "darkred",
                                        margin = margin(b=0.2, unit = "in")),
            legend.text = element_text(size = 13),
            legend.box.margin = margin(l=0.2, unit = "in"))+
      labs(title = "GO Enriched Terms")
  })
  
  output$plt6 <- renderPlot({
    req(go.dplot())
    go.dplot()
  }, width = 800, height = 700)
  
  # Saving parameters:
  observeEvent(input$download.go.dplot, {
    showModal(
      modalDialog(
        title = "Save Dotplot as PDF",
        numericInput("width_godotplot", "Width (in inches):", value = 6, min = 4, step = 0.5),
        numericInput("height_godotplot", "Height (in inches):", value = 8, min = 4, step = 0.5),
        textInput("go.dplot_name", "Filename:", value = ""),
        footer = tagList(
          modalButton("Cancel"),
          downloadButton("download.go.dplot.modal", "Download", class = "btn-success")
        )
      )
    )
  })
  
  # Apply pairwise term-sim to GOobj:
  GOobject_paired <- eventReactive(GOobject(), {
    req(GOobject())
    GOobject <- GOobject()
    pairwise_termsim(GOobject, method = "JC", showCategory = 250)
  })
  
  #### p3: Correlation plot   ----------------------------------------------------------
  corrplot_sel_terms <- reactiveVal(character(0))
  
  observeEvent(input$corrplot.show.cat, {
    selected_terms <- input$corrplot.show.cat
    corrplot_sel_terms(selected_terms)
  })
  
  # Selected terms:
  output$corrp_selected_terms <- renderPrint({
    selected_terms <- input$corrplot.show.cat
    cat(selected_terms, sep = "\n")
  })
  
  Corr_plot <- eventReactive(input$go.corrplot, {
    req(GOobject_paired(), corrplot_sel_terms())
    GOobject_paired <- GOobject_paired()
    term_sim_mat <- GOobject_paired@termsim[corrplot_sel_terms(),corrplot_sel_terms()]
    colnames(term_sim_mat) <- ifelse(nchar(colnames(term_sim_mat)) <= 25, 
                                     colnames(term_sim_mat), 
                                     paste0(substr(colnames(term_sim_mat),1,23),"...."))
    rownames(term_sim_mat) <- ifelse(nchar(rownames(term_sim_mat)) <= 60, 
                                     rownames(term_sim_mat), 
                                     paste0(substr(rownames(term_sim_mat),1,58),"...."))
    
    corrplot::corrplot(term_sim_mat,
                       type = "upper",
                       method = input$corr.methode,
                       tl.col = "black",
                       hclust.method = input$corr.hclust,
                       order = "hclust",
                       tl.cex = input$tl.cex,
                       tl.srt = input$tl.srt,
                       cl.cex = 1.1,
                       cl.ratio = 0.16,
                       mar = c(0,0,0,0))
  })
  
  output$plt.corrplot <- renderPlot({
    req(Corr_plot())
    Corr_plot()
  }, height = 1000,  width = 1100)
  
  #### p4: Cnet plot   ----------------------------------------------------------
  cnetplot_sel_terms <- reactiveVal(character(0))
  
  observeEvent(input$cnetplot.show.cat, {
    selected_terms <- input$cnetplot.show.cat
    cnetplot_sel_terms(selected_terms)
  })
  
  # Selected terms:
  output$cnetp_selected_terms <- renderPrint({
    selected_terms <- input$cnetplot.show.cat
    cat(selected_terms, sep = "\n")
  })
  
  # Construct the cnet plot:
  CNetPlot <- eventReactive(input$go.cnetplot, {
    req(GOobject(),  List_markers())
    GOobject <- GOobject()
    selected_terms <- cnetplot_sel_terms()
    List_markers <- List_markers()
    fc <- List_markers[[input$go.list.genes]][["All_FC"]]
    gene_lbl_size <- input$cnet.genesize
    cat_lbl_size <- input$cnet.catsize
    selected_layout <- input$cnet.layout
    
    col_gradient <- switch(input$go.list.genes.reg,
                           "All" = scale_color_gradientn(
                             colours = c("midnightblue", "white", "darkred"),
                             values = c(0, 0.5, 1),
                             limits = c(-5, 5)
                           ),
                           "Up" = scale_color_gradientn(
                             colours = c("#ffffff", "#ffcc22", "#991111", "#500000"),
                             values = c(0, 0.3, 0.7, 1),
                             limits = c(0, 5)
                           ),
                           "Down" = scale_color_gradientn(
                             colours = c("#000050", "#111199", "#22ccff", "#ffffff"),
                             values = c(0, 0.3, 0.7, 1),
                             limits = c(-5, 0)
                           ))
    
    cnetplot(x = GOobject,
             showCategory = selected_terms,
             layout = selected_layout,
             cex.params = list(gene_node = 0.7, 
                               gene_label = gene_lbl_size, 
                               category_node = 1.5,
                               category_label = cat_lbl_size),
             color.params = list(category = "#2277cc",
                                 gene = "#552299",
                                 edge = T,
                                 foldChange = fc))+
      labs(color = "logFC")+
      theme(legend.box.margin = margin(l=0.3, unit = "in"),
            legend.title = element_text(size = 14, face = "bold", colour = "darkred", 
                                        margin = margin(t=0.3,b=0.1, unit = "in")),
            legend.text = element_text(size = 10, face = "bold"))+
      col_gradient
  })
  
  # Plot cnet:
  output$plt.cnetplot <- renderPlot({
    req(CNetPlot())
    CNetPlot()
  }, height = 1000, width = 1200)
  
  # Saving parameters:
  observeEvent(input$download.go.cnetplot, {
    showModal(
      modalDialog(
        title = "Save Cnetplot as PDF",
        numericInput("width_gocnettplot", "Width (in inches):", value = 10, min = 4, step = 0.5),
        numericInput("height_gocnettplot", "Height (in inches):", value = 10, min = 4, step = 0.5),
        textInput("go.cnetplot_name", "Filename:", value = ""),
        footer = tagList(
          modalButton("Cancel"),
          downloadButton("download.go.cnetplot.modal", "Download", class = "btn-success")
        )
      )
    )
  })
  
  #### p5: Tree plot   ----------------------------------------------------------
  treeplot_sel_terms <- reactiveVal(character(0))
  
  observeEvent(input$treeplot.show.cat, {
    selected_terms <- input$treeplot.show.cat
    treeplot_sel_terms(selected_terms)
  })
  
  # Selected terms:
  output$treep_selected_terms <- renderPrint({
    selected_terms <- input$treeplot.show.cat
    cat(selected_terms, sep = "\n")
  })
  
  # Construct the treeplot:
  Treeplot <- eventReactive(input$go.treeplot, {
    req(GOobject_paired(), List_markers())
    GOobject_paired <- GOobject_paired()
    selected_terms <- treeplot_sel_terms()
    lab_format <- switch(input$treeplot_lab_format,
                         "1" = 10,
                         "2" = 20,
                         "3" = 30,
                         "4" = 40)
    treeplot(GOobject_paired, 
             showCategory = selected_terms,
             fontsize = input$treeplot_fontsize,
             color = "p.adjust",
             cex_category = 1,
             cluster.params = list(n= input$treeplot_nCluster,
                                   method = input$tree.clust.method,
                                   color = MyPalette[1:input$treeplot_nCluster],
                                   label_words_n = 3,
                                   label_format = input$treeplot_lab_format),
             offset.params = list(bar_tree = input$treeplot_bar_tree, 
                                  tiplab = input$treeplot_tiplab, 
                                  extend = 0.5, 
                                  hexpand = input$treeplot_hexpand),
             highlight.params = list(align = "both"), hilight = input$treeplot_Clust_highlight)+
      scale_color_gradientn(colours = c("darkred","gold"),
                            values = c(0,1),
                            limits = c(0,0.05),
                            guide = guide_colorbar(title = "pValue.adj"))+
      scale_size_continuous(range = c(1,7))+
      labs(colour = "p.adjusted")+
      theme(legend.title = element_text(face = "bold", 
                                        colour = "darkred",
                                        size = 14,
                                        margin = margin(b=0.3, unit = "in")))
  })
  
  # Plot treeplot:
  output$plt.treeplot <- renderPlot({
    req(Treeplot())
    Treeplot()
  }, height = 900, width = 900)
  
  # Saving parameters:
  observeEvent(input$download.go.treeplot, {
    showModal(
      modalDialog(
        title = "Save treeplot as PDF",
        numericInput("width_gotreeplot", "Width (in inches):", value = 10, min = 4, step = 0.5),
        numericInput("height_gotreeplot", "Height (in inches):", value = 10, min = 4, step = 0.5),
        textInput("go.treeplot_name", "Filename:", value = ""),
        footer = tagList(
          modalButton("Cancel"),
          downloadButton("download.go.treeplot.modal", "Download", class = "btn-success")
        )
      )
    )
  })
  
  #### p6: Emap plot   ----------------------------------------------------------
  emapplot_sel_terms <- reactiveVal(character(0))
  
  observeEvent(input$emapplot.show.cat, {
    selected_terms <- input$emapplot.show.cat
    emapplot_sel_terms(selected_terms)
  })
  
  # Selected terms:
  output$emapp_selected_terms <- renderPrint({
    selected_terms <- input$emapplot.show.cat
    cat(selected_terms, sep = "\n")
  })
  
  # Construct the emapplot:
  emapplot <- eventReactive(input$go.emapplot, {
    req(GOobject_paired(), List_markers())
    GOobject_paired <- GOobject_paired()
    selected_terms <- emapplot_sel_terms()
    display_params <- input$emap.displayparams
    show_edges <- "Show edges" %in% display_params
    group_clusters <- "Group into clusters" %in% display_params
    show_legend <- "Show cluster legend" %in% display_params
    
    enrichplot::emapplot(GOobject_paired,
                         showCategory = selected_terms,
                         repel = T,
                         edge.params = list(show = show_edges, min = input$emap.edge.min),
                         cex.params = list(category_node = input$emap.cex.node, 
                                           category_label = input$emap.cex.label, 
                                           line = input$emap.cex.line,
                                           label_group = 1),
                         cluster.params = list(cluster = group_clusters, 
                                               method = stats::kmeans, 
                                               n = input$emap.clust.n, 
                                               legend = show_legend, 
                                               label_style = "shadowtext", 
                                               label_words_n = input$emap.clust.labs, 
                                               label_format = 30),
                         layout.params = list(layout = input$emap.layout),
                         node_label = "category")+
      theme(legend.box.margin = margin(l=0.2, unit = "in"),
            plot.margin = margin(t= 0.1, b=0.1, unit = "in"),
            legend.title = element_text(face = "bold", 
                                        colour = "darkred",
                                        size = 14,
                                        margin = margin(b=0.3, unit = "in")))+
      scale_fill_gradientn(colours = c("darkred","gold"),
                           values = c(0,1),
                           limits = c(0, 0.05))+
      labs(fill= "p.adjust")
  })
  
  # Plot emapplot:
  output$plt.emapplot <- renderPlot({
    req(emapplot())
    emapplot()
  }, height = 900, width = 900)
  
  # Saving parameters:
  observeEvent(input$download.go.emapplot, {
    showModal(
      modalDialog(
        title = "Save Emapplot as PDF",
        numericInput("width_goemapplot", "Width (in inches):", value = 10, min = 4, step = 0.5),
        numericInput("height_goemapplot", "Height (in inches):", value = 10, min = 4, step = 0.5),
        textInput("go.emapplot_name", "Filename:", value = ""),
        footer = tagList(
          modalButton("Cancel"),
          downloadButton("download.go.emapplot.modal", "Download", class = "btn-success")
        )
      )
    )
  })
  
  
  
  
  
  
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
      print(Dplot())
      dev.off()
    }
  )
  
  # Heatmap plot:
  output$download.hmplot.modal <- downloadHandler(
    filename = function(){paste0(input$hmplot_name,".pdf")},
    content = function(file){
      pdf(file, width = input$width_hmplot, height = input$height_hmplot)
      hm <- hm()
      
      ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom")
      dev.off()
    }
  )
  
  # List of markers:
  output$download.markers <- downloadHandler(
    filename = function() {paste0("Markers-",Sys.Date(),".rds")},
    content = function(file) {saveRDS(List_markers(), file)}
  )
  
  # VennDiagram plot:
  output$download.vplot.modal <- downloadHandler(
    filename = function(){paste0(input$vplot_name,".pdf")},
    content = function(file){
      pdf(file, width = input$width_vplot, height = input$height_vplot)
      print(Venn_plot())
      dev.off()
    }
  )
  
  # GO data table:
  output$download.tab.go <- downloadHandler(
    filename = function(){paste0("GO_res_", input$go.list.genes,"_", input$go.list.genes.reg,".xlsx")},
    content = function(file){write_xlsx(go.tab(), path = file)}
  )
  
  # List of GO objects:
  output$dwnload.GoObjList <- downloadHandler(
    filename = function() {paste0("GO_objects-",Sys.Date(),".rds")},
    content = function(file) {saveRDS(List_GO(), file)}
  )
  
  # GO Dot plot:
  output$download.go.dplot.modal <- downloadHandler(
    filename = function(){paste0(input$go.dplot_name,".pdf")},
    content = function(file){
      pdf(file, width = input$width_godotplot, height = input$height_godotplot)
      print(go.dplot())
      dev.off()
    }
  )
  
  # GO Cnetwork plot:
  output$download.go.cnetplot.modal <- downloadHandler(
    filename = function(){paste0(input$go.cnetplot_name,".pdf")},
    content = function(file){
      pdf(file, width = input$width_gocnettplot, height = input$height_gocnettplot)
      print(CNetPlot())
      dev.off()
    }
  )
  
  # GO Tree plot:
  output$download.go.treeplot.modal <- downloadHandler(
    filename = function(){paste0(input$go.treeplot_name,".pdf")},
    content = function(file){
      pdf(file, width = input$width_gotreeplot, height = input$height_gotreeplot)
      print(Treeplot())
      dev.off()
    }
  )
  
  # GO Emap plot:
  output$download.go.emapplot.modal <- downloadHandler(
    filename = function(){paste0(input$go.emapplot_name,".pdf")},
    content = function(file){
      pdf(file, width = input$width_goemapplot, height = input$height_goemapplot)
      print(Emapplot())
      dev.off()
    }
  )
  
  
}




#  ShinyApp run     ########################
shinyApp(ui, server)














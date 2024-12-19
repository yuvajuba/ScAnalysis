#   Setting     #################

# library(knitr)
library(Seurat)
library(dplyr)
# library(tidyr)
library(ggplot2)
library(tibble)
library(writexl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
# library(fgsea)
library(enrichplot)
# library(msigdbr)
library(ComplexHeatmap)
library(DT)
library(stringr)
library(shiny)
library(bslib)
# library(readxl)
# library(Polychrome)
library(circlize)
# library(RColorBrewer)
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


















#   Shiny Application   ##################

options(shiny.maxRequestSize = 900*1024^2)

#  User interface    ####################


ui <- fluidPage(
  
  theme = bs_theme(bootswatch = "flatly"),
  h2("1- Importing the data", 
     style = "color:gold ; font-weight:700 ; background-color:black ; text-align:left ; margin-top:20px"),
  
  ##  I- Import data   --------------------------------------------------------
  
  fluidRow(
    column(
      width = 6,
      fileInput(inputId = "seurat",
                label = p("Load the seurat object", 
                          style = "color:darkred ; font-weight:600 ; font-size:120%"),
                accept = ".rds")
    ),
    column(
      width = 5,
      offset = 1,
      verbatimTextOutput("text")
    )
  ),
  
  h2("2- Clustering selection",
     style = "color:gold ; font-weight:700 ; background-color:black ; text-align:left ; margin-top:120px ; margin-bottom:20px"),
  
  ##  II- Cluster Selection   -------------------------------------------------
  
  fluidRow(
    column(
      width = 6,
      selectInput(inputId = "clusters",
                  label = p("Select your metadata",
                            style = "color:darkred ; font-weight:600 ; font-size:120%"),
                  choices = c()),
      p("From the metadata, select the column you want to label with your plot 
        (for example, clusters at different resolutions, your global conditions ...)",
        style = "font-style:italic ; color:darkgray ; font-size:80%"),
      plotOutput("plt1")
    ),
    column(
      width = 5,
      offset = 1,
      actionButton(inputId = "act",
                   label = "View the summary",
                   class = "btn-lg btn-success"),
      p("", style = "margin-top:20px"),
      tableOutput("tab1")
    )
  ),
  
  h2("3- Differential expressed analysis",
     style = "color:gold ; font-weight:700 ; background-color:black ; text-align:left ; margin-top:120px ; margin-bottom:20px"),
  
  ##  III- DEA   -------------------------------------------------------------------
  
  ###  1- Initiate   ---------------------------------------------------------
  
  fluidRow(
    column(
      width = 5,
      p("From the column selected above, you'll select the 2 clusters from which you want to find 
        the differential expressed genes. It would be cluster1 vs cluster2",
        style = "font-size:120% ; font-weight:500 ; color:navy"),
      
      selectInput(inputId = "ident.1",
                  label = "select cluster 1",
                  choices = c()),
      selectInput(inputId = "ident.2",
                  label = "select cluster 2",
                  choices = c()),
      
      p("Then we'll set up some parameters necessary to perform the DEA from the 
        FindMarkers() function",
        style = "font-size:120% ; font-weight:500 ; color:navy"),
      
      sliderInput(inputId = "log2fc",
                  label = "logfc.threshold",
                  min = 0,
                  max = 3,
                  step = 0.1,
                  value = 0.5),
      sliderInput(inputId = "min.pct",
                  label = "min percentage of cell expressing the gene",
                  min = 0,
                  max = 1,
                  step = 0.05,
                  value = 0),
      checkboxInput(inputId = "only.pos",
                    label = "Select only Up genes",
                    value = F)
    ),
    column(
      width = 6,
      offset = 1,
      textOutput(outputId = "txt2"),
      plotOutput(outputId = "plt2")
    )
  ),
  
  ###  2- Additional filtering   ---------------------------------------------
  
  fluidRow(
    p("Additional filtering",
      style = "color:darkred ; font-weight:600 ; font-size:120% ; background-color:gold ; margin-top:10px"),
    
    column(
      width = 3,
      sliderInput("pct.1_filt", "Filter pct.1", min = 0, max = 1, value = c(0,1)),
      p("Define the minimum and maximum percentage of cells expressing the gene in condition 1. 
        This allows you to filter genes based on their expression prevalence within the selected group.",
        style = "font-style:italic ; color:darkgray ; font-size:80% ; margin-bottom:40px"),
      actionButton(inputId = "act2",
                   label = "Run DEA",
                   class = "btn-lg btn-primary"),
      textOutput("dea_log")
    ),
    
    column(
      width = 3,
      sliderInput("pct.2_filt", "Filter pct.2", min = 0, max = 1, value = c(0,1)),
      p("Define the minimum and maximum percentage of cells expressing the gene in condition 2.",
        style = "font-style:italic ; color:darkgray ; font-size:80%")
    ),
    
    column(
      width = 3,
      sliderInput("diffpct_filt", "Filter Diff_pct", min = 0, max = 1, value = 0),
      p("Define the threshold for the absolute difference in gene expression percentages between the 
        two conditions. This helps to identify genes with a substantial expression difference 
        (|pct.1 - pct.2|).",
        style = "font-style:italic ; color:darkgray ; font-size:80%")
    ),
    
    column(
      width = 3,
      sliderInput("logfc_filt", "Up range avg_log2FC", min = 0, max = 10, value = c(0,10), step = 0.1),
      sliderInput("logfc_filt2", "Down range avg_log2FC", min = -10, max = 0, value = c(-10,0), step = 0.1),
      p("Define the ranges for filtering genes based on their log2 fold change. 
        Use the 'Range Up avg_log2FC' slider to select positively regulated genes, 
        and the 'Range Down avg_log2FC' slider for negatively regulated genes.",
        style = "font-style:italic ; color:darkgray ; font-size:80%")
    )
  ),
  
  
  ###  3- Results   ----------------------------------------------------------
  
  navbarPage(
    "DEA results",
    
    #### p1: Datatable   ------------------------------------------------------
    tabPanel("Table",
             fluidRow(
               dataTableOutput("tab2"),
               p("Saving parameters", 
                 style = "font-size:120% ; font-weight:600 ; background-color:gold ; color:darkred"),
               textInput("filename", "File name (.xlsx)", value = "DEA_results"),
               downloadButton("download","Download...", class = "btn-danger")
             )),
    
    #### p2: Gene lists    ----------------------------------------------------
    
    ##### a: Displaying lists     ----------------------------------------------
    tabPanel("Genes",
             p("Displaying the list of the differential expressed genes with 
               the parameters above",
               style = "color:darkgreen ; font-weight:600 ; font-size:130%"),
             p("Up regulated genes",
               style = "color:darkred ; font-weight:600 ; font-size:120% ; background-color:gold"),
             verbatimTextOutput("genes1"),
             p("Down regulated genes",
               style = "color:darkred ; font-weight:600 ; font-size:120% ; background-color:gold"),
             verbatimTextOutput("genes2"),
             
             ##### b: Cumulative list of markers     ----------------------------
             h3("Cumulative marker list for comparison",
                style = "color:white ; font-weight:600 ; background-color:red ; margin-top:30px"),
             p("In this section, you can build a cumulative list of markers from multiple differential 
               expression analyses. This allows you to compare the markers identified across different 
               conditions or contrasts.",
               style = "color:darkgreen ; font-weight:600 ; font-size:120%"),
             p("For each analysis, enter the name of the contrast (e.g., 'Cluster1_vs_Cluster5'), 
               then click the 'Add to main list' button. The list will automatically update with the 
               Up and Down-regulated genes from the current comparison.",
               style = "color:navy ; font-weight:600 ; font-size:120%"),
             p("NOTE: Do not close the app between analyses to preserve the list.",
               style = "color:red ; font-weight:600 ; font-size:100%"),
             
             fluidRow(
               column(
                 width = 4,
                 textInput("txt4", "name your contrast"),
                 actionButton("act8","Add to main list", class = "btn-lg btn-success")
               ),
               column(
                 width = 7,
                 offset = 1,
                 verbatimTextOutput("showMarkers"),
                 downloadButton("download_all_markers", "Download the whole list", class = "btn-lg btn-primary"),
                 p("", style = "margin-top:30px"),
                 actionButton("ac","Clear you list", class = "btn-lg btn-danger")
               ),
               
               ##### c: Comparion between contrasts     ----------------------------
               h3("Compare the markers",
                  style = "color:white ; font-weight:600 ; background-color:red ; margin-top:30px"),
               p("Finding common or specific genes among the contrasts !",
                 style = "color:darkgreen ; font-weight:600 ; font-size:120%"),
               
               column(
                 width = 3,
                 selectInput("contrasts", "Your selection", 
                             choices = c(), 
                             multiple = T, selectize = F, width = "450px", size = 6),
                 p("You can make up to 4 comparison, more will throw you an error !",
                   style = "color:darkgray ; font-weight:600 ; font-size:80% ; font-style:italic"),
                 
                 radioButtons("up_down", "Which genes to compare ?", 
                              choices = c("Up genes","Down genes","All genes"), inline = F, selected = ""),
                 
                 actionButton("act9", "Display Venn diagram", class = "btn-sm btn-primary")
               ),
               
               column(
                 width = 8,
                 offset = 1,
                 
                 plotOutput("plt6", height = 500),
                 
                 p("See the list of the common genes",
                   style = "color:purple ; font-weight:600 ; font-size:120% ; margin-top:60px"),
                 
                 selectInput("sel.set", "Select you intersection", choices = c(), width = "650px"),
                 verbatimTextOutput("nb_common_genes"),
                 verbatimTextOutput("common_genes"),
                 
                 p("Genes specific to a condition", 
                   style = "color:purple ; font-weight:600 ; font-size:120% ; margin-top:20px"),
                 
                 selectInput("specific.genes", "Select the contrast", choices = c(), width = "650px"),
                 selectInput("specific.genes2", "Select the condition", choices = c(), width = "650px"),
                 verbatimTextOutput("nb_genes"),
                 verbatimTextOutput("specific_genes")
               )
             )),
    
    
    
    #### p3: Bar plot    ------------------------------------------------------
    tabPanel("Barplot",
             fluidRow(
               column(
                 width = 3,
                 p("Displaying parameters", style = "font-weight:600 ; color:darkred ; font-size:120%"),
                 numericInput("nUp", "Up genes to display", value = 10),
                 numericInput("nDown", "Down genes to display", value = 10)
               ),
               column(
                 width = 8,
                 offset = 1,
                 plotOutput("plt3", width = "100%", height = "700px")
               )
             ),
             p("Saving parameters", 
               style = "font-size:120% ; font-weight:600 ; background-color:gold ; color:darkred"),
             fluidRow(
               column(width = 4,
                      textInput("barplot_name", "Save as ...", value = "DEGs")),
               column(width = 4,
                      numericInput("width_bplot", "width", value = 8)),
               column(width = 4,
                      numericInput("height_bplot", "height", value = 10)),
               
               downloadButton("savebarplot","Save as pdf...", class = "btn-danger")
             )),
    
    #### p4: Dot plot   -------------------------------------------------------
    tabPanel("Dotplot",
             fluidRow(
               column(
                 width = 3,
                 p("Displaying parameters", style = "font-weight:600 ; color:darkred ; font-size:120%"),
                 selectInput("select.meta",
                             "Select your metadata",
                             choices = c()),
                 radioButtons("radio.select",
                              "Select your list of genes",
                              choices = c("Up Genes", "Down Genes"),
                              selected = "",
                              inline = T),
                 p("From the 'Genes' tabset resulting of the analysis you made",
                   style = "font-style:italic ; color:darkgray ; font-size:80%"),
                 sliderInput("slider.nbGenes",
                             "Range of genes to display",
                             min = 1,
                             max = 50,
                             value = c(1,10)),
                 p("The range here is set from the top (most up-regulated for Up genes, 
                   and most down-regulated for Down genes) to the bottom (least up-regulated 
                   for Up genes, and least down-regulated for Down genes)",
                   style = "font-style:italic ; color:darkgray ; font-size:80%"),
                 actionButton("act4", "Display", class = "btn-lg btn-success")
               ),
               column(
                 width = 9,
                 plotOutput("plt4", width = "100%", height = "700px")
               )
             ),
             p("Saving parameters", 
               style = "font-size:120% ; font-weight:600 ; background-color:gold ; color:darkred"),
             fluidRow(
               column(width = 4,
                      textInput("dotplot_name", "Save as ...", value = "DEGs")),
               column(width = 4,
                      numericInput("width_dplot", "width", value = 8)),
               column(width = 4,
                      numericInput("height_dplot", "height", value = 10)),
               
               downloadButton("savedotplot","Save as pdf...", class = "btn-danger")
             )),
    
    
    #### p5: Heatmap   --------------------------------------------------------
    tabPanel("Heatmap",
             p("Heatmap annotation", 
               style = "font-size:120% ; font-weight:600 ; background-color:gold ; color:darkred"),
             
             fluidRow(
               column(
                 width = 4,
                 p("Up genes selection", 
                   style = "font-weight:600 ; color:darkred ; font-size:130%"),
                 selectInput("hm.up.select", "Select up genes to display", choices = c(), 
                             multiple = T, selectize = F, width = "300px", size = 15),
                 p("", style = "margin-top:30px"),
                 actionButton("act6", "Select / Update", width = "200px"),
                 p("", style = "margin-top:40px"),
                 textOutput("txt3")
               ),
               column(
                 width = 4,
                 p("Down genes selection", 
                   style = "font-weight:600 ; color:darkred ; font-size:130%"),
                 selectInput("hm.down.select", "Select down genes to display", choices = c(), 
                             multiple = T, selectize = F, width = "300px", size = 15)
               ),
               column(
                 width = 4,
                 p("Set heatmap annotation", 
                   style = "font-weight:600 ; color:darkred ; font-size:130%"),
                 selectInput("hm.topannot.sel",
                             "Add annotation on the top",
                             choices = c(), multiple = T, selectize = F, width = "300px", size = 10),
                 p("", style = "margin-top:20px"),
                 actionButton("act5",
                              "metadata overview"),
                 p("", style = "margin-top:20px"),
                 tableOutput("tab3"),
                 p("", style = "margin-top:40px")
               )
             ),
             
             p("Heatmap plot", 
               style = "font-size:120% ; font-weight:600 ; background-color:gold ; color:darkred ; margin-top:30px"),
             
             fluidRow(
               column(
                 width = 3,
                 textInput("hm.title",
                           "Set a heatmap title",
                           value = "Heatmap representation"),
                 radioButtons("hm.rownames",
                              "Show gene names",
                              choices = c(TRUE, FALSE),
                              inline = T,
                              selected = T),
                 radioButtons("hm.rowdend",
                              "Show row dend",
                              choices = c(TRUE, FALSE),
                              inline = T,
                              selected = F),
                 p("", style = "margin-top:40px"),
                 actionButton("act7",
                              "Display the heatmap",
                              class = "btn-lg btn-success")
               ),
               column(
                 width = 9,
                 plotOutput("plt5", width = "100%", height = "800px")
               )
             ),
             p("Saving parameters", 
               style = "font-size:120% ; font-weight:600 ; background-color:gold ; color:darkred"),
             fluidRow(
               column(width = 4,
                      textInput("hm_name", "Save as ...", value = "heatmap")),
               column(width = 4,
                      numericInput("width_hm", "width", value = 8)),
               column(width = 4,
                      numericInput("height_hm", "height", value = 10)),
               
               downloadButton("savehm","Save as pdf...", class = "btn-danger")
             ))
  ),
  
  
  
  ##  IV- Discriminant markers   -------------------------------------------------------------------
  
  h2("4- Finding discriminant markers",
     style = "color:gold ; font-weight:700 ; background-color:black ; text-align:left ; margin-top:120px ; margin-bottom:20px"),
  
  p("Now, for each selected group (cluster), we'll apply the ROC (Receive Operating Characteristics) test to find 
    its discriminant genes according to their expression",
    style = "font-size:120% ; font-weight:600 ; color:darkgreen ; margin-bottom:50px"),
  
  fluidRow(
    column(
      width = 3,
      selectInput("ident.roc", 
                  "Select your cluster",
                  choices = c()),
      sliderInput(inputId = "log2fc.roc",
                  label = "logfc.threshold",
                  min = 0,
                  max = 3,
                  step = 0.1,
                  value = 1),
      sliderInput(inputId = "power.roc",
                  label = "Power cutoff",
                  min = 0,
                  max = 1,
                  value = 0.6,
                  step = 0.02),
      p("It is the statistical power of the ROC test, it represents how confident you must be about the results 
        (like the pvalue in the other test). 
        The closer it is to 1, the more confident you should be. A power under 0.5 is considered as not strong",
        style = "font-size:80% ; font-weight:500 ; color:gray ; font-style:italic ; margin-bottom:50px"),
      
      actionButton("act10", "Run the test" , class = "btn-lg btn-primary")
    ),
    
    column(
      width = 6,
      offset = 2,
      plotOutput("plt7")
    )
  ),
  
  fluidRow(
    column(
      width = 8,
      offset = 4,
      p("The results",
        style = "font-size:120% ; font-weight:600 ; background-color:gold ; 
          color:darkred ; margin-top:40px ; margin-bottom:20px"),
      dataTableOutput("table3"),
      p("'myAUC' is the column of interest here, it represent the area under the curve of the ROC test. 
          It mesures the ability of a gene to discriminate between the specified group (cluster) and all the 
          other groups. A value closer to 1 indicates perfect discimination, whereas a value close to 0.5 
          suggest no discrimination. Of course a value close to 0 means the gene may be discriminant elsewhere.",
        style = "font-size:90% ; font-weight:600 ; 
          color:black ; margin-top:10px ; margin-bottom:20px")
    )
  ),
  
  
  
  
  ##  V- PEA   -------------------------------------------------------------------
  
  ###  1- GO   -------------------------------------------------------------------
  
  h2("5- Pathway Enrichment Analysis",
     style = "color:gold ; font-weight:700 ; background-color:black ; text-align:left ; margin-top:120px ; margin-bottom:20px"),
  
  p("Now with the list of markers we've filled, let's proceed to the PE Analysis",
    style = "font-size:120% ; font-weight:600 ; color:navy"),
  
  navlistPanel(
    id = "PEA",
    widths = c(2,10),
    tabPanel("GO analysis",
             p("Let's proceed to GO analysis",
               style = "font-size:130% ; font-weight:600 ; color:darkred ; margin-bottom:20px"),
             
             fluidRow(
               column(
                 width = 3,
                 selectInput("go.list.genes", "Select gene list", choices = c()),
                 selectInput("go.list.genes.reg", "Select the regulation", choices = c("Up","Down","All")),
                 textInput("go.keytype", "Key Type", value = "SYMBOL"),
                 selectInput("go.ont", "Ontology", choices = c("BP","MF", "CC"), selected = "BP"),
                 radioButtons("go.maxgeneset", "Limit the geneset size", 
                              choices = c("<500","<1000","<5000","Inf"), inline = T),
                 actionButton("act11", "RUN GO", class = "btn-lg btn-primary"),
                 p("", style = "margin-top:20px"),
                 actionButton("addGO", "Add to list", class = "btn-lg btn-success"),
                 p("Saved objects :", style = "margin-top:20px ; margin-bottom:10px ; color:darkred"),
                 verbatimTextOutput("showGOobj")
               ),
               column(
                 width = 9,
                 dataTableOutput("tab4"),
                 p("You can access to the gene list for each term by entering the GO ID as a key down here !",
                   style = "font-size:100% ; font-weight:600 ; color:darkgreen ; margin-bottom:20px ; margin-top:30px"),
                 
                 textInput("goID", "Enter GO ID", value = ""),
                 actionButton("act12", "Find term genes", class = "btn-sm btn-success"),
                 verbatimTextOutput("term_genes")
               )
             ),
             
             fluidRow(
               p("And now, let's visualize the results",
                 style = "font-size:130% ; font-weight:600 ; color:darkred ; margin-bottom:20px ; margin-top:40px"),
               
               ####  p1: Dot plot   --------------------------------------------
               tabsetPanel(
                 tabPanel("Dotplot",
                          fluidRow(
                            column(
                              width = 6,
                              selectInput("dotplot.show.cat", "Select terms", choices = c(), 
                                          size = 15, multiple = T, selectize = F, width = "550px")
                            ),
                            
                            column(
                              width = 5,
                              p("Selected Pathways are:",
                                style = "font-size:120% ; font-weight:600 ; margin-bottom:20px ; color:darkgreen"),
                              verbatimTextOutput("dotp_selected_terms"),
                              
                              p("", style = "margin-top:40px"),
                              
                              actionButton("act14", "Plot Dotplot", class = "btn-lg btn-success")
                            )
                          ),
                          
                          p("Dot plot",
                            style = "font-weight:600 ; color:darkred ; background-color:gold ; 
                                margin-top:50px ; margin-bottom:30px ; font-size:130% ; text-align:center"),
                          plotOutput("plt8")),
                 
                 ####  p2: CNetwork plot   -------------------------------------
                 
                 tabPanel("Network plot",
                          p("Now, to concider the potentially biological complexities in which a gene may belong 
                              to multiple annotation categories and provide information of numeric changes if 
                              available, the cnetplot() function will be used to extract complex associations. 
                              It depicts the linkages of genes and biological concepts as a network.",
                            style = "font-weight:600 ; color:black ; margin-top:20px ; margin-bottom:30px"),
                          
                          fluidRow(
                            column(
                              width = 6,
                              selectInput("go.cnet.sel", "Select terms", choices = c(), 
                                          size = 15, multiple = T, selectize = F, width = "550px"),
                              
                              radioButtons("go.cnet.layout", "Select the layout", selected = "kk",
                                           choices = c("kk","circle","grid","fr"), inline = T),
                              
                              numericInput("gene_lbl_size", label = "Genes label size", value = 1, step = 0.1),
                              numericInput("categ_node_size", label = "Terms label size", value = 1.4, step = 0.1)
                            ),
                            column(
                              width = 5,
                              p("Selected Pathways are:",
                                style = "font-size:120% ; font-weight:600 ; margin-bottom:20px ; color:darkgreen"),
                              verbatimTextOutput("cnet_selected_terms"),
                              
                              p("", style = "margin-top:40px"),
                              
                              actionButton("act13", "Plot Cnetplot", class = "btn-lg btn-success")
                              
                              
                            )
                          ),
                          
                          p("Network plot",
                            style = "font-weight:600 ; color:darkred ; background-color:gold ;
                                margin-top:50px ; margin-bottom:30px ; font-size:130% ; text-align:center"),
                          
                          plotOutput("plt9")),
                 
                 
                 ####  p3: Corr plot   ---------------------------------------------
                 
                 tabPanel("Correlation plot",
                          p("A correlation score between the terms is calculated using Jaccard's similarity 
                            index (JC). A hierarchical clustering will be then made using that JC scores, and 
                            also for the clustering, several method can be applied !",
                            style = "font-weight:600 ; color:black ; margin-top:20px ; margin-bottom:30px"),
                          
                          fluidRow(
                            column(
                              width = 6,
                              selectInput("terms_corrplot", "Select terms", choices = c(), 
                                          size = 20, multiple = T, selectize = F, width = "550px")
                            ),
                            column(
                              width = 4,
                              offset = 1,
                              radioButtons("corr.methode", "Select correlation method", 
                                           choices = c("circle","square","ellipse","number","shade","color","pie")),
                              selectInput("corr.hclust","Select clustering method",
                                          choices = c("complete", "ward", "ward.D", "ward.D2", "single", "average",
                                                      "mcquitty", "median", "centroid")),
                              p("", style = "margin-top:30px"),
                              actionButton("act14", "Plot the correlation", class = "btn-lg btn-success")
                            )
                          ),
                          
                          fluidRow(
                            
                            p("Correlation plot",
                              style = "font-weight:600 ; color:darkred ; margin-top:20px ; 
                            margin-bottom:30px ; font-size:120% ; background-color:gold"),
                            
                            plotOutput("plt10")
                          )),
                 
                 
                 
                 ####  p4: Treeplot   ------------------------------------------
                 
                 tabPanel("Tree plot",
                          p("A hierarchical clustering between the selected terms below to see the connections !",
                            style = "font-weight:600 ; color:black ; margin-top:20px ; margin-bottom:30px"),
                          
                          fluidRow(
                            column(
                              width = 6,
                              selectInput("terms_treeplot", "Select terms", choices = c(), 
                                          size = 20, multiple = T, selectize = F, width = "550px")
                            ),
                            column(
                              width = 5,
                              p("Selected Pathways are:",
                                style = "font-size:120% ; font-weight:600 ; margin-bottom:20px ; color:darkgreen"),
                              verbatimTextOutput("cnet_selected_terms")
                            )
                          ),
                          
                          fluidRow(
                            column(
                              width = 3,
                              numericInput("treeplot_cex_category", "Hexpand params", value = 0.1, step = 0.01),
                              p("", style = "margin-bottom:30px"),
                              actionButton("act15", "Plot the correlation", class = "btn-lg btn-success")
                            ),
                            column(
                              width = 3,
                              numericInput("treeplot_nCluster", "Split into n clusters", value = 5, step = 1)
                            ),
                            column(
                              width = 3,
                              selectInput("treeplot_Clust_method", "Clustering method", 
                                          choices = c("ward.D","complete","ward.D2","median","single"),
                                          selected = "ward.D")
                            ),
                            column(
                              width = 3,
                              checkboxInput("treeplot_Clust_highlight", "Highlight clusters", value = T)
                            )
                          ),
                          
                          fluidRow(
                            p("Tree plot",
                              style = "font-weight:600 ; color:darkred ; margin-top:20px ; 
                            margin-bottom:30px ; font-size:120% ; background-color:gold"),
                            
                            plotOutput("plt11")
                          )),
                 
                 
                 
                 
                 
                 ####  p5: Enrichment map   ------------------------------------------
                 
                 tabPanel("Enrich map",
                          p("Enrichment map organizes enriched terms into a network with edges connecting 
                            overlapping gene sets. In this way, mutually overlapping gene sets are tend to 
                            cluster together, making it easy to identify functional module.",
                            style = "font-weight:600 ; color:black ; margin-top:20px ; margin-bottom:30px"))
               )
             ))
  )
  
  
  
)














#  Function    ##########################

server <- function(input, output, session){
  
  ##  I-  Import data   ------------------------------------------------------
  
  seurat <- reactive({
    req(input$seurat)
    seurat <- readRDS(input$seurat$datapath)
  })
  
  output$text <- renderPrint({
    req(seurat())
    print(seurat())
  })
  
  ##  II- Data manipulation   ------------------------------------------------
  
  observeEvent(seurat(),{
    seurat <- seurat()
    # clustering_columns <- c(grep("res.0",names(seurat@meta.data), value = T),
    #                         "hash.ID", "Conditions", "Clusters")
    
    clustering_columns <- names(seurat@meta.data)[sapply(seurat@meta.data, function(x) !is.numeric(x))]
    
    updateSelectInput(session, "clusters", choices = clustering_columns, selected = "Conditions")
    updateSelectInput(session, "select.meta", choices = clustering_columns, selected = "Conditions")
  })
  
  summary_table <- eventReactive(input$act, {
    req(seurat(), input$clusters)
    seurat <- seurat()
    seurat@meta.data %>% 
      dplyr::group_by(!!sym(input$clusters)) %>% 
      dplyr::summarise(Nb_of_cells = n())
  })
  
  
  output$tab1 <- renderTable({
    summary_table()
  })
  
  
  output$plt1 <- renderPlot({
    req(seurat(), input$clusters)
    seurat <- seurat()
    seurat@meta.data %>%
      ggplot()+
      geom_point(aes(x= umap_1,
                     y= umap_2,
                     colour= .data[[input$clusters]]),
                 size = 1)+
      labs(colour="Clusters")+
      scale_colour_manual(values = MyPalette)+
      theme_bw()+
      theme(
        axis.title = element_text(face= "bold", size= 15),
        axis.text = element_text(size= 12),
        legend.title = element_text(face= "bold", size= 15, colour= "navy"),
        legend.text = element_text(size= 13, colour= "black")
      )+
      guides(colour = guide_legend(override.aes = list(size=3)))
  })
  
  ##  III-  DEA   ------------------------------------------------------------
  
  ### 1-  Setting the contrasts   --------------------------------------------
  
  observeEvent(input$clusters, {
    req(seurat(), input$clusters)
    seurat <- seurat()
    
    cluster_levels <- unique(seurat@meta.data[[input$clusters]])
    updateSelectInput(session, "ident.1", choices = cluster_levels, selected = cluster_levels[1])
    updateSelectInput(session, "ident.2", choices = cluster_levels, selected = cluster_levels[2])
    updateSelectInput(session, "ident.roc", choices = cluster_levels, selected = cluster_levels[1])
  })
  
  
  output$txt2 <- renderText({
    req(seurat(), input$ident.1, input$ident.2)
    paste0("You are comparing [", input$ident.1, "] vs [", 
           input$ident.2, "]")
  })
  
  
  output$plt2 <- renderPlot({
    req(seurat(), input$clusters, input$ident.1, input$ident.2)
    seurat <- seurat()
    
    seurat@meta.data %>% 
      mutate(cond = ifelse(.data[[input$clusters]] == input$ident.1 | 
                             .data[[input$clusters]] == input$ident.2, T, F)) %>% 
      ggplot()+
      geom_point(aes(x= umap_1,
                     y= umap_2,
                     colour= .data[[input$clusters]],
                     alpha= cond),
                 size= 0.9)+
      labs(colour="Clusters",
           alpha="")+
      scale_colour_manual(values = MyPalette)+
      scale_alpha_manual(values = c("TRUE"=1, "FALSE"=0.02), guide = "none")+
      theme_bw()+
      theme(
        axis.title = element_text(face= "bold", size= 15),
        axis.text = element_text(size= 12),
        legend.title = element_text(face= "bold", size= 15, colour= "navy"),
        legend.text = element_text(size= 13, colour= "black")
      )+
      guides(colour = guide_legend(override.aes = list(size=3)))
    
  })
  
  DEA_message <- reactiveVal("")
  
  observeEvent(input$act2, {
    
    showNotification("DEA is running ...", type = "message")
    
    message("The following parameters are used: ")
    
    message(paste("clustering = ", input$clusters))
    message(paste("ident.1 = ", input$ident.1))
    message(paste("ident.2 = ", input$ident.2))
    message(paste("log2fc = ", input$log2fc))
    message(paste("min.pct = ", input$min.pct))
    message(paste("only.pos = ", input$only.pos))
    
  })
  
  
  ### 2-  Running DEA   ------------------------------------------------------
  
  DEA_results <- eventReactive(input$act2, {
    
    seurat <- seurat()
    Idents(seurat) <- as.factor(seurat@meta.data[[input$clusters]])
    
    seurat %>% 
      FindMarkers(
        ident.1 = input$ident.1,
        ident.2 = input$ident.2,
        only.pos = input$only.pos,
        min.pct = input$min.pct,
        logfc.threshold = input$log2fc,
        min.diff.pct = input$diffpct_filt,
        verbose = TRUE
      ) %>%
      dplyr::mutate(Diff_pct = pct.1 - pct.2,
                    avg_log2FC = round(avg_log2FC,3)) %>%
      dplyr::filter(p_val_adj < 0.05) %>% 
      dplyr::filter(pct.1 >= input$pct.1_filt[1] & pct.1 <= input$pct.1_filt[2],
                    pct.2 >= input$pct.2_filt[1] & pct.2 <= input$pct.2_filt[2],
                    avg_log2FC >= input$logfc_filt[1] & avg_log2FC <= input$logfc_filt[2] |
                      avg_log2FC >= input$logfc_filt2[1] & avg_log2FC <= input$logfc_filt2[2]) %>%
      dplyr::arrange(desc(avg_log2FC)) %>%
      dplyr::select(avg_log2FC, pct.1, pct.2, Diff_pct, p_val_adj) %>%
      rownames_to_column(var = "Genes") 
    
  })
  
  
  observeEvent(DEA_results(),{
    
    res <- DEA_results()
    max_fc <- max(res[["avg_log2FC"]])
    min_fc <- min(res[["avg_log2FC"]])
    
    DEA_message("Analysis completed succesfully !")
    
    message("Analysis completed succesfully !")
    
    updateSliderInput(session, "logfc_filt", min = 0, max = 10, value = c(0, max_fc), step = 0.1)
    updateSliderInput(session, "logfc_filt2", min = -10, max = 0, value = c(min_fc, 0), step = 0.1)
  })
  
  
  output$dea_log <- renderText({
    req(DEA_message())
    DEA_message <- DEA_message()
    print(DEA_message)
  })
  
  
  
  ### 3-  Displaying results    ----------------------------------------------
  
  ####  p1: Datatable ------------------------------------------------------------
  
  output$tab2 <- renderDataTable({
    DEA_results() %>%
      datatable()
  })
  
  
  
  
  ####  p2: Genes list  ----------------------------------------------------------
  List_markers <- reactiveVal(list())
  Up_Genes <- reactiveVal(character(0))
  Down_Genes <- reactiveVal(character(0))
  
  ##### Displaying the lists    ------------------------------------------------
  observeEvent(DEA_results(), {
    req(DEA_results())
    results <- DEA_results() 
    
    Up_Genes(results %>% 
               dplyr::filter(avg_log2FC > 0) %>% 
               dplyr::pull(Genes))
    Down_Genes(results %>% 
                 dplyr::filter(avg_log2FC < 0) %>% 
                 dplyr::arrange(avg_log2FC) %>% 
                 dplyr::pull(Genes))
  })
  
  
  output$genes1 <- renderPrint({
    Up_Genes()
  })
  
  output$genes2 <- renderPrint({
    Down_Genes()
  })
  
  
  ##### Filling the cumulative list   --------------------------------------------
  observeEvent(input$act8, {
    req(Up_Genes(), Down_Genes(), DEA_results())
    
    res <- DEA_results()
    fc <- res[["avg_log2FC"]] %>% setNames(res[["Genes"]])
    
    tmp_list <- list(
      Up = Up_Genes(),
      Down = Down_Genes(),
      All = c(Up_Genes(), Down_Genes()),
      All_FC = fc
    )
    
    current_list <- List_markers()
    contrast_name <- input$txt4
    current_list[[contrast_name]] <- tmp_list
    List_markers(current_list)
    
    showNotification(paste0("Contrast ",contrast_name, " added to the list"), type = "message")
  })
  
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
  
  output$showMarkers <- renderPrint({
    req(List_markers())
    List_markers <- List_markers()
    str(List_markers)
  })
  
  
  ##### Comparison part     ----------------------------------------------------
  observeEvent(List_markers(), {
    req(List_markers())
    List_markers <- List_markers()
    
    updateSelectInput(session, "contrasts", choices = names(List_markers))
    
  })
  
  ######  Venn diagram    ------------------------------------------------------
  venn_obj <- eventReactive(input$act9, {
    req(List_markers())
    List_markers <- List_markers()
    genes <- switch(input$up_down,
                    "Up genes" = "Up",
                    "Down genes" = "Down",
                    "All genes" = "All")
    
    sets <- list()
    for(set in input$contrasts){
      sets[[set]] <- List_markers[[set]][[genes]]
    }
    
    Display_Venn(sets, set.names = input$contrasts)
    
  })
  
  venn_plot <- eventReactive(venn_obj(), {
    req(venn_obj())
    venn_obj <- venn_obj()
    
    venn_obj[["plot"]]
  })
  
  output$plt6 <- renderPlot({
    req(venn_plot())
    venn_plot()
  }, height = 500)
  
  
  ######  Displaying genes    --------------------------------------------------
  observeEvent(venn_obj(), {
    req(venn_obj())
    venn_obj <- venn_obj()
    
    updateSelectInput(session, "sel.set", 
                      choices = names(venn_obj[["intersections"]]), selected = "")
    updateSelectInput(session, "specific.genes", 
                      choices = names(venn_obj[["group_specific"]]), selected = "")
  })
  
  observeEvent(input$specific.genes, {
    req(venn_obj())
    venn_obj <- venn_obj()
    
    updateSelectInput(session, "specific.genes2", 
                      choices = names(venn_obj[["group_specific"]][[input$specific.genes]]), 
                      selected = "")
  })
  
  output$nb_common_genes <- renderPrint({
    req(venn_obj())
    venn_obj <- venn_obj()
    print(length(venn_obj[["intersections"]][[input$sel.set]]))
  })
  
  output$common_genes <- renderPrint({
    req(venn_obj())
    venn_obj <- venn_obj()
    print(venn_obj[["intersections"]][[input$sel.set]])
  })
  
  output$nb_genes <- renderPrint({
    req(venn_obj())
    venn_obj <- venn_obj()
    print(length(venn_obj[["group_specific"]][[input$specific.genes]][[input$specific.genes2]]))
  })
  
  output$specific_genes <- renderPrint({
    req(venn_obj())
    venn_obj <- venn_obj()
    print(venn_obj[["group_specific"]][[input$specific.genes]][[input$specific.genes2]])
  })
  
  
  
  #### p3: Bar plot  -------------------------------------------------------------
  
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
        legend.title = element_text(face= "bold", size= 15, colour= "navy"),
        legend.text = element_text(size= 13, colour= "black"),
        plot.title = element_text(size = 17, colour = "darkred", face = "bold")
      )+
      guides(colour = guide_legend(override.aes = list(size=3)))
  }, width = 700, height = 700)
  
  
  #### p4: Dot plot --------------------------------------------------------------
  
  observeEvent(Up_Genes(), {
    req(Up_Genes(), Down_Genes())
    len_UpGenes <- length(Up_Genes())
    len_DownGenes <- length(Down_Genes())
    
    message("Nb of Up Genes = ", 
            len_UpGenes, " | ",
            Up_Genes()[1], ", ", Up_Genes()[2], ", ", Up_Genes()[3], " ...")
    message("Nb of Down Genes = ", 
            len_DownGenes, " | ",
            Down_Genes()[1], ", ", Down_Genes()[2], ", ", Down_Genes()[3], " ...")
  })
  
  
  observeEvent(input$radio.select, {
    req(Up_Genes(), Down_Genes())
    len_UpGenes <- length(Up_Genes())
    len_DownGenes <- length(Down_Genes())
    
    if(input$radio.select == "Up Genes"){
      updateSliderInput(session, "slider.nbGenes", max = len_UpGenes)
    } else {
      updateSliderInput(session, "slider.nbGenes", max = len_DownGenes)
    }
  })
  
  
  dot.plot <- eventReactive(input$act4, {
    req(seurat(), Up_Genes(), Down_Genes(), DEA_results())
    seurat <- seurat()
    Up_Genes <- Up_Genes()
    Down_Genes <- Down_Genes()
    
    genes <- if(input$radio.select == "Up Genes"){
      Up_Genes[input$slider.nbGenes[1]:input$slider.nbGenes[2]]
    } else {
      Down_Genes[input$slider.nbGenes[1]:input$slider.nbGenes[2]]
    }
    
    n <- length(unique(seurat[[input$select.meta]]))
    m <- sapply(unique(seurat[[input$select.meta]]),nchar)
    angle = 55
    
    if(n <= 7 && median(m) <= 8){
      angle = 0
    }
    
    if(length(genes) > 0){
      DotPlot(seurat,
              features = genes,
              group.by = input$select.meta)+
        RotatedAxis()+
        coord_flip()+
        theme_bw()+
        labs(x="",
             y="")+
        scale_color_gradient(low = "orange", high = "navy")+
        theme(
          axis.title = element_text(face= "bold", size= 15),
          axis.text.y = element_text(size= 16, face = "bold"),
          axis.text.x = element_text(size= 16, face = "bold", angle = angle),
          legend.title = element_text(face= "bold", size= 15, colour= "navy"),
          legend.text = element_text(size= 13, colour= "black"),
          plot.title = element_text(size = 17, colour = "darkred", face = "bold")
        )
    } else {
      showNotification("No genes to display in this selected range",
                       type = "error")
    }
  })
  
  
  output$plt4 <- renderPlot({
    dot.plot()
  })
  
  
  #### p5: Heatmap  --------------------------------------------------------------
  
  ##### - Expression Assay preparation :    --------------------------------------
  
  observeEvent(DEA_results(), {
    req(Up_Genes(), Down_Genes())
    up <- Up_Genes()
    down <- Down_Genes()
    
    updateSelectInput(session, "hm.up.select", choices = up)
    updateSelectInput(session, "hm.down.select", choices = down)
  })
  
  
  ExpressionAssay <- eventReactive(input$act6, {
    req(seurat(), Up_Genes(), Down_Genes())
    seurat <- seurat()
    up <- Up_Genes()
    down <- Down_Genes()
    
    FetchData(seurat, vars = c(input$hm.up.select, input$hm.down.select)) %>% 
      t() %>% 
      as.matrix() %>% 
      log1p()
    
  })
  
  output$txt3 <- renderText({
    req(ExpressionAssay())
    
    exp <- ExpressionAssay()
    if(nrow(exp) > 0){
      print(paste0("The expression Assay has been successfully created 
                   with the selected amount of genes : ", nrow(exp), " genes"))
    } else {
      print("None of the genes have been found !!")
    }
  })
  
  
  ##### - Heatmap annotations :   ----------------------------------------------
  
  observeEvent(DEA_results(), {
    req(seurat())
    seurat <- seurat()
    topannot <- names(seurat@meta.data)[sapply(seurat@meta.data, function(x) !is.numeric(x))]
    
    updateSelectInput(session, "hm.topannot.sel", choices = topannot)
  })
  
  ## 2-1  metadata :
  
  hm.metadata <- eventReactive(input$act5, {
    req(seurat())
    seurat <- seurat()
    
    MyVars <- input$hm.topannot.sel
    
    
    if(length(MyVars) == 1){
      meta <- as.data.frame(setNames(seurat@meta.data[[input$hm.topannot.sel]], 
                                      rownames(seurat@meta.data)))
      colnames(meta) <- input$hm.topannot.sel
      
    } else if(length(MyVars) > 1){
      meta <- seurat@meta.data[,MyVars]
    }
    
    showNotification("Metadata created !", type = "message")
    meta
  })
  
  output$tab3 <- renderTable({
    req(hm.metadata())
    head(hm.metadata(),5)
  })
  
  ## 2-2  Top annotations :
  
  TopAnnotation <- eventReactive(hm.metadata(), {
    req(hm.metadata())
    hm.metadata <- hm.metadata()
    
    ### Setting colors & top annotations
    if(ncol(hm.metadata) == 1){
      
      col1 <- Assign_colors(object = hm.metadata[[input$hm.topannot.sel]],
                            pal = MyPalette)
      
      HeatmapAnnotation(df = hm.metadata,
                        col = setNames(list(col1),
                                       c(input$hm.topannot.sel)),
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
      
      col1 <- Assign_colors(object = hm.metadata[[input$hm.topannot.sel[1]]],
                            pal = MyPalette)
      col2 <- Assign_colors(object = hm.metadata[[input$hm.topannot.sel[2]]],
                            pal = MyPalette[which(!(MyPalette %in% col1))])
      
      HeatmapAnnotation(df = hm.metadata,
                        col = setNames(list(col1, col2),
                                       c(input$hm.topannot.sel[1],
                                         input$hm.topannot.sel[2])),
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
      
      
    } else if(ncol(hm.metadata) == 3){
      
      col1 <- Assign_colors(object = hm.metadata[[input$hm.topannot.sel[1]]],
                            pal = MyPalette)
      col2 <- Assign_colors(object = hm.metadata[[input$hm.topannot.sel[2]]],
                            pal = MyPalette[which(!(MyPalette %in% col1))])
      col3 <- Assign_colors(object = hm.metadata[[input$hm.topannot.sel[3]]],
                            pal = MyPalette[which(!(MyPalette %in% c(col1,col2)))])
      
      HeatmapAnnotation(df = hm.metadata,
                        col = setNames(list(col1, col2, col3),
                                       c(input$hm.topannot.sel[1],
                                         input$hm.topannot.sel[2],
                                         input$hm.topannot.sel[3])),
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
      
    } else if(ncol(hm.metadata) > 3){
      showNotification("The heatmap can't take more than 3 top annotations ! selecting the 3 first !",
                       type = "warning", duration = 10)
      
      hm.metadata <- hm.metadata[,1:3]
      col1 <- Assign_colors(object = hm.metadata[[input$hm.topannot.sel[1]]],
                            pal = MyPalette)
      col2 <- Assign_colors(object = hm.metadata[[input$hm.topannot.sel[2]]],
                            pal = MyPalette[which(!(MyPalette %in% col1))])
      col3 <- Assign_colors(object = hm.metadata[[input$hm.topannot.sel[3]]],
                            pal = MyPalette[which(!(MyPalette %in% c(col1,col2)))])
      
      HeatmapAnnotation(df = hm.metadata,
                        col = setNames(list(col1, col2, col3),
                                       c(input$hm.topannot.sel[1],
                                         input$hm.topannot.sel[2],
                                         input$hm.topannot.sel[3])),
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
    } 
    
    
    
  })
  
  ## 2-3  Right annotations :
  
  RightAnnotation <- eventReactive(ExpressionAssay(), {
    req(ExpressionAssay(), DEA_results())
    exp <- ExpressionAssay()
    dea_res <- DEA_results()
    
    rowAnnotation(
      "log2FC" = anno_barplot(dea_res$avg_log2FC[which(dea_res$Genes %in% rownames(exp))], 
                              axis = TRUE, 
                              bar_width = 1, 
                              border = TRUE,
                              cex = 1))
    
  })
  
  
  ##### - Heatmap plot :    -----------------------------------------------------
  hm <- eventReactive(input$act7, {
    req(ExpressionAssay(), RightAnnotation(), TopAnnotation())
    exp <- ExpressionAssay()
    RightAnnotation <- RightAnnotation()
    TopAnnotation <- TopAnnotation()
    
    showNotification("Heatmap on its way !!!", duration = 10, type = "message")
    message("Heatmap : ", dim(exp)[1]," genes | ", dim(exp)[2], " cells")
    
    ColExpr <- colorRamp2(breaks = c(0, max(exp)/3, max(exp)*2/3, max(exp)), 
                          colors = c("white", "gold", "cadetblue","darkred"))
    
    ComplexHeatmap::Heatmap(exp,
                            top_annotation = TopAnnotation,
                            right_annotation = RightAnnotation,
                            cluster_rows = TRUE,
                            cluster_columns = TRUE,
                            show_row_names = ifelse(input$hm.rownames == T, T, F),
                            show_column_names = F,
                            col = ColExpr,
                            show_row_dend = ifelse(input$hm.rowdend == T, T, F),
                            row_dend_side = "left",
                            row_dend_width = unit(2,"cm"),
                            show_column_dend = F,
                            row_names_side = "left",
                            row_names_gp = gpar(fontface = "bold",
                                                fontsize = 12,
                                                col = "darkred"),
                            column_title = input$hm.title,
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
  
  
  output$plt5 <- renderPlot({
    req(hm())
    hm <- hm()
    ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom")
  })
  
  
  
  
  ##  IV- Discriminant markers   ---------------------------------------------------------------
  
  
  output$plt7 <- renderPlot({
    req(seurat(), input$clusters, input$ident.roc)
    seurat <- seurat()
    
    seurat@meta.data %>% 
      mutate(cond = ifelse(.data[[input$clusters]] == input$ident.roc, T, F)) %>% 
      ggplot()+
      geom_point(aes(x= umap_1,
                     y= umap_2,
                     colour= .data[[input$clusters]],
                     alpha= cond),
                 size= 0.9)+
      labs(colour="Clusters",
           alpha="")+
      scale_colour_manual(values = MyPalette)+
      scale_alpha_manual(values = c("TRUE"=1, "FALSE"=0.02), guide = "none")+
      theme_bw()+
      theme(
        axis.title = element_text(face= "bold", size= 15),
        axis.text = element_text(size= 12),
        legend.title = element_text(face= "bold", size= 15, colour= "navy"),
        legend.text = element_text(size= 13, colour= "black")
      )+
      guides(colour = guide_legend(override.aes = list(size=3)))
    
  })
  
  
  
  Results.disc.genes <- eventReactive(input$act10, {
    req(seurat())
    seurat <- seurat()
    Idents(seurat) <- as.factor(seurat@meta.data[[input$clusters]])
    
    FindMarkers(seurat,
                ident.1 = input$ident.roc, 
                logfc.threshold = input$log2fc.roc, 
                test.use = "roc") %>% 
      dplyr::filter(power >= input$power.roc) %>% 
      dplyr::arrange(desc(myAUC)) %>% 
      head(20)
    
  })
  
  observeEvent(Results.disc.genes(), {
    req(Results.disc.genes())
    Results.disc.genes <- Results.disc.genes()
    
    if(nrow(Results.disc.genes) == 0){
      showNotification("No discriminant genes found !!", type = "warning")
    } else {
      showNotification("Analysis completed successfully !!", type = "message")
    }
  })
  
  output$table3 <- renderDataTable({
    req(Results.disc.genes())
    Results.disc.genes() 
  }, options = list(pageLength = 5, scrollX = T))
  
  
  
  
  ##  V- PEA   -------------------------------------------------------------------
  
  ###  1- GO   -------------------------------------------------------------------
  
  List_GO <- reactiveVal(list())
  
  observeEvent(List_markers(), {
    
    req(List_markers())
    List_markers <- List_markers()
    
    updateSelectInput(session, "go.list.genes", choices = names(List_markers))
  })
  
  
  en.go <- eventReactive(input$act11, {
    
    req(List_markers())
    List_markers <- List_markers()
    MyGenes <- List_markers[[input$go.list.genes]][[input$go.list.genes.reg]]
    
    maxGSSize <- switch(input$go.maxgeneset,
                        "<500" = 500,
                        "<1000" = 1000,
                        "<5000" = 5000,
                        "Inf" = Inf)
    
    showNotification("GO is running .....", type = "message")
    
    enrichGO(gene = MyGenes,
             OrgDb = org.Hs.eg.db,
             keyType = input$go.keytype,
             ont = input$go.ont,
             qvalueCutoff = 0.1,
             pvalueCutoff = 0.05,
             readable = T,
             maxGSSize = maxGSSize,
             minGSSize = 10)
    
    
  })
  
  observeEvent(input$addGO, {
    req(en.go(), input$go.list.genes, input$go.list.genes.reg)
    
    current_list <- List_GO()
    contrast_name <- paste0(input$go.list.genes,"_",input$go.list.genes.reg)
    
    if(contrast_name %in% names(current_list)){
      showNotification(paste0("The contrast ", contrast_name, " already exists in the list!"), 
                       type = "warning")
    } else {
      current_list[[contrast_name]] <- en.go()
      List_GO(current_list)
      showNotification(paste0("GO analysis for ", contrast_name, " added to the list!"), type = "message")
    }
  })
  
  
  output$showGOobj <- renderPrint({
    req(List_GO())
    List_GO <- List_GO()
    
    if(length(List_GO) > 0){
      for(i in names(List_GO)){
        print(i)
      }
    }
    
  })
  
  go.results <- eventReactive(en.go(), {
    
    req(en.go())
    
    en.go <- en.go()
    
    en.go@result %>% 
      dplyr::filter(p.adjust < 0.05) %>% 
      dplyr::mutate(RichFactor = round(Count / as.numeric(sub("/\\d+","",BgRatio)),6)) %>% 
      dplyr::select(ID, Description, GeneRatio, BgRatio, Count, RichFactor, p.adjust, geneID)
    
  })
  
  
  
  
  output$tab4 <- renderDataTable({
    req(go.results())
    go.results() %>% 
      dplyr::select(ID, Description, GeneRatio, BgRatio, RichFactor, p.adjust) %>% 
      datatable(options = list(pageLength = 5, scrollX = TRUE))
  })
  
  
  term_genes <- reactiveVal(character(0))
  
  observeEvent(input$act12, {
    req(go.results())
    go.results <- go.results()
    
    term_genes(extract_term_genes(go.results, input$goID))
    
    message("GO ID entered : ", input$goID)
    message("Description of the term : ", go.results[["Description"]][which(go.results[["ID"]] == input$goID)])
    message("Number of genes : ", length(term_genes()))
    
  })
  
  output$term_genes <- renderPrint({
    req(term_genes())
    print(term_genes())
  })
  
  
  
  
  
  
  #### update selection      ----------------------------------------------------
  
  observeEvent(go.results(), {
    req(en.go(), go.results())
    
    go.results <- go.results()
    
    updateSelectInput(session, "go.cnet.sel", choices = go.results[["Description"]])
    updateSelectInput(session, "dotplot.show.cat", choices = go.results[["Description"]])
    updateSelectInput(session, "terms_corrplot", choices = go.results[["Description"]], selected = "ward.D")
    updateSelectInput(session, "terms_treeplot", choices = go.results[["Description"]], selected = "ward.D")
    
  })
  
  
  
  #### Dot plot    --------------------------------------------------------------
  
  output$dotp_selected_terms <- renderPrint({
    req(input$dotplot.show.cat)
    
    selected_terms <- input$dotplot.show.cat
    
    cat(selected_terms, sep = "\n")
  })
  
  
  Dplot <- eventReactive(input$act14, {
    req(go.results())
    
    go.results() %>% 
      dplyr::filter(Description %in% input$dotplot.show.cat) %>% 
      ggplot(aes(x= RichFactor, y= fct_reorder(Description, RichFactor)))+
      geom_segment(aes(xend= 0, yend= Description))+
      geom_point(aes(color= p.adjust, size= Count))+
      scale_color_viridis_c(guide = guide_colorbar(reverse = T))+
      scale_size_continuous(range = c(3,12))+
      theme_linedraw()+
      theme(panel.grid = element_blank(),
            panel.border = element_blank(),
            plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
            plot.margin = unit(c(0.2,0.1,0.1,0), "in"),
            axis.title.x = element_text(size = 15, face = "bold", margin = margin(t=0.2, unit = "in")),
            axis.title.y = element_blank(),
            axis.text = element_text(size = 14, face = "bold"),
            legend.title = element_text(size = 15, face = "bold", margin = margin(b=0.2, unit = "in")),
            legend.text = element_text(size = 12),
            legend.box.margin = margin(l=0.2, unit = "in"))+
      labs(title = "GO Enriched Terms")
    
  })
  
  output$plt8 <- renderPlot({
    req(Dplot())
    Dplot()
    
    
  },height = 700, width = 900)
  
  
  
  #### Cnet plot    --------------------------------------------------------------
  
  output$cnet_selected_terms <- renderPrint({
    req(input$go.cnet.sel)
    
    selected_terms <- input$go.cnet.sel
    
    cat(selected_terms, sep = "\n")
  })
  
  CNetPlot <- eventReactive(input$act13, {
    req(en.go(), List_markers())
    en.go <- en.go()
    List_markers <- List_markers()
    
    fc <- List_markers[[input$go.list.genes]][["All_FC"]]
    selected_terms <- input$go.cnet.sel
    selected_layout <- input$go.cnet.layout
    gene_label_size <- input$gene_lbl_size
    term_label_size <- input$categ_node_size
    
    if(is.null(fc) | is.null(selected_terms) | is.null(selected_layout)){
      warning("Either the fold change values, the selected terms or the selected layout 
              is NULL, check out your code !")
    }
    
    if(class(en.go) == "enrichResult"){
      message("checking your object !")
      message("enrichResult object class confirmed !")
    }
    
    if(!is.null(List_markers)){
      message("List of markers used : ", input$go.list.genes)
      message("Regulation selected : ", input$go.list.genes.reg)
      message("Nb of markers in it : ", length(List_markers[[input$go.list.genes]][[input$go.list.genes.reg]]))
    }
    
    if(input$go.list.genes.reg == "All"){
      
      cnetplot(x = en.go,
               showCategory = selected_terms,
               layout = selected_layout,
               cex.params = list(gene_node = 0.6, 
                                 gene_label = gene_label_size, 
                                 category_node = 1.6,
                                 category_label = term_label_size),
               color.params = list(category = "#2277cc",
                                   gene = "#552299",
                                   edge = T,
                                   foldChange = fc))+
        labs(color = "logFC")+
        theme(legend.box.margin = margin(l=0.3, unit = "in"),
              legend.title = element_text(size = 14, face = "bold", colour = "darkred", 
                                          margin = margin(t=0.3,b=0.1, unit = "in")),
              legend.text = element_text(size = 10, face = "bold"))+
        scale_color_gradientn(colours = c("midnightblue", "white", "darkred"),
                              values = c(0, 0.5, 1),
                              limits = c(-4, 4))
      
      
    } else if(input$go.list.genes.reg == "Up"){
      
      cnetplot(x = en.go,
               showCategory = selected_terms,
               layout = selected_layout,
               cex.params = list(gene_node = 0.6, 
                                 gene_label = gene_label_size, 
                                 category_node = 1.6,
                                 category_label = term_label_size),
               color.params = list(category = "#2277cc",
                                   gene = "#552299",
                                   edge = T,
                                   foldChange = fc))+
        labs(color = "logFC")+
        theme(legend.box.margin = margin(l=0.3, unit = "in"),
              legend.title = element_text(size = 14, face = "bold", colour = "darkred", 
                                          margin = margin(t=0.3,b=0.1, unit = "in")),
              legend.text = element_text(size = 10, face = "bold"))+
        scale_color_gradientn(colours = c("white", "gold", "darkred"),
                              values = c(0, 0.5, 1),
                              limits = c(0, 4)) 
      
      
    } else if(input$go.list.genes.reg == "Down"){
      
      cnetplot(x = en.go,
               showCategory = selected_terms,
               layout = selected_layout,
               cex.params = list(gene_node = 0.6, 
                                 gene_label = gene_label_size, 
                                 category_node = 1.6,
                                 category_label = term_label_size),
               color.params = list(category = "#2277cc",
                                   gene = "#552299",
                                   edge = T,
                                   foldChange = fc))+
        labs(color = "logFC")+
        theme(legend.box.margin = margin(l=0.3, unit = "in"),
              legend.title = element_text(size = 14, face = "bold", colour = "darkred", 
                                          margin = margin(t=0.3,b=0.1, unit = "in")),
              legend.text = element_text(size = 10, face = "bold"))+
        scale_color_gradientn(colours = c("midnightblue", "gold", "white"),
                              values = c(0, 0.5, 1),
                              limits = c(-4, 0))
      
      
    }
    
  })
  
  
  output$plt9 <- renderPlot({
    req(CNetPlot())
    CNetPlot()
    
  }, height = 800, width = 1000)
  
  
  ####  Corr plot   ------------------------------------------------------------
  
  en.go.paired <- eventReactive(en.go(), {
    req(en.go())
    en.go <- en.go()
    pairwise_termsim(en.go, method = "JC", showCategory = 250)
  })
  
  Cplor <- eventReactive(input$act14, {
    req(en.go.paired())
    en.go.paired <- en.go.paired()
    selected_terms <- input$terms_corrplot
    
    term_sim_mat <- en.go.paired@termsim[selected_terms,selected_terms]
    
    corrplot::corrplot(term_sim_mat,
                       type = "upper",
                       method = input$corr.methode,
                       tl.col = "black",
                       hclust.method = input$corr.hclust)
  })
  
  output$plt10 <- renderPlot({
    req(Cplor())
    Cplor()
  }, height = 1000, width = 1000)
  
  
  
  
  #### Tree plot   ---------------------------------------------------------------
  
  output$cnet_selected_terms <- renderPrint({
    req(input$terms_treeplot)
    
    selected_terms <- input$terms_treeplot
    
    cat(selected_terms, sep = "\n")
  })
  
  
  Tplot <- eventReactive(input$act15, {
    req(en.go.paired())
    en.go <- en.go.paired()
    
    treeplot(en.go, 
             showCategory = input$terms_treeplot,
             fontsize = 4,
             color = "p.adjust",
             cex_category = input$treeplot_cex_category,
             cluster.params = list(n= input$treeplot_nCluster,
                                   method = input$treeplot_Clust_method,
                                   color = NULL,
                                   label_words_n = 3,
                                   label_format = 10),
             offset.params = list(bar_tree = rel(2), tiplab = rel(4), extend = 0.5, hexpand = 0.15),
             highlight.params = list(align = "both"), hilight = input$treeplot_Clust_highlight)+
      scale_color_gradientn(colours = c("darkred","gold"),
                            values = c(0,1),
                            limits = c(0,0.05),
                            guide = guide_colorbar(title = "pValue.adj"))+
      scale_size_continuous(range = c(1,7))+
      labs(colour = "p.adjusted")+
      theme(legend.title = element_text(face = "bold", 
                                        colour = "darkred",
                                        size = 12,
                                        margin = margin(b=0.2, unit = "in")))
    
  })
  
  
  
  output$plt11 <- renderPlot({
    req(Tplot())
    Tplot()
  }, height = 1000, width = 850)
  
  
  
  
  
  
  
  ##  VI- Downloading   ------------------------------------------------------
  
  output$download <- downloadHandler(
    filename = function(){
      paste0(input$filename, ".xlsx")
    },
    content = function(file){
      res <- DEA_results() %>% 
        dplyr::filter(pct.1 >= input$pct.1_filt[1] & pct.1 <= input$pct.1_filt[2],
                      pct.2 >= input$pct.2_filt[1] & pct.2 <= input$pct.2_filt[2],
                      avg_log2FC >= input$logfc_filt[1] & avg_log2FC <= input$logfc_filt[2] |
                        avg_log2FC >= input$logfc_filt2[1] & avg_log2FC <= input$logfc_filt2[2],
                      abs(Diff_pct) >= input$diffpct_filt) %>%
        dplyr::arrange(desc(avg_log2FC))
      
      writexl::write_xlsx(res, path = file)
    }
  )
  
  
  output$savebarplot <- downloadHandler(
    filename = function() {
      paste0(input$barplot_name, ".pdf")
    },
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
        theme(
          axis.title = element_text(face= "bold", size= 15),
          axis.text.y = element_text(size= 15, face = "bold"),
          axis.text.x = element_text(size= 13, face = "bold"),
          legend.title = element_text(face= "bold", size= 15, colour= "navy"),
          legend.text = element_text(size= 13, colour= "black"),
          plot.title = element_text(size = 17, colour = "darkred", face = "bold")
        ) -> Plot
      
      print(Plot)
      
      dev.off()
    }
  )
  
  
  output$savedotplot <- downloadHandler(
    filename = function() {
      paste0(input$dotplot_name, ".pdf")
    },
    content = function(file){
      pdf(file, width = input$width_dplot, height = input$height_dplot)
      
      seurat <- seurat()
      genes <- if (input$radio.select == "Up Genes") {
        Up_Genes()[input$slider.nbGenes[1]:input$slider.nbGenes[2]]
      } else {
        Down_Genes()[input$slider.nbGenes[1]:input$slider.nbGenes[2]]
      }
      
      DotPlot(seurat,
              features = genes,
              group.by = input$select.meta)+
        RotatedAxis()+
        coord_flip()+
        theme_bw()+
        labs(x="",
             y="")+
        scale_color_gradient(low = "orange", high = "navy")+
        theme(
          axis.title = element_text(face= "bold", size= 15),
          axis.text.y = element_text(size= 16, face = "bold"),
          axis.text.x = element_text(size= 16, face = "bold"),
          legend.title = element_text(face= "bold", size= 15, colour= "navy"),
          legend.text = element_text(size= 13, colour= "black"),
          plot.title = element_text(size = 17, colour = "darkred", face = "bold")
        ) -> Plot2
      
      print(Plot2)
      
      dev.off()
    }
  )
  
  
  output$savehm <- downloadHandler(
    filename = function() {
      paste0(input$hm_name, ".pdf")
    },
    content = function(file){
      pdf(file, width = input$width_hm, height = input$height_hm)
      
      hm <- hm()
      
      ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom")
      
      dev.off()
    }
  )
  
  
  output$download_all_markers <- downloadHandler(
    filename = function() {
      paste("List_markers", Sys.Date(), ".rds", sep = "")
    },
    content = function(file) {
      saveRDS(List_markers(), file)
    }
  )
  
}











#  Run App   ############################

shinyApp(ui = ui,
         server = server)
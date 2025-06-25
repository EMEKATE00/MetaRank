################################################################################
#                                                                              #
#                               Master’s Thesis                                #
#                                                                              #
#     Title: TFM: Development of the Meta-Analysis Tool MetaRank               #
#                                                                              #
#     Author: Maksym Kupchyk Tiurin                                            #
#     Supervisor: Francisco García García                 Year: 2025           #
#     Co-supervisor: Rubén Grillo Risco                                        #
#                                                                              #
#     Academic Institution: CIPF - Centro de Investigación Príncipe Felipe     #
#                                                                              #
################################################################################

#-------------------------------------------------------------------------------
# Loading of the required dependencies 
#-------------------------------------------------------------------------------

library(shiny)                  # enables building interactive web app
library(shinycssloaders)        # add a loading animation to outputs instead
library(shinyWidgets)           # custom input controls and user interface components
library(shinyBS)                # adds additional Twitter Bootstrap components 
library(shinyjs)                # perform common useful JavaScript operations
library(bslib)                  # toolkit for Shiny and R based on Bootstrap
library(dplyr)                  # tool for working with data.frame objects
library(DT)                     # interactive representation of data.frames
library(ggplot2)                # system for declaratively creating graphics
library(processx)               # run system processes in the background
library(plotly)                 # makes interactive, publication-quality graphs
library(zip)                    # allows dwnload more than one file at the same time

#------------------------------------------------------------------------------- 
# The user interface (ui) object controls the layout and appearance of your app
#-------------------------------------------------------------------------------

ui <- page_sidebar(
  
  # Add a CSS file to change aesthetics and layout
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "style.css")),
  
  tags$head(
    tags$style(HTML("
    .dt-header-with-icon { display: inline-flex; align-items: center; }
    .dt-header-label { margin-right: 4px; }
    .dt-info-icon {
      display: inline-block; width: 16px; height: 16px; line-height: 16px;
      text-align: center; border-radius: 50%; background: #888; color: white;
      font-size: 12px; cursor: help;
    }
    .dt-info-icon:hover { background: #555; }
  "))),
  
  
  title = "MetaRank",
  
  sidebar = sidebar(width = "380px", useShinyjs(),
                    
            # Select the type of analysis: 
            div(selectInput(
                inputId = "selectedPackage",
                label   = "Select Meta-Analysis type",
                choices = c(
                  "Weighted (ranked by scores)"       = "Rankprodpack",
                  "Unweighted (simple ranked lists)"   = "RRApack"),
                selected = "Rankprodpack"),
              title = HTML(
                "Choose between:\n- Weighted: Gene lists are ranked based on numerical scores that represent their importance or association with a given condition. These scores could be p-values, fold-changes, or other statistical measures that quantify the relevance of each gene.\n- Unweighted: Gene lists are simply ranked by their order in different datasets, without assigning any numerical score or weight. This method assumes that the genes’ relative positions across datasets are the key factor, rather than a statistical score.\n- The Weighted method is more suitable when you have numerical data that reflects the significance of each gene. The Unweighted method is ideal for cases where you are only interested in the relative rankings of genes across different datasets."),
              style = "cursor: help; display: inline-block; width: 100%;"),

            
            tabsetPanel(
              tabPanel(title = "Analysis", br(),
                
                #--------------------------------------------------------------
                # Data input methods
                #--------------------------------------------------------------
                
                tags$h5("Input"),
                
                radioButtons(inputId  = "inputMethod",
                             label    = "Input Method",
                             choices  = c("Upload Files" = "files", "Paste Genes" = "paste"),
                             selected = "files"),
                
                
                # Conditional: Files + Rankprodpack
                #--------------------------------------------------------------
                conditionalPanel(condition = "input.inputMethod == 'files' && input.selectedPackage == 'Rankprodpack'",
                  fluidRow(
                    column(10, fileInput(inputId   = "files",
                                         label     = "Upload gene rankings (RP)",
                                         multiple  = TRUE,
                                         accept    = c(".txt", ".tsv", ".csv"))),
                    column(2, actionButton(inputId = "showExample", label = "", 
                                          icon     = icon("info-circle", lib = "font-awesome"),
                                          class    = "btn-info"))),
                  
                  tags$p("Use example data?"),
                  fluidRow(
                    column(10,
                      div(
                        switchInput(
                          inputId  = "useExampleData",
                          onLabel  = "Yes",
                          offLabel = "No",
                          value    = FALSE),
                        title = HTML(
                          "These example human gene lists (4 studies related to lung cancer: GSE10072, GSE19188, GSE63459, GSE75037) contain gene symbols and p-values, have no headers, and include duplicates and blank lines to simulate raw data."),
                        style = "cursor: help;")),
                    column(2,
                      downloadLink(
                        outputId = "downloadExampleRank",
                        label    = icon("download"),
                        style    = "margin-top: 10px;")))),
                
                
                # Condicional: Files + RRApack
                #--------------------------------------------------------------
                conditionalPanel(
                  condition = "input.inputMethod == 'files' && input.selectedPackage == 'RRApack'",
                  fluidRow(
                    column(10, fileInput(inputId  = "files_rra",
                                         label    = "Upload gene rankings (RRA)",
                                         multiple = TRUE,
                                         accept   = c(".txt"))),
                    column(2, actionButton(inputId = "showExampleRRA", label = "", 
                                          icon    = icon("info-circle", lib = "font-awesome"),
                                          class   = "btn-info"))),
                  
                  tags$p("Use example data?"),
                  fluidRow(
                    column(10,
                      div(
                        switchInput(
                          inputId  = "useExampleData_rra",
                          onLabel  = "Yes",
                          offLabel = "No",
                          value    = FALSE),
                        title = HTML(
                          "These example human gene lists (4 studies related to lung cancer: GSE10072, GSE19188, GSE63459, GSE75037) contain gene symbols only (no p-values), have no headers, include duplicates and blank lines to simulate raw data."),
                        style = "cursor: help;")),
                    column(2,
                      downloadLink(
                        outputId = "downloadExampleRRA",
                        label    = icon("download"),
                        style    = "margin-top: 10px;")))),
                
                
                # Condicional: Paste + Rankprodpack
                #--------------------------------------------------------------
                conditionalPanel(
                  condition = "input.inputMethod == 'paste' && input.selectedPackage == 'Rankprodpack'",
                  radioButtons(inputId  = "pasteFormat",
                               label    = "Select Paste Format",
                               choices  = c("CSV" = "csv", "TSV" = "tsv"),
                               selected = "csv"),
                  
                  conditionalPanel(
                    condition = "input.pasteFormat == 'csv'",
                    div(
                      textAreaInput(
                        inputId     = "pasteText",
                        label       = "Paste your data",
                        width       = "300px",
                        height      = "200px",
                        placeholder = "Gene,Stat.data\nMyc,1.2e-06\nGapdh,3.5e-04\nIl2ra,0.0021\n...\n###\nGene,Stat.data\nGapdh,5.4e-07\nIl6,9.8e-04\nNanog,0.0062\n..."),
                      title = HTML("CSV Paste Format: \n- Each block must start with Gene,Stat.data header. \n- Use comma ',' as field separator and '.' as decimal point.\n- Separate multiple lists with '###' on its own line.\n- Do not use commas as decimal separators.\n- For full details, see the format instructions."),
                      style = "cursor: help;")),
                  
                  conditionalPanel(
                    condition = "input.pasteFormat == 'tsv'",
                    div(
                      textAreaInput(
                        inputId     = "pasteText1",
                        label       = "Paste your data",
                        width       = "300px",
                        height      = "200px",
                        placeholder = "Gene\tStat.data\nMyc\t1.2e-06\nGapdh\t3.5e-04\nIl2ra\t0.0021\n...\n###\nGene\tStat.data\nGapdh\t5.4e-07\nIl6\t9.8e-04\nNanog\t0.0062\n..."),
                      title = HTML(
                        "TSV Paste Format:\n- Each block must start with 'Gene\\tStat.data' header.\n- Use tabulation '\\t' as field separator and '.' as decimal point.\n- Separate multiple lists with '###' on its own line.\n- Do not use commas or semicolons as separators.\n- For full details, see the format instructions."),
                      style = "cursor: help;"))),
                
                
                # Condicional: Paste + RRApack
                #--------------------------------------------------------------
                conditionalPanel(
                  condition = "input.inputMethod == 'paste' && input.selectedPackage == 'RRApack'",
                  div(
                    textAreaInput(
                      inputId     = "pasteText_rra",
                      label       = "Paste your ranked gene lists",
                      width       = "300px",
                      height      = "200px",
                      placeholder = "TP53\nMYC\nEGFR\nBRCA1\n...\n###\nBRCA1\nPTEN\nMYC\n..."),
                    title = HTML(
                      "RRA Paste Format:\n- No header required; only gene symbols.\n- One gene per line, separated by '\\n'.\n- Separate multiple lists with '###' on its own line.\n- Do not include numeric values or additional columns.\n- Raw, unweighted gene rankings only."),
                    style = "cursor: help;")),
                
                
                tags$hr(style = "margin-top: 40px; margin-bottom: 40px;"),
                
                tags$h5("Parameters"),
                
                # ------------------------------------------------------------
                # Conditional: Rankprodpack
                # ------------------------------------------------------------
                conditionalPanel(
                  condition = "input.selectedPackage == 'Rankprodpack'",
                  
                  div(
                    selectInput(
                      inputId = "RankBasedMethod",
                      label   = "Rank-Based Meta-Analysis Methods",
                      choices = c("RankProd basic" = "rp", "RankProd advance" = "rp.adv"),
                      selected = "rp"),
                    title = HTML(
                      "Choose your RankProd method:\n- RankProd basic: Basic analysis without additional grouping.\n- RankProd advance: Advanced RankProd allowing an 'origin' vector to group your lists by metadata (e.g., data source, date of acquisition, technology)."),
                    style = "cursor: help; display: inline-block;"),
                  
                  conditionalPanel(
                    condition = "input.RankBasedMethod == 'rp.adv'",
                    fluidRow(
                      column(10, textInput(
                                 inputId     = "origin",
                                 label       = "Origin",
                                 placeholder = "1,1,2,2,1")),
                      column(2,  actionButton(
                                 inputId     = "showOriginInfo",
                                 label       = "",
                                 icon        = icon("info-circle", lib = "font-awesome"),
                                 class       = "btn-info")))),
                  
                  # Minimum number of datasets
                  div(sliderInput(
                      inputId = "gene_app_filter_rp",
                      label   = "Minimum number of datasets",
                      min     = 1,
                      max     = 4,
                      value   = 4,
                      step    = 1,
                      ticks   = TRUE),
                    title = HTML(
                      "Gene Appearance Filter:\n - Require a gene to appear in at least this many lists.\n- E.g., if BRCA1 appears in 2/4 lists and you set this to 3, it will be excluded.\n- Setting to 2 will keep it."),
                    style = "cursor: help;"),
                  
                  # Ranking direction
                  div(radioButtons(
                      inputId  = "ranking_direction",
                      label    = "Ranking Direction",
                      choices  = c(
                        "Ascending (lower values = better)"  = "asc",
                        "Descending (higher values = better)" = "desc"),
                      selected = "asc",
                      inline   = TRUE),
                    title = HTML(
                      "Ranking Direction:\n- Ascending: Best for p-values (smaller = more significant).\n- Descending: Best for fold-changes (larger = more significant).\n- Preserves the order in your input lists."),
                    style = "cursor: help;"),
                  
                  # NA management
                  div(selectInput(
                      inputId  = "NA_management",
                      label    = "NA Management",
                      choices  = c(
                        "Impute NA"   = "impute",
                        "Ignore NA"   = "ignore",
                        "Penalize NA" = "penalize"),
                      selected = "impute"),
                    title = HTML(
                      "Missing Value Handling:\n- Impute NA: Replace missing with the gene’s median across lists (extra penalty optional).\n- Ignore NA: Analyze only present values (extra penalty optional).\n- Penalize NA: Replace missing with the worst value (no extra penalty)."),
                    style = "cursor: help;"),
                  
                  # Extra penalization (only when Impute or Ignore is selected)
                  conditionalPanel(
                    condition = "input.NA_management == 'impute' || input.NA_management == 'ignore'",
                    tags$p("Apply an extra penalization?"),
                    div(switchInput(
                        inputId  = "extra_penalty_switch",
                        onLabel  = "Yes",
                        offLabel = "No",
                        value    = FALSE),
                      title = HTML(
                        "Extra Penalization:\n- After initial ranking, adds a penalty based on missing appearances.\n- AdjustedRank = Rank + ((total_lists − Count) × (max_rank / total_lists)).\n- Genes in fewer lists receive a larger penalty."),
                      style = "cursor: help;"))),
                
                
                # -----------------------------
                # Conditional: RRApack
                # -----------------------------
                conditionalPanel(
                  condition = "input.selectedPackage == 'RRApack'",
                  
                  div(selectInput(
                      inputId = "rra_method",
                      label   = "Aggregation Method",
                      choices = c(
                        "Robust Rank Aggregation" = "RRA",
                        "Min"                      = "min",
                        "Geometric Mean"           = "geom.mean",
                        "Arithmetic Mean"          = "mean",
                        "Median"                   = "median",
                        "Stuart"                   = "stuart"),
                      selected = "RRA"),
                    title = HTML(
                      "Aggregation Methods:\n- RRA: Probabilistic rank aggregation robust to noise.\n- Min: Takes the minimum rank across lists.\n- Geometric Mean: Geometric mean of ranks.\n- Arithmetic Mean: Arithmetic mean of ranks.\n- Median: Median of ranks.\n- Stuart: Stuart’s method for order-statistics aggregation."),
                    style = "cursor: help;"),
                  
                  div(sliderInput(
                      inputId = "gene_app_filter_rra",
                      label   = "Minimum number of datasets",
                      min     = 1,
                      max     = 4,
                      value   = 4,
                      step    = 1,
                      ticks   = TRUE),
                    title = HTML(
                      "Gene Appearance Filter:\n- Require a gene to appear in at least this many lists.\n- E.g., setting to 2 keeps genes in ≥2/4 lists."),
                    style = "cursor: help;"),
                  
                  div(selectInput(
                      inputId = "rra_na_handling",
                      label   = "Missing Value Handling",
                      choices = c(
                        "Ignore missing values"   = "ignore",
                        "Penalize missing values" = "penalize"),
                      selected = "ignore"),
                    title = HTML(
                      "Missing Value Handling:\n- Ignore: Proceed with only present values.\n- Penalize: Substitute missing with the worst observed rank."),
                    style = "cursor: help;")),
                
                
                tags$hr(style = "margin-top: 40px; margin-bottom: 40px;"),
                
                #--------------------------------------------------------------
                # Analysis options
                #--------------------------------------------------------------
                
                tags$h5("Enrichment analysis"),
                
                div(checkboxInput(inputId = "enableEnrichment",
                              label = "Perform Enrichment Analysis?",
                              value = FALSE),
                title = HTML("Check this box to run an Over-Representation Analysis (ORA) on your meta-analysis results. It will take the top N genes you’ve selected (up to 100) from the consensus ranking and test for pathway or GO term enrichment."),
                style = "cursor: help;"),
                
                conditionalPanel(
                  condition = "input.enableEnrichment == true",
                      
                      div(sliderInput(inputId = "topGenes",
                                  label = "Number of Top Genes",
                                  min = 10,
                                  max = 100,
                                  value = 20,
                                  step = 1),
                      title = HTML("Select how many of the top-ranked genes (between 10 and 100) from the meta-analysis consensus list you want to use for the over-representation analysis (ORA)."), 
                      style = "cursor: help;"),
                  
                      selectInput(inputId = "oraDatabase",
                                  label = "ORA Database",
                                  choices = c("Gene Ontology (GO)" = "go", "KEGG" = "kegg", "Reactome" = "reactome"),
                                  selected = "go"),
                      
                      conditionalPanel(
                        condition = "input.oraDatabase == 'go'",
                        
                        radioButtons(inputId = "ontology",
                                     label = "Ontology",
                                     choices = c("Biological Process (BP)" = "bp",
                                                 "Molecular Function (MF)" = "mf",
                                                 "Cellular Component (CC)" = "cc"),
                                     selected = "bp")),
                  
                      selectInput(inputId = "species",
                                  label = "Species",
                                  choices = c("Homo sapiens" = "Hsa",
                                              "Mus musculus" = "Mmu",
                                              "Rattus norvegicus" = "Rno"),
                                  selected = "Hsa"),
                  
                      selectInput(inputId = "IDtype",
                                  label = "Gene ID",
                                  choices = c("SYMBOL", "ENTREZID", "ENSEMBL"),
                                  selected = "SYMBOL"))),
        
              
              
        tabPanel(
          title = "Data Visualization", br(),
          
          #--------------------------------------------------------------------
          # Meta-analysis table setting
          #--------------------------------------------------------------------
          
          tags$h5("Meta-analysis table setting"),
          
          conditionalPanel(
            condition = "input.selectedPackage == 'Rankprodpack'",
            checkboxGroupInput(
              inputId = "selectedColumnsResultsRP",
              label = "Select Columns to Display",
              choices = c("GeneID", "Rank", "FileCount", "FileNames", "GenePositions", 
                          "RP_stat", "PFP", "pvalue", "p.adjust"),
              selected = c("GeneID", "Rank", "FileCount", "FileNames", "GenePositions", 
                           "RP_stat", "PFP", "pvalue", "p.adjust"),
              inline = FALSE)), 
          
          conditionalPanel(
            condition = "input.selectedPackage == 'RRApack'",
            checkboxGroupInput(
              inputId = "selectedColumnsResultsRRA",
              label = "Select Columns to Display",
              choices = c("GeneID", "Rank", "Score", "p.adjust", "FileCount", "FileNames", "GenePositions"),
              selected = c("GeneID", "Rank", "Score", "p.adjust", "FileCount", "FileNames", "GenePositions"),
              inline = FALSE)), 
          
          tags$hr(style = "margin-top: 40px; margin-bottom: 40px;"),
          
          # -------------------------------------------------------------------
          # Upset Plot Settings
          # -------------------------------------------------------------------

          tags$h5("Upset Plot Settings"),
          
          sliderInput(inputId = "upsetTextSize",
                      label = "Upset Text Size",
                      min = 15, max = 40, value = 25),
          colorPickr(inputId = "upsetSetsColor",
                     label = "Sets Bar Color (horizontal)",
                     selected = "#2ba915",
                     opacity = TRUE),
          colorPickr(inputId = "upsetIntersectionColor",
                     label = "Intersection Bar Color (vertical)",
                     selected = "#0838a0",
                     opacity = TRUE), 
          
          tags$hr(style = "margin-top: 40px; margin-bottom: 40px;"),
          
          # -------------------------------------------------------------------
          # Heatmap Settings
          # -------------------------------------------------------------------
          
          tags$h5("Heatmap Settings"),
          
          sliderInput(inputId = "heatmapTickSize",
                      label = "Heatmap Tick Size",
                      min = 5, max = 30, value = 15),
          sliderInput(inputId = "heatmapTitleSize",
                      label = "Heatmap Title Size",
                      min = 10, max = 30, value = 20),
          selectInput(inputId = "heatmapColorscale",
                      label = "Heatmap Colorscale",
                      choices = c("Viridis", "Cividis", "Portland"),
                      selected = "Viridis"),
          
          tags$hr(style = "margin-top: 40px; margin-bottom: 40px;"),
          
          # -------------------------------------------------------------------
          # Enrichment Plot Settings (Dot/Bar Plot)
          # -------------------------------------------------------------------
          
          tags$h5("Enrichment table Settings"),
          
          checkboxGroupInput(
            inputId = "selectedColumnsRank",
            label = "Select Columns to Display in the Table",
            choices = c("ID", "Description", "GeneRatio", "BgRatio", 
                        "pvalue", "p.adjust", "Count", "geneID"),
            selected = c("ID", "Description", "GeneRatio", "BgRatio", 
                         "pvalue", "p.adjust", "Count", "geneID"),
            inline = FALSE),
          
          tags$hr(style = "margin-top: 40px; margin-bottom: 40px;"),
          
          tags$h5("Enrichment Plot Settings"),
          
          sliderInput(inputId = "nPlotTerms",
                      label = "Number of Terms to Show",
                      min = 1, max = 20, value = 10),
          
          selectInput(inputId = "plotLabel",
                      label = "Y-Axis Label",
                      choices = c("Term" = "ID", 
                                  "Description" = "Description", 
                                  "Term + Description" = "ID_Description"),
                      selected = "ID_Description"),
          
          selectInput(inputId = "plotType",
                      label = "Plot Type",
                      choices = c("Dot Plot" = "dotplot",
                                  "Bar Plot" = "barplot"),
                      selected = "dotplot"),
          
          colorPickr(inputId = "colorLow",
                     label = "Low Color",
                     selected = "#0838a0",
                     opacity = TRUE),
          
          colorPickr(inputId = "colorHigh",
                     label = "High Color",
                     selected = "#2ba915",
                     opacity = TRUE),
          
          sliderInput(inputId = "textSize",
                      label = "Text Size",
                      min = 5, max = 20, value = 12))),
        
        actionButton(inputId = "Run", label = "Run Analysis", style = "width: 100%; margin-top: 15px;")),
    
  
    # -------------------------------------------------------------------------
    # DISPLAY RESULTS IN THE FORM OF GRAPHS AND TABLES 
    # -------------------------------------------------------------------------

    card(tabsetPanel(id = "mainTabs",
        
        tabPanel("MetaRanking Results", br(),
          
          shinycssloaders::withSpinner(
            dataTableOutput("Analysis_results_table")), br(), 
          
          uiOutput("downloadButtonUI"), br(),
          
          shinyjs::hidden(
            div(id = "excludedGenesSection", br(),
                
                h5("Excluded terms (below frequency threshold):"),
                dataTableOutput("excludedGenesTable"))), br(),
          
          shinycssloaders::withSpinner(plotOutput("upsetPlot", height = "600px")), br(),
          
          conditionalPanel(
            condition = "output.Hidebuttons",
            uiOutput("downloadButtonUpset")), br(),
          
          shinycssloaders::withSpinner(plotlyOutput("heatmapPlot", height = "800px")), br(),
          
          uiOutput("downloadHeatmapUI")),
        
        tabPanel("Enrichment Results", br(), 
          
          shinycssloaders::withSpinner(DT::dataTableOutput("enrichmentTable")), br(),
          
          uiOutput("downloadButtonUI2"), br(),
          
          shinycssloaders::withSpinner(plotlyOutput("enrichmentDotPlot", height = "600px")),
          
          uiOutput("downloadPlotUI3")))))




server <- function(input, output, session) {
  
  # NO NEED, BUT CHARGE THE FUNCTIONS
  source("./R/MetaRank_Functions.R")
  source("./R/ORA.R")

  # CHARGE THE FILES USED AS EXAMPLE
  example_files <- c(
    "./example_data/data1.tsv",
    "./example_data/data2.tsv",
    "./example_data/data3.tsv",
    "./example_data/data4.tsv")
  
  example_files_rra <- c(
    "./example_data/gene1.txt",
    "./example_data/gene2.txt",
    "./example_data/gene3.txt",
    "./example_data/gene4.txt")
  
  # SET ANNOTATION DIR WITH ALL THE FILES
  annotation_dir <- file.path(getwd(), "database_annotations")
  
  
  #----------------------------------------------------------------------------
  # DATA READING
  #----------------------------------------------------------------------------
  gene_data <- eventReactive(input$Run, {
    
    withProgress(message = "Loading gene data...", value = 0.1, {
      
      # -------------------- RRA --------------------
      
      if (input$selectedPackage == "RRApack") {
        
        if (isTRUE(input$useExampleData_rra)) {
          lists <- read_rra_files(example_files_rra)
          
        } else if (input$inputMethod == "paste") {
          if (!validate_rra_format(input, session)) return(NULL)
          req(input$pasteText_rra)
          lists <- paste_rra_lists(input$pasteText_rra)
          
        } else if (input$inputMethod == "files") {
          if (!validate_file_format(input, session)) return(NULL)
          req(input$files_rra)
          lists <- read_rra_files(input$files_rra$datapath)
          names(lists) <- basename(input$files_rra$name)}
        
        return(lists)}
      
      
      # ---------------- RankProd con método "paste" ----------------
      
      if (input$selectedPackage == "Rankprodpack" &&
          input$inputMethod     == "paste") {
        
        text_input <- if (input$pasteFormat == "csv") input$pasteText else input$pasteText1
        
        if (is.null(text_input) || !nzchar(text_input)) {
          showModal(modalDialog(
            title = "⚠️ No data pasted",
            "Please paste your gene lists before running the analysis.",
            easyClose = TRUE,
            footer = modalButton("Close")))
          return(NULL)}
        
        is_csv <- grepl(",",  text_input)
        is_tsv <- grepl("\t", text_input)
        format_ok <- (input$pasteFormat == "csv" && is_csv) ||
          (input$pasteFormat == "tsv" && is_tsv)
        
        if (!format_ok) {
          showModal(modalDialog(
            title = tags$div("⚠️ Format Mismatch Detected", style = "text-align:center;font-weight:bold;"),
            easyClose = TRUE,
            size = "l",
            fluidRow(
              column(6,
                     tags$h5("CSV Example"),
                     HTML("<pre>Gene,Stat.data\nBRCA1,1.5\nTP53,2.3\nMYC,4.2\n...</pre>")),
              column(6,
                     tags$h5("TSV Example"),
                     HTML("<pre>Gene\tStat.data\nBrca1\t1.5\nTp53\t2.3\nMYC\t4.2\n...</pre>"))),
            tags$hr(),
            tags$h4("Recommendations"),
            tags$ul(
              tags$li("Select correct paste format (CSV or TSV)."),
              tags$li("Use comma (,) for CSV; tab (\\t) for TSV."),
              tags$li("Include a header row: Gene,Stat.data"),
              tags$li("Use '###' only if splitting multiple lists.")),
            footer = modalButton("Close")))
          return(NULL)}}
      
      
      # ---------------- RankProd input (example, paste, or files) ----------------
      
      if (isTRUE(input$useExampleData)) {
        dfs <- read_rankprod_files(example_files)
        
      } else if (input$inputMethod == "paste") {
        if (input$pasteFormat == "csv") {
          req(input$pasteText)
          dfs <- paste_rankprod_lists(input$pasteText, format = "csv")
          
        } else if (input$pasteFormat == "tsv") {
          req(input$pasteText1)
          dfs <- paste_rankprod_lists(input$pasteText1, format = "tsv")}
        
      } else if (input$inputMethod == "files") {
        if (!validate_file_format(input, session)) return(NULL)
        req(input$files)
        dfs <- read_rankprod_files(input$files$datapath)
        names(dfs) <- basename(input$files$name)}
      
      return(dfs)})})
  
  
  #----------------------------------------------------------------------------
  # DATA PREPARING
  #----------------------------------------------------------------------------
  
  prepared_data <- reactive({
    raw <- gene_data()
    req(raw)
    
    is_rra <- input$selectedPackage == "RRApack"
    is_valid <- if (is_rra) function(x) length(x) > 0 else function(x) is.data.frame(x) && nrow(x) > 0
    data_list <- Filter(is_valid, raw)
    
    req(length(data_list) > 0, "No valid input data.")
    names_list <- names(data_list)
    
    all_genes <- if (is_rra) data_list else lapply(data_list, function(df) if ("Gene" %in% colnames(df)) df$Gene else character(0))
    
    all_genes <- Filter(function(x) length(x) > 0, all_genes)
    req(length(all_genes) > 0, "No valid genes after filtering.")
    
    app_df <- count_gene_appearance_rra(all_genes, names_list)
    
    if (nrow(app_df) == 0) {
      showModal(modalDialog(
        title = "⚠️ No Repeated Genes Found",
        "There are no genes that appear in more than one list. Please check your input data or lower the appearance threshold.",
        easyClose = TRUE,
        footer = modalButton("Close")))
      return(NULL)}
    
    app_df <- app_df[app_df$Count > 1, , drop = FALSE]
    
    max_appearance <- if (nrow(app_df) > 0) max(app_df$Count) else 0
    
    min_appearance <- if (input$selectedPackage == "Rankprodpack") {
      req(input$gene_app_filter_rp)
      input$gene_app_filter_rp
    } else {
      req(input$gene_app_filter_rra)
      input$gene_app_filter_rra}
    
    if (min_appearance > max_appearance && max_appearance > 0) {
      showModal(modalDialog(
        title = "⚠️ Appearance Threshold Too High",
        paste0("You have set a minimum appearance filter of ", min_appearance, 
               ", but the maximum gene appearance is ", max_appearance, ". ",
               "Please adjust the minimum appearance threshold to be less than or equal to ", max_appearance, 
               " and try again."),
        easyClose = TRUE,
        footer = modalButton("Close")))
      return(NULL)
    } else if (max_appearance == 0) {
      showModal(modalDialog(
        title = "⚠️ No Genes Meet the Appearance Threshold",
        "No genes were found that appear in more than one list. Please check your input or lower the minimum appearance threshold.",
        easyClose = TRUE,
        footer = modalButton("Close")))
      return(NULL)}
    
    genes_keep <- app_df$Gene[app_df$Count >= min_appearance]
    
    valid <- if (is_rra) {
      lapply(data_list, function(lst) intersect(lst, genes_keep))
    } else {
      lapply(data_list, function(df) if ("Gene" %in% colnames(df)) filter(df, Gene %in% genes_keep) else df)}
    
    excluded <- if (is_rra) {
      lapply(data_list, function(lst) data.frame(Gene = setdiff(lst, genes_keep)))
    } else {
      lapply(data_list, function(df) if ("Gene" %in% colnames(df)) filter(df, !Gene %in% genes_keep) else df)}
    
    list(valid = valid, excluded = excluded, files = names_list)})
  
  
  
  
  #----------------------------------------------------------------------------
  # PERFORM THE META-ANALYSIS
  #----------------------------------------------------------------------------
  
  combined_results <- eventReactive(input$Run, {
    raw <- gene_data()
    data <- prepared_data()
    req(data$valid)
    
    # ORIGIN FORMAT VALIDATION-------------------------------------------------
    if (!validate_origin(input, session, data)) {
      return(NULL)}
    
    withProgress(message = "Running analysis...", value = 0, {
      if (input$selectedPackage == "Rankprodpack") {
        direction <- if (input$ranking_direction == "asc") "ascending" else "descending"
        na_rm <- input$NA_management != "ignore"
        
        statdata <- generate_statdata_matrix(
          gene_data = data$valid,
          na_management = input$NA_management,
          ranking_direction = direction)
        
        cl <- rep(1, ncol(statdata))
        res <- NULL
        
        if (input$RankBasedMethod == "rp") {
          res <- perform_RP(statdata, cl,
                            logged = TRUE,
                            na.rm = na_rm,
                            ranking_direction = direction,
                            progress = incProgress)
        } else {
          origin <- as.numeric(strsplit(gsub("\\s+", "", input$origin), ",")[[1]])
          req(length(origin) == ncol(statdata))
          req(all(!is.na(origin)))
          
          res <- perform_RP_advance(statdata, cl, origin,
                                    logged = TRUE,
                                    na.rm = na_rm,
                                    ranking_direction = direction,
                                    progress = incProgress)}
        
        app <- count_gene_appearance_with_files(data$valid, data$files)
        
        apply_penalty <- input$extra_penalty_switch && input$NA_management %in% c("impute", "ignore")
        tbl <- if (apply_penalty) {
          create_combined_rank_table2(res, app)
        } else {
          create_combined_rank_table(res, app)}
        
      } else {
        full_flag <- input$rra_na_handling == "ignore"
        tbl <- run_rra_analysis(data$valid, method = input$rra_method, full = full_flag)}
      
      pos_df <- get_gene_positions(raw)
      tbl <- left_join(tbl, pos_df, by = "GeneID")
      
      incProgress(1, detail = "Analysis complete.")
      tbl})})
  
  #----------------------------------------------------------------------------
  # DATA SELECT FOR THE ENRICHMENT ANALYSIS
  #----------------------------------------------------------------------------
  
  selected_genes <- eventReactive(input$Run, {
    req(combined_results(), input$enableEnrichment)
    withProgress(message = "Selecting the top genes...", detail = "Remember you can select how much genes you want to analyze", value = 0, {
      result <- get_top_genes(combined_results(), input$topGenes)
      incProgress(1, detail = "We got the most relevant genes!")
      result})})
  
  
  #----------------------------------------------------------------------------
  # ENRICHMENT ANALYSIS OF THE FIRST 100 GENES
  #----------------------------------------------------------------------------
  
  enrichment_results <- reactiveVal(NULL)
  
  observeEvent(input$Run, {
    req(input$species, input$oraDatabase, input$ontology, selected_genes)
    
    if (input$enableEnrichment) {
      genes_to_test <- selected_genes()
      
      # GENEID VALIDATION------------------------------------------------------
      valid <- validate_gene_ids(selected_genes(), input$IDtype)
      if (!valid) return()
      
      # ORGANISM VALIDATION----------------------------------------------------
      if (!validate_organism(genes_to_test, input$species, input$IDtype, session)) {
        return()}
      
      withProgress(message = "Ejecutando análisis de enriquecimiento", value = 0.1, {
        
        result <- switch(input$oraDatabase,
                         "go" = ORA_GO(
                           gene_list = selected_genes(),
                           organism  = input$species,
                           ontology  = toupper(input$ontology),
                           pvalue    = 1,
                           ID_type   = input$IDtype),
                         
                         "kegg" = ORA_KEGG(
                           gene_list = selected_genes(),
                           organism  = input$species,
                           pvalue    = 1,
                           ID_type   = input$IDtype),
                         
                         "reactome" = ORA_REACTOME(
                           gene_list = selected_genes(),
                           organism  = input$species,
                           pvalue    = 1,
                           ID_type   = input$IDtype),
                         
                         stop("Base de datos seleccionada no soportada."))
        
        if (!check_empty_results(result, session)) {
          return()}
        
        enrichment_results(result)})
    } else {
      enrichment_results(NULL)}})
  
  
  #----------------------------------------------------------------------------
  # DATA REPRESENTATION OF META-ANALYSIS (CARD)
  #----------------------------------------------------------------------------
  
  # RANKPROD ANALYSIS (TABLE + 2 DOWNLOAD BUTTONS + 2 BUTTONS OF EXCLUDED GENES)
  #----------------------------------------------------------------------------
  output$Analysis_results_table <- DT::renderDataTable({
    req(combined_results()) 
    req(input$selectedPackage) 
    
    selected_columns <- if (input$selectedPackage == "Rankprodpack") {
      input$selectedColumnsResultsRP
    } else if (input$selectedPackage == "RRApack") {
      input$selectedColumnsResultsRRA
    } else {
      return(data.frame(Message = "No package selected"))}
    
    data <- combined_results()
    if (nrow(data) == 0) return(data.frame(Message = "No data available"))
    
    selected_data <- data[, selected_columns, drop = FALSE]
    if (nrow(selected_data) == 0) return(data.frame(Message = "No data available"))
    
    num_rows <- input$Analysis_results_table_length
    table_height <- paste0(num_rows * 40, "px")
    
    colTooltips4 <- list(
      GeneID = "Unique gene identifier, which can be a HUGO symbol (e.g., TP53), an Entrez ID (e.g., 7157), or an Ensembl ID (e.g., ENSG00000141510), depending on the selected input format.",
      Rank = "Consensus ranking of the gene across all input lists; lower values indicate higher overall relevance or consistency among the lists.",
      FileCount = "Number of input gene lists in which this gene appears; a higher count suggests greater consistency across datasets.",
      FileNames = "Names of the input files where the gene was found, separated by spaces; useful for identifying the sources supporting the gene's relevance.",
      RP_stat = "Rank Product statistic calculated to assess the significance of the gene's ranking across multiple lists; lower values suggest higher significance.",
      PFP = "Estimated Proportion of False Positives, analogous to the False Discovery Rate (FDR); lower values indicate more reliable findings.",
      pvalue = "Raw p-value from the meta-analysis, indicating the probability of observing the gene's ranking by chance; lower values suggest higher significance.",
      p.adjust = "Adjusted p-value accounting for multiple hypothesis testing using the Benjamini-Hochberg method; helps control the FDR.",
      GenePositions = "Rank positions of the gene in each individual input list; provides insight into the gene's performance across different datasets.",
      Score = "Aggregate score computed to reflect the gene's overall ranking consistency; lower scores indicate higher agreement across lists.")
    
    dt <- datatable(
      selected_data,
      rownames  = FALSE,
      selection = "multiple",
      filter    = "top",
      escape    = FALSE,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        scrollY = table_height,
        autoWidth = FALSE,
        columnDefs = list(list(
          targets = "_all",
          className = "dt-wrap",
          width = paste0(round(100 / max(1, ncol(selected_data))), "%")
        )),
        headerCallback = htmlwidgets::JS(
          "function(thead, data, start, end, display){",
          "  var tooltips = ", jsonlite::toJSON(colTooltips4), ";",
          "  $('th', thead).each(function(){",
          "    var txt = $(this).text().trim();",
          "    if(tooltips[txt]){",
          "      var label = $('<span>').addClass('dt-header-label').text(txt);",
          "      var icon  = $('<span>?</span>')",
          "                     .addClass('dt-info-icon')",
          "                     .attr('title', tooltips[txt]);",
          "      $(this).empty().append($('<div>').addClass('dt-header-with-icon')",
          "                               .append(label).append(icon));",
          "    }",
          "  });",
          "}")))
    
    if (input$selectedPackage == "Rankprodpack") {
      dt <- dt %>% formatSignif(columns = intersect(c("RP_stat", "PFP", "pvalue", "p.adjust"), colnames(selected_data)), digits = 3)
    } else if (input$selectedPackage == "RRApack") {
      dt <- dt %>% formatSignif(columns = intersect(c("Score", "p.adjust"), colnames(selected_data)), digits = 3)}
    
    return(dt)})
  
  
  output$downloadButtonUI <- renderUI({
    req(combined_results())
    fluidRow(
      column(5, downloadButton("downloadResultsCSV", "Download CSV", style = "width: 100%;")),
      column(5, downloadButton("downloadResultsTSV", "Download TSV", style = "width: 100%;")),
      column(1, actionButton("toggleExcludedGenes", 
                label = NULL, 
                icon = icon("eye"), 
                style = "width: 100%;", 
                title = "Show/hide table of excluded genes with unique appearance")),
      column(1, downloadButton("downloadExcludedGenesList", 
                label = NULL, 
                icon = icon("download"), 
                style = "width: 100%;", 
                title = "Download excluded gene list (.txt)")))})
  
  # DOWNLOAD CSV
  output$downloadResultsCSV <- downloadHandler(
    filename = function() {
      paste0("MetaRank_Results_", Sys.Date(), ".csv")},
    content = function(file) {
      req(combined_results(), input$selectedPackage, input$Analysis_results_table_rows_all)
      selected_columns <- switch(input$selectedPackage,
                                 "Rankprodpack" = input$selectedColumnsResultsRP,
                                 "RRApack" = input$selectedColumnsResultsRRA, NULL)
      req(selected_columns)
      data <- combined_results()
      filtered_data <- data[input$Analysis_results_table_rows_all, selected_columns, drop = FALSE]
      write.csv(filtered_data, file, row.names = FALSE)})
  
  
  # DOWNLOAD TSV
  output$downloadResultsTSV <- downloadHandler(
    filename = function() {
      paste0("MetaRank_Results_", Sys.Date(), ".tsv")},
    content = function(file) {
      req(combined_results(), input$selectedPackage, input$Analysis_results_table_rows_all)
      selected_columns <- switch(input$selectedPackage,
                                 "Rankprodpack" = input$selectedColumnsResultsRP,
                                 "RRApack" = input$selectedColumnsResultsRRA, NULL)
      req(selected_columns)
      data <- combined_results()
      filtered_data <- data[input$Analysis_results_table_rows_all, selected_columns, drop = FALSE]
      write.table(filtered_data, file, sep = "\t", row.names = FALSE, quote = FALSE)})
  
  
  # DOWNLOAD EXCLUDED
  output$downloadExcludedGenesList <- downloadHandler(
    filename = function() paste0("Excluded_Genes_List_", Sys.Date(), ".txt"),
    content = function(file) {
      excluded <- prepared_data()$excluded
      file_names <- prepared_data()$files
      combined <- do.call(rbind, lapply(seq_along(excluded), function(i) {
        df <- excluded[[i]]
        if (nrow(df) > 0) df$Source <- file_names[i]
        df}))
      if (!"Gene" %in% colnames(combined)) return()
      
      combined <- combined |>
        group_by(Gene) |>
        summarise(Source = paste(unique(Source), collapse = ", "), .groups = "drop")
      
      writeLines(sort(unique(combined$Gene)), con = file)})

  
  # EXCLUDED TABLE
  output$excludedGenesTable <- DT::renderDT({
    prep <- prepared_data()
    req(prep$excluded, prep$files)
    
    combined <- do.call(rbind, lapply(seq_along(prep$excluded), function(i) {
      df <- prep$excluded[[i]]
      if (nrow(df) > 0) df$Source <- prep$files[i]
      df}))
    req(nrow(combined) > 0, "No excluded genes.")
    
    combined <- combined |>
      group_by(Gene) |>
      summarise(Source = paste(unique(Source), collapse = ", "), .groups = "drop")
    
    raw   <- gene_data()
    names <- prep$files
    all_genes <- if (input$selectedPackage=="RRApack") raw else lapply(raw, `[[`, "Gene")
    app_df <- count_gene_appearance_rra(all_genes, names)
    
    combined <- combined |>
      left_join(app_df[, c("Gene","Count")], by = "Gene")
    
    datatable(
      setNames(combined[, c("Gene", "Count", "Source")], c("GeneID", "FileCount", "FileNames")),
      options   = list(pageLength = 5, autoWidth = TRUE),
      rownames  = FALSE,
      selection = "multiple"
    )})
  
  
  
  
  # UPSET PLOT (PLOT + 2 DOWNLOAD BUTTONS)
  #----------------------------------------------------------------------------
  upsetPlot <- reactive({
    req(gene_data())
    create_upset_plot(
      gene_data = gene_data(),
      text_size = input$upsetTextSize,
      sets_color = input$upsetSetsColor,
      intersection_color = input$upsetIntersectionColor)})
  
  output$upsetPlot <- renderPlot({
    req(upsetPlot())
    upsetPlot()})
  
  output$downloadButtonUpset <- renderUI({
    req(upsetPlot())
    fluidRow(
      column(6, downloadButton("downloadUpsetPNG", "Download PNG", style = "width: 100%;")),
      column(6, downloadButton("downloadUpsetJPG", "Download JPG", style = "width: 100%;")))})
  
  # DOWNLOAD PNG
  output$downloadUpsetPNG <- downloadHandler(
    filename = function() {
      paste0("UpSetR_Plot_", Sys.Date(), ".png")},
    content = function(file) {
      save_upset_plot_transparent(upsetPlot(), file)},
    contentType = "image/png")
  
  # DOWNLOAD JPG
  output$downloadUpsetJPG <- downloadHandler(
    filename = function() {
      paste0("UpSetR_Plot_", Sys.Date(), ".jpg")},
    content = function(file) {
      jpeg(file, width = 14.5 * 96, height = 8.33 * 96, res = 96)  # Abrir dispositivo JPEG
      print(create_upset_plot(
        gene_data = gene_data(),
        text_size = input$upsetTextSize,
        sets_color = input$upsetSetsColor,
        intersection_color = input$upsetIntersectionColor))
      dev.off()},
    contentType = "image/jpg")
  
  
  
  # HEATMAP (PLOT + MULTIPLE DOWNLOAD BUTTON (PNG, JPG, HTML))
  #----------------------------------------------------------------------------
  heatmap_obj <- reactive({
    req(gene_data())
    generate_heatmap(
      gene_data = gene_data(),
      text_size = input$heatmapTickSize,
      title_size = input$heatmapTitleSize,
      tick_size = input$heatmapTickSize,
      colorscale = input$heatmapColorscale)})
  
  heatmap_flat_obj <- reactive({
    req(gene_data())
    generate_flat_heatmap(
      gene_data = gene_data(),
      text_size = input$heatmapTickSize,
      title_size = input$heatmapTitleSize,
      tick_size = input$heatmapTickSize,
      colorscale = input$heatmapColorscale)})
  
  output$heatmapPlot <- renderPlotly({
    req(heatmap_obj())
    ggplotly(heatmap_obj())})
  
  output$downloadHeatmapUI <- renderUI({
    req(combined_results())
    dropdownButton(
      label = "Download Heatmap",
      icon = icon("download"),
      circle = FALSE,
      status = "secondary",
      width = "300px",
      tags$h4("Download Format"),
      downloadButton("downloadHeatmapHTML", "HTML", class = "btn-block"),
      downloadButton("downloadHeatmapPNG", "PNG", class = "btn-block"),
      downloadButton("downloadHeatmapJPG", "JPG", class = "btn-block"))})

  # DOWNLOAD HTML 
  output$downloadHeatmapHTML <- downloadHandler(
  filename = function() {
    paste("MetaRank_Heatmap_", Sys.Date(), ".html", sep = "")},
  content = function(file) {
    htmlwidgets::saveWidget(ggplotly(heatmap_obj()), file, selfcontained = TRUE)})

  # DOWNLOAD PNG
  output$downloadHeatmapPNG <- downloadHandler(
    filename = function() {
      paste("MetaRank_Heatmap_", Sys.Date(), ".png", sep = "")},
    content = function(file) {
      ggsave(file, plot = heatmap_flat_obj(), device = "png", width = 10, height = 8, dpi = 300)})
  
  # DOWNLOAD JPG
  output$downloadHeatmapJPG <- downloadHandler(
    filename = function() {
      paste("MetaRank_Heatmap_", Sys.Date(), ".jpg", sep = "")},
    content = function(file) {
      ggsave(file, plot = heatmap_flat_obj(), device = "jpeg", width = 10, height = 8, dpi = 300)})
  
  
  
  #----------------------------------------------------------------------------
  # DATA REPRESENTATION OF ENRICHMENT ANALYSIS (CARD)
  #----------------------------------------------------------------------------
  
  # ENRICHMENT ANALYSIS (TABLE + 2 BUTTONS + PLOT + 3 BUTTONS)
  #----------------------------------------------------------------------------
  output$enrichmentTable <- DT::renderDataTable({
    req(enrichment_results())
    req(input$selectedColumnsRank)
    
    data <- as.data.frame(enrichment_results(), stringsAsFactors = FALSE)
    if ("ID" %in% colnames(data)) data$ID <- as.character(data$ID)
    
    selected_data <- data[, input$selectedColumnsRank, drop = FALSE]
    order_idx <- which(colnames(selected_data) == "p.adjust")
    order_list <- if (length(order_idx) > 0) {
      list(list(order_idx - 1, "asc"))
    } else {
      list()}
    
    num_rows    <- input$enrichmentTable_length
    table_h     <- paste0(num_rows * 40, "px")
    
    colTooltips3 <- list(
      ID = paste(
        "Stable identifier of the enriched term or pathway.",
        "Examples: GO accessions (e.g. GO:0008150),",
        "KEGG codes (e.g. 00010), Reactome IDs (e.g. R-HSA-199420)."),
      Description = paste(
        "Descriptive name conveying the biological role of the term or pathway.",
        "E.g. “cellular response to chemical stimulus” (GO),",
        "“Glycolysis / Gluconeogenesis” (KEGG), etc."),
      GeneRatio = paste(
      "Proportion of your input genes associated with this term.",
      "Calculated as k/n, where:",
        "- k = number of genes in your input that are associated with this term",
        "- n = total number of genes in your input list",
      "This indicates how enriched the term is among your selected genes."),
      
      BgRatio = paste(
      "Proportion of background genes associated with this term.",
      "Calculated as K/N, where:",
        "- K = number of background genes associated with this term",
        "- N = total number of genes in the background",
      "This reflects how common the term is in the overall gene universe used for the analysis."),
      
      pvalue = paste(
        "Raw p-value from the hypergeometric test,",
        "the probability of ≥ k genes in the term by chance,",
        "where k = ListCount, n = input size, M = GeneCount, N = background size."),
      p.adjust = paste(
        "Adjusted p-value controlling the False Discovery Rate",
        "via Benjamini–Hochberg correction."),
      Count = paste(
        "Total number of background genes annotated to this term (M in the test).",
        "Defines the expected gene set size under the null."),
      geneID = paste(
        "List of gene symbols from the input that map to this term,",
        "showing exactly which genes drive the enrichment."))
    
    DT::datatable(
      selected_data,
      rownames  = FALSE,
      selection = "multiple",
      filter    = "top",
      escape    = FALSE,
      options   = list(
        pageLength = 10,
        scrollX    = TRUE,
        scrollY    = table_h,
        autoWidth  = FALSE,
        columnDefs = list(list(
          targets   = "_all",
          className = "dt-wrap",
          width     = paste0(round(100 / max(1, ncol(selected_data))), "%")
        )),
        order         = order_list,
        headerCallback = htmlwidgets::JS(
          "function(thead, data, start, end, display){",
          "  var tooltips = ", jsonlite::toJSON(colTooltips3), ";",
          "  $('th', thead).each(function(){",
          "    var txt = $(this).text().trim();",
          "    if(tooltips[txt]){",
          "      var label = $('<span>').addClass('dt-header-label').text(txt);",
          "      var icon  = $('<span>?</span>')",
          "                     .addClass('dt-info-icon')",
          "                     .attr('title', tooltips[txt]);",
          "      $(this).empty().append($('<div>').addClass('dt-header-with-icon')",
          "                               .append(label).append(icon));",
          "    }",
          "  });",
          "}"
        ))) %>%
      formatSignif(
        columns = intersect(c("pvalue", "p.adjust"), colnames(selected_data)),
        digits  = 4)})
  
  output$downloadButtonUI2 <- renderUI({
    req(combined_results())
    fluidRow(
      column(6, downloadButton("downloadEnrichmentCSV", "Download CSV", style = "width: 100%;")),
      column(6, downloadButton("downloadEnrichmentTSV", "Download TSV", style = "width: 100%;")))})

  # DOWNLOAD CSV
  output$downloadEnrichmentCSV <- downloadHandler(
    filename = function() {
      paste("MetaRank_results_Enrich", Sys.Date(), ".csv", sep = "")},
    content = function(file) {
      req(enrichment_results(), input$selectedColumnsRank, input$enrichmentTable_rows_all)
      data <- as.data.frame(enrichment_results())
      filtered_data <- data[input$enrichmentTable_rows_all, input$selectedColumnsRank, drop = FALSE]
      write.csv(filtered_data, file, row.names = FALSE)})
  
  # DOWNLOAD TSV
  output$downloadEnrichmentTSV <- downloadHandler(
    filename = function() {
      paste("MetaRank_results_Enrich", Sys.Date(), ".tsv", sep = "")},
    content = function(file) {
      req(enrichment_results(), input$selectedColumnsRank, input$enrichmentTable_rows_all)
      data <- as.data.frame(enrichment_results())
      filtered_data <- data[input$enrichmentTable_rows_all, input$selectedColumnsRank, drop = FALSE]
      write.table(filtered_data, file, row.names = FALSE, sep = "\t")})
  
  
  # REACTIVE PLOT FOR ENRICHMENT ANALYSIS
  filtered_enrichment <- reactive({
    req(enrichment_results())
    if (!is.null(input$enrichmentTable_rows_all)) {
      enrichment_results()[as.numeric(input$enrichmentTable_rows_all), ]
    } else {
      enrichment_results()}})
  
  plot_obj <- reactive({
    req(filtered_enrichment(), input$plotLabel)
    enriched_data <- filtered_enrichment()
    num_terms <- min(input$nPlotTerms, nrow(enriched_data))
    text_size <- input$textSize
    
    data_plot <- enriched_data[1:num_terms, ]
    
    if (input$plotLabel == "ID_Description") {
      data_plot$ID_Description <- paste(data_plot$ID, data_plot$Description, sep = " - ")
      y_aes <- "ID_Description"
    } else {
      y_aes <- input$plotLabel}
    
    if (input$plotType == "dotplot") {
      p <- ggplot(data_plot, aes(x = Count, y = reorder(.data[[y_aes]], -p.adjust))) +
        geom_point(aes(size = Count, color = p.adjust)) +
        scale_color_gradient(low = input$colorLow, high = input$colorHigh) +
        labs(title = "ORA Dot Plot", x = "Count", y = NULL, color = "p.adjust") +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = text_size),
          axis.text.x = element_text(size = text_size),
          axis.title = element_text(size = text_size + 1),
          plot.title = element_text(size = text_size + 2, face = "bold", hjust = 0.5),
          legend.text = element_text(size = text_size - 2),
          legend.title = element_text(size = text_size - 1))
    } else {
      p <- ggplot(data_plot, aes(x = Count, y = reorder(.data[[y_aes]], -p.adjust))) +
        geom_bar(stat = "identity", aes(fill = p.adjust)) +
        scale_fill_gradient(low = input$colorLow, high = input$colorHigh) +
        labs(title = "ORA Bar Plot", x = "Count", y = NULL, fill = "p.adjust") +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = text_size),
          axis.text.x = element_text(angle = 45, hjust = 1, size = text_size),
          axis.title = element_text(size = text_size + 1),
          plot.title = element_text(size = text_size + 2, face = "bold", hjust = 0.5),
          legend.text = element_text(size = text_size - 2),
          legend.title = element_text(size = text_size - 1))}
    
    return(p)})
  
  output$enrichmentDotPlot <- renderPlotly({
    req(plot_obj())
    ggplotly(plot_obj())})
  
  output$downloadPlotUI3 <- renderUI({
    req(enrichment_results())
    dropdownButton(
      label = "Download Plot",
      icon = icon("download"),
      circle = FALSE,
      status = "secondary",
      width = "300px",
      tags$h4("Download Format"),
      downloadButton("downloadPlotHTML", "HTML", class = "btn-block"),
      downloadButton("downloadPlotPNG", "PNG", class = "btn-block"),
      downloadButton("downloadPlotJPG", "JPG", class = "btn-block"))})
  
  # DOWNLOAD HTML
  output$downloadPlotHTML <- downloadHandler(
    filename = function() {
      paste("MetaRank_Enrich_plot_", Sys.Date(), ".html", sep = "")},
    content = function(file) {
      htmlwidgets::saveWidget(ggplotly(plot_obj()), file, selfcontained = TRUE)})
  
  # DOWNLOAD PNG 1200x800px
  output$downloadPlotPNG <- downloadHandler(
    filename = function() {
      paste("MetaRank_Enrich_plot_", Sys.Date(), ".png", sep = "")},
    content = function(file) {
      ggsave(file, plot = plot_obj(), device = "png",
             width = 14.5, height = 8.33, units = "in", dpi = 96)})
  
  # DOWNLOAD JPG 1200x800px
  output$downloadPlotJPG <- downloadHandler(
    filename = function() {
      paste("MetaRank_Enrich_plot_", Sys.Date(), ".jpg", sep = "")},
    content = function(file) {
      ggsave(file, plot = plot_obj(), device = "jpeg",
             width = 14.5, height = 8.33, units = "in", dpi = 96)})
  
  


  
  #----------------------------------------------------------------------------
  # NAVIGATION EXTRAS (INFORMATION BUTTONS + DATA UPDATE)
  #----------------------------------------------------------------------------
  
  # HIDE BUTTONS BEFORE THE ANALYSIS
  output$Hidebuttons <- reactive({
    !is.null(combined_results())})
  outputOptions(output, "Hidebuttons", suspendWhenHidden = FALSE)
  
  
  # EXAMPLE OF FILE FORMAT INFORMATION (I BUTTON NEXT TO UPLOAD A FILE)
  observeEvent(input$showExample, {
    showModal(
      modalDialog(
        title = tags$div("Example File Contents for RP (Weighted)", style = "text-align: center; font-weight: bold; font-size: 20px;"),
        easy_close = TRUE,
        size = "l",
        fluidRow(
          column(6, 
                 tags$h5("Tab-Separated Format (.tsv)"),
                 HTML("<pre style='background-color:#f8f9fa; padding:10px; border-radius:5px;'>COL10A1\t12.01\nCOL11A1\t11.74\nGREM1\t10.74\nMMP1\t10.11\nMMP12\t9.29\nSPINK1\t8.97\n...</pre>")),
          column(6, 
                 tags$h5("Comma-Separated Format (.csv)"),
                 HTML("<pre style='background-color:#f8f9fa; padding:10px; border-radius:5px;'>COL10A1,12.01\nCOL11A1,11.74\nGREM1,10.74\nMMP1,10.11\nMMP12,9.29\nSPINK1,8.97\n...</pre>"))),
        br(),
        tags$h5("Important Notes", style = "margin-bottom: 20px;"),
        tags$div(style = "margin-top: 10px;",
                 tags$p(HTML("<b>• Accepted formats:</b> <code>.tsv</code>, <code>.csv</code>, and <code>.txt</code>. Files with the <code>.txt</code> extension are automatically interpreted as tab-separated."),
                        style = "margin-bottom: 6px;"),
                 tags$p(HTML("<b>• One gene per line:</b> Each row should contain a gene name, optionally followed by a numeric value (e.g., expression score or fold change)."),
                        style = "margin-bottom: 6px;"),
                 tags$p(HTML("<b>• Decimal format:</b> Please use a <code>dot (.)</code> for decimal values. Commas <code>(,)</code> are not allowed as decimal separators."),
                        style = "margin-bottom: 6px;"),
                 tags$p(HTML("<b>• Pre-sorting is optional:</b> Gene lists do not need to be ordered by value, but doing so can improve meta-analysis accuracy."),
                        style = "margin-bottom: 6px;"),
                 tags$p(HTML("<b>• Unlimited lists:</b> You can upload or paste as many gene lists as you want — there is no hard limit."),
                        style = "margin-bottom: 6px;"),
                 tags$p(HTML("<b>⚠ No headers:</b> Do <u>not</u> include any header row in your files. Input should start directly with gene names."),
                        style = "margin-bottom: 0px;")),
        footer = modalButton("Close")))})
  
  
  # EXAMPLE OF FILE FORMAT INFORMATION RRA (I BUTTON NEXT TO UPLOAD A FILE)
  observeEvent(input$showExampleRRA, {
    showModal(
      modalDialog(
        title = tags$div("Example File Contents for RRA (Unweighted)", style = "text-align: center; font-weight: bold; font-size: 20px;"),
        easy_close = TRUE,
        size = "l",
        fluidRow(
          column(12, 
                 tags$h5("Text Format (.txt)"),
                 HTML("<pre style='background-color:#f8f9fa; padding:10px; border-radius:5px;'>BRCA1\nEGFR\nMYC\nPTEN\nTP53\n...</pre>"))),
        br(),
        tags$h5("Important Notes", style = "margin-bottom: 20px;"),
        tags$div(style = "margin-top: 10px;",
                 tags$p(HTML("<b>• Accepted format:</b> Only <code>.txt</code> files are supported."),
                        style = "margin-bottom: 6px;"),
                 tags$p(HTML("<b>• One gene per line:</b> Each row must contain exactly one gene name."),
                        style = "margin-bottom: 6px;"),
                 tags$p(HTML("<b>• Ranked input required:</b> Gene lists must be <u>pre-ranked</u>. The algorithm uses the gene order for aggregation."),
                        style = "margin-bottom: 6px;"),
                 tags$p(HTML("<b>• Unlimited lists:</b> You can upload as many gene lists as needed — there is no upper limit."),
                        style = "margin-bottom: 6px;"),
                 tags$p(HTML("<b>⚠ No headers:</b> Do <u>not</u> include any header row. Input should start directly with gene names."),
                        style = "margin-bottom: 0px;")),
        footer = modalButton("Close")))})
  
  
  # EXAMPLE OF ORIGIN MANAGMENT (I BUTTON NEXT TO ORIGIN) 
  observeEvent(input$showOriginInfo, {
    showModal(
      modalDialog(
        title = tags$div("Origin Information", style = "text-align: center; font-weight: bold; font-size: 20px;"),
        size = "l",
        div(style = "text-align: center;",
            img(src = "org.png", height = "500px", width = "500px")),
        br(),
        tags$h5("Important Notes", style = "margin-bottom: 20px;"),
        tags$div(style = "margin-top: 10px;",
                 tags$p(HTML("<b>• Purpose:</b> This input allows grouping your gene lists by metadata such as <i>laboratory</i>, <i>technology</i>, <i>acquisition date</i>, etc."),
                        style = "margin-bottom: 6px;"),
                 tags$p(HTML("<b>• Format:</b> Enter a numeric vector separated by commas <code>(,)</code>, such as <code>1,2,1,2,2</code>."),
                        style = "margin-bottom: 6px;"),
                 tags$p(HTML("<b>• Required length:</b> The number of values <u>must match</u> the number of gene lists you’ve uploaded or pasted."),
                        style = "margin-bottom: 6px;"),
                 tags$p(HTML("<b>• Group replication:</b> Each group must contain at least two lists. Single-list groups are not supported."),
                        style = "margin-bottom: 6px;"),
                 tags$p(HTML("<b>•⚠ Important:</b> Valid numeric vector for the example data is <code>1,2,1,2.</code>"),
                        style = "margin-bottom: -6px;")),
        easy_close = TRUE,
        footer = modalButton("Close")))})
  
  
  # SHOW/HIDE EXCLUDED GENES TABLE 
  observeEvent(input$toggleExcludedGenes, {
    toggle("excludedGenesSection")})
  
  
  # SLIDER AUTO-UPDATE 
  observeEvent({
    input$files
    input$useExampleData
    input$pasteText
    input$pasteText1
    input$files_rra
    input$useExampleData_rra
    input$pasteText_rra}, 
    
    {total_lists <- NULL
    
    if (input$selectedPackage == "Rankprodpack") {
      if (input$useExampleData && !is.null(example_files)) {
        total_lists <- length(example_files)
      } else if (input$inputMethod == "files" && !is.null(input$files)) {
        total_lists <- length(input$files$datapath)
      } else if (input$inputMethod == "paste") {
        if (input$pasteFormat == "csv" && !is.null(input$pasteText) && nzchar(input$pasteText)) {
          total_lists <- length(strsplit(input$pasteText, "###")[[1]])
        } else if (input$pasteFormat == "tsv" && !is.null(input$pasteText1) && nzchar(input$pasteText1)) {
          total_lists <- length(strsplit(input$pasteText1, "###")[[1]])}}
      
      if (!is.null(total_lists) && total_lists > 0) {
        updateSliderInput(session, "gene_app_filter_rp",
                          min = 1,
                          max = total_lists,
                          step = 1,
                          value = total_lists)}
      
    } else if (input$selectedPackage == "RRApack") {
      if (input$useExampleData_rra && !is.null(example_files_rra)) {
        total_lists <- length(example_files_rra)
      } else if (input$inputMethod == "files" && !is.null(input$files_rra)) {
        total_lists <- length(input$files_rra$datapath)
      } else if (input$inputMethod == "paste" && !is.null(input$pasteText_rra) && nzchar(input$pasteText_rra)) {
        total_lists <- length(strsplit(input$pasteText_rra, "###")[[1]])}
      
      if (!is.null(total_lists) && total_lists > 0) {
        updateSliderInput(session, "gene_app_filter_rra",
                          min = 1,
                          max = total_lists,
                          step = 1,
                          value = total_lists)}}})
  
  selected_gene_app_filter <- reactive({
    if (input$selectedPackage == "Rankprodpack") {
      input$gene_app_filter_rp
    } else {
      input$gene_app_filter_rra}})
  
  
  # RANKPRODPACK EXAMPLE ZIP
  output$downloadExampleRank <- downloadHandler(
    filename = function() { "MetaRank_example_Rankprodpack.zip" },
    content = function(zipfile) {
      # copiar a temp
      tmp <- tempdir()
      files_cp <- file.path(tmp, basename(example_files))
      file.copy(example_files, files_cp, overwrite = TRUE)
      # zip
      owd <- setwd(tmp); on.exit(setwd(owd), add=TRUE)
      utils::zip(zipfile, basename(files_cp), extras = "-j")},
    contentType = "application/zip")
  
  # ROBUSTRANKAGGREGPACK EXAMPLE ZIP
  output$downloadExampleRRA <- downloadHandler(
    filename = function() { "MetaRank_example_RRApack.zip" },
    content = function(zipfile) {
      tmp <- tempdir()
      files_cp <- file.path(tmp, basename(example_files_rra))
      file.copy(example_files_rra, files_cp, overwrite = TRUE)
      owd <- setwd(tmp); on.exit(setwd(owd), add=TRUE)
      utils::zip(zipfile, basename(files_cp), extras = "-j")},
    contentType = "application/zip")
  
  
  prev_input_method <- reactiveVal("files")
  
  observeEvent(input$inputMethod, {
    prev <- prev_input_method()
    
    if (input$inputMethod == "paste" && prev == "files") {
      shinyjs::reset("files")
      shinyjs::reset("files_rra")
      
      updateTextInput(session, "pasteText", value = "")
      updateTextInput(session, "pasteText1", value = "")
      updateTextInput(session, "pasteText_rra", value = "")
      
      # Reset switches
      updateSwitchInput(session, "useExampleData", value = FALSE)
      updateSwitchInput(session, "useExampleData_rra", value = FALSE)
      
    } else if (input$inputMethod == "files" && prev == "paste") {
      updateTextInput(session, "pasteText", value = "")
      updateTextInput(session, "pasteText1", value = "")
      updateTextInput(session, "pasteText_rra", value = "")}
    
    prev_input_method(input$inputMethod)})
  
  
  
 
}

shinyApp(ui = ui, server = server)
# MetaRank_Functions.R
library(dplyr)
library(ggplot2)
library(RankProd)
library(bslib)
library(DT)
library(shinycssloaders)
library(shinyWidgets)
library(shiny)
library(shinyBS)
library(plotly)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(msigdbr)
library(UpSetR)
library(ReactomePA)
library(reactome.db)
library(RobustRankAggreg)
library(viridis)
library(RColorBrewer)
library(magick)

###############################################################################

#----------------------- NAVIGATION FUNCTIONS --------------------------------#

###############################################################################


#' Check if enrichment results are empty and show modal if so
#' @param results  Data.frame or object with enrichment results
#' @param session  Shiny session object
#' @return TRUE if results are non-empty; FALSE (and shows modal) if empty

check_empty_results <- function(results, session) {
  is_empty <- is.null(results) || (is.data.frame(results) && nrow(results) == 0)
  
  if (is_empty) {
    showModal(modalDialog(
      title = tags$div("⚠️ No Enrichment Results Found", 
                       style = "text-align: center; font-weight: bold; font-size: 20px;"),
      easyClose = TRUE,
      size = "l",
      tags$h4("No biological terms matched your input genes"),
      tags$p("The enrichment analysis did not find any relevant pathways, terms or functions associated with the provided gene list."),
      tags$hr(),
      tags$h4("Suggestions to improve your input"),
      tags$ul(
        tags$li("Check the organism selection (e.g., Human, Mouse or Rat)."),
        tags$li("Verify gene identifiers and formats (e.g., SYMBOL, ENTREZID or ENSEMBL)."),
        tags$li("Ensure your genes are annotated in the selected database."),
        tags$li("Use more genes if possible — small lists often yield no results.")),
      tags$p("If the problem persists, try example data to confirm the tool is working."),
      footer = modalButton("Close")))
    return(FALSE)}
  return(TRUE)}


#' Validate pasted RRA format (one column only, lists separated by ###)
#' @param input Shiny input object
#' @param session Shiny session object
#' @return TRUE if valid; FALSE (and shows modal) if any invalid genes found

validate_rra_format <- function(input, session) {
  if (!(input$selectedPackage == "RRApack" &&
        input$inputMethod     == "paste"   &&
        nzchar(input$pasteText_rra))) {
    return(TRUE)}
  
  txt <- input$pasteText_rra
  has_sep   <- grepl("###", txt, fixed = TRUE)
  bad_chars <- grepl("[,;|:\\\\t\\]", txt)
  
  if (!has_sep || bad_chars) {
    showModal(modalDialog(
      title    = tags$div("⚠️ Unweighted Paste Format Error",
                          style = "text-align:center;font-weight:bold;font-size:20px;"),
      easyClose = TRUE, size = "l",
      tags$h5("Correct Format Examples", style = "margin-bottom:10px;"),
      fluidRow(
        column(4, tags$h6("SYMBOL"),
               HTML("<pre style='background-color:#f8f9fa; padding:10px; border-radius:5px;'>TP53\nMYC\n###\nBRCA1\nPTEN\n...</pre>")),
        column(4, tags$h6("ENTREZID"),
               HTML("<pre style='background-color:#f8f9fa; padding:10px; border-radius:5px;'>7157\n4609\n###\n672\n5728\n...</pre>")),
        column(4, tags$h6("ENSEMBL"),
               HTML("<pre style='background-color:#f8f9fa; padding:10px; border-radius:5px;'>ENSG00000141510\nENSG00000136997\n###\nENSG00000012048\nENSG00000171862\n...</pre>"))),
      br(), tags$h5("Possible Issues"),
      tags$div(style = "margin-top:10px;",
               tags$p("- You must include '###' on its own line to separate multiple lists."),
               tags$p("- Do NOT use commas (,), tabs, periods (.), hyphens (-), semicolons (;), or pipes (|) in gene lines."),
               tags$p("- Paste only raw gene identifiers, one per line.")),
      
      br(), tags$h5("Recommendations"),
      tags$div(style = "margin-top:10px;",
               tags$p("- Use '###' to split each ranked list."),
               tags$p("- Avoid CSV/TSV formats or pasted tables."),
               tags$p("- Each line must contain exactly one gene identifier.")),
      
      footer = modalButton("Close")))
    shinyjs::disable("Run")
    shinyjs::delay(100, shinyjs::enable("Run"))
    return(FALSE)}
  return(TRUE)}


#' Validate the origin vector for rp.adv
#' @param input Shiny input object
#' @param session Shiny session object
#' @param pd List returned por prepared_data(), debe incluir pd$valid
#' @return TRUE if valid; FALSE (and shows modal) if any invalid genes found

validate_origin <- function(input, session, pd) {
  if (!(input$selectedPackage == "Rankprodpack" && input$RankBasedMethod == "rp.adv")) {
    return(TRUE)}
  req(pd$valid)
  n_lists <- length(pd$valid)
  origin_raw <- input$origin
  
  if (is.null(origin_raw) || !nzchar(origin_raw)) {
    showModal(modalDialog(
      title = "⚠️ Origin Required",
      "You must provide an origin vector for 'RankProd.Advanced' before running the analysis. Please check the correct usage format in the information button next to origin input.",
      easyClose = TRUE,
      footer = modalButton("Close")))
    return(FALSE)}
  
  if (grepl("[\\.\\s;\\t/:\\-\\|_\\\\]", origin_raw)) {
    showModal(modalDialog(
      title = "⚠️ Invalid Separator in Origin",
      "Please use commas (,) only—no spaces, semicolons or tabs—in your origin vector. Please check the correct usage format in the information button next to origin input.",
      easyClose = TRUE,
      footer = modalButton("Close")))
    return(FALSE)}
  
  if (!grepl("^[0-9,]+$", origin_raw)) {
    showModal(modalDialog(
      title = tags$div("⚠️ Invalid Input Format",
                       style = "text-align: center; font-weight: bold; font-size: 20px;"),
      easyClose = TRUE,
      size = "l",
      tags$h4("A numeric vector separated by commas is required.", style = "margin-bottom: 10px;"),
      tags$p("Please enter a list of numeric values separated by commas, without any letters or special characters."),
      tags$pre("Numeric vector of the example data: 1,2,1,2",
               style = "background-color:#f8f9fa; padding:10px; border-radius:5px;"),
      footer = modalButton("Close")))
    return(FALSE)}
  
  origin_vals <- strsplit(gsub("\\s+", "", origin_raw), ",")[[1]]
  if (length(origin_vals) != n_lists) {
    showModal(modalDialog(
      title = "⚠️ Origin Length Mismatch",
      paste0("Expected ", n_lists, " values (one per list), but got ", length(origin_vals), 
             ". Please check the correct usage format in the information button next to origin input."),
      easyClose = TRUE,
      footer = modalButton("Close")))
    return(FALSE)}
  
  counts <- table(origin_vals)
  if (any(counts == 1)) {
    showModal(modalDialog(
      title = "⚠️ Insufficient Replicates in Origin",
      "Each origin label must appear at least twice. Check your origin vector. Please check the correct usage format in the information button next to origin input.",
      easyClose = TRUE,
      footer = modalButton("Close")))
    return(FALSE)}
  return(TRUE)}


#' Validate gene nomenclature against selected organism
#' @param gene_list Character vector of genes to validate
#' @param organism  One of "Hsa", "Mmu", "Rno"
#' @param session   Shiny session (para showModal)
#' @return TRUE if valid; FALSE (and shows modal) if any invalid genes found

validate_organism <- function(gene_list, organism, IDtype, session) {
  if (length(gene_list) == 0) return(TRUE)
  
  # Define pattern based on ID_type and organism
  pattern <- switch(IDtype,
                    SYMBOL = switch(organism,
                                    Hsa = "^[A-Z][A-Z0-9._/\\-]*$",         # Human SYMBOL
                                    Mmu = "^[A-Z][a-z0-9._/\\-]*$",         # Mouse SYMBOL
                                    Rno = "^[A-Z][a-z0-9._/\\-]*$",         # Rat SYMBOL
                                    NULL),
                    ENTREZID = "^[0-9]+$",                                  # Only digits
                    ENSEMBL = switch(organism,
                                     Hsa = "^ENSG[0-9]{11}$",               # Human ENSEMBL
                                     Mmu = "^ENSMUSG[0-9]{11}$",            # Mouse ENSEMBL
                                     Rno = "^ENSRNOG[0-9]{11}$",            # Rat ENSEMBL
                                     NULL),
                    NULL)
  
  if (is.null(pattern)) return(TRUE)
  
  # Organism display names for modal
  organism_names <- list(
    Hsa = "<i>Homo sapiens</i> (Human)",
    Mmu = "<i>Mus musculus</i> (Mouse)",
    Rno = "<i>Rattus norvegicus</i> (Rat)")
  
  invalid <- unique(gene_list[!grepl(pattern, gene_list)])
  if (length(invalid) > 0) {
    show_genes <- head(invalid, 5)
    extra_note <- if (length(invalid) > 5) paste0("… and ", length(invalid) - 5, " more") else ""
    
    showModal(modalDialog(
      title = tags$div("⚠️ Organism–Gene Nomenclature Issue",
                       style = "text-align:center; font-weight:bold; font-size:20px;"),
      easyClose = TRUE, size = "l",
      
      tags$h5(HTML(paste("Selected organism:", organism_names[[organism]])), style = "margin-bottom:10px;"),
      tags$p("The following gene names do not match the expected format for the selected organism:"),
      tags$pre(paste(c(show_genes, extra_note), collapse = "\n"),
               style = "background-color:#f8f9fa; padding:10px; border-radius:5px;"),
      tags$h5("Expected Formats:", style = "margin-top:20px;"),
      fluidRow(
        column(4, tags$h6(HTML("Format for <i>Homo sapiens</i> (Human)"), style = "margin-top:15px;"),
               HTML("<pre style='background-color:#f8f9fa; padding:10px; border-radius:5px;'>BRCA1\n672\nENSG00000012048</pre>"),
               tags$p(HTML("<b>• SYMBOL:</b> All uppercase letters (e.g., <code>BRCA1</code>, <code>TP53</code>).")),
               tags$p(HTML("<b>• ENTREZID:</b> Numerical only, specific to NCBI Entrez (e.g., <code>672</code>).")),
               tags$p(HTML("<b>• ENSEMBL:</b> Ensembl human genes start with <code>ENSG</code> (e.g., <code>ENSG00000012048</code>)."))),
        column(4, tags$h6(HTML("Format for <i>Mus musculus</i> (Mouse)"), style = "margin-top:15px;"),
               HTML("<pre style='background-color:#f8f9fa; padding:10px; border-radius:5px;'>Brca1\n12189\nENSMUSG00000017146</pre>"),
               tags$p(HTML("<b>• SYMBOL:</b> Capitalized first letter, rest lowercase (e.g., <code>Brca1</code>, <code>Tp53</code>).")),
               tags$p(HTML("<b>• ENTREZID:</b> Numerical only (e.g., <code>12189</code>).")),
               tags$p(HTML("<b>• ENSEMBL:</b> Ensembl mouse genes start with <code>ENSMUSG</code> (e.g., <code>ENSMUSG00000017146</code>)."))),
        column(4, tags$h6(HTML("Format for <i>Rattus norvegicus</i> (Rat)"), style = "margin-top:15px;"),
               HTML("<pre style='background-color:#f8f9fa; padding:10px; border-radius:5px;'>Brca1\n24169\nENSRNOG00000020708</pre>"),
               tags$p(HTML("<b>• SYMBOL:</b> Capitalized first letter, rest lowercase (e.g., <code>Brca1</code>, <code>Tp53</code>).")),
               tags$p(HTML("<b>• ENTREZID:</b> Numerical only (e.g., <code>24169</code>).")),
               tags$p(HTML("<b>• ENSEMBL:</b> Ensembl rat genes start with <code>ENSRNOG</code> (e.g., <code>ENSRNOG00000020708</code>).")))),
      
      footer = modalButton("Close")))
    return(FALSE)}
  return(TRUE)}


#' Validate geneid format
#' @param gene_list Character vector of genes to validate
#' @param ID_type  One of "SYMBOL", "ENSEMBL", "ENTREZID"
#' @param session   Shiny session (para showModal)
#' @return TRUE if valid; FALSE (and shows modal) if any invalid genes found

validate_gene_ids <- function(gene_list, ID_type) {
  gene_list <- as.character(gene_list)
  
  ensembl_pattern <- "^ENS[A-Z]{0,4}[0-9]{6,}$"
  pattern <- switch(ID_type,
                    SYMBOL   = "^[A-Za-z][A-Za-z0-9._/\\-]*$",
                    ENTREZID = "^[0-9]+$",
                    ENSEMBL  = ensembl_pattern,
                    stop("Unsupported ID type"))
  
  is_valid <- grepl(pattern, gene_list)
  
  if (ID_type == "SYMBOL") {
    is_ensembl_like <- grepl(ensembl_pattern, gene_list)
    is_valid <- is_valid & !is_ensembl_like}
  
  invalid <- gene_list[!is_valid]
  
  if (length(invalid) > 0) {
    show_genes <- head(invalid, 10)
    extra_note <- if (length(invalid) > 10) paste0("… and ", length(invalid) - 10, " more") else ""
    
    showModal(modalDialog(
      title = tags$div("⚠️ Gene ID Format Error",
                       style = "text-align: center; font-weight: bold; font-size: 20px;"),
      easyClose = TRUE,
      size = "l",
      
      tags$h5(paste("Selected ID type:", ID_type), style = "margin-bottom: 10px;"),
      tags$p("The following gene identifiers do not match the expected format:"),
      tags$pre(paste(c(show_genes, extra_note), collapse = "\n"),
               style = "background-color:#f8f9fa; padding:10px; border-radius:5px;"),
      
      tags$h5("Expected Formats:", style = "margin-top:20px;"),
      fluidRow(
        column(4, tags$h6(HTML("Format for <i>Homo sapiens</i> (Human)"), style = "margin-top:10px;"),
               HTML("<pre style='background-color:#f8f9fa; padding:10px; border-radius:5px;'>BRCA1\n672\nENSG00000012048</pre>"),
               tags$p(HTML("<b>• SYMBOL:</b> All uppercase letters (e.g., <code>BRCA1</code>, <code>TP53</code>).")),
               tags$p(HTML("<b>• ENTREZID:</b> Numerical only, e.g., <code>672</code>.")),
               tags$p(HTML("<b>• ENSEMBL:</b> Starts with <code>ENSG</code> (e.g., <code>ENSG00000012048</code>)."))),
        column(4, tags$h6(HTML("Format for <i>Mus musculus</i> (Mouse)"), style = "margin-top:10px;"),
               HTML("<pre style='background-color:#f8f9fa; padding:10px; border-radius:5px;'>Brca1\n12189\nENSMUSG00000017146</pre>"),
               tags$p(HTML("<b>• SYMBOL:</b> First letter uppercase, rest lowercase (e.g., <code>Brca1</code>).")),
               tags$p(HTML("<b>• ENTREZID:</b> Numerical only, e.g., <code>12189</code>.")),
               tags$p(HTML("<b>• ENSEMBL:</b> Starts with <code>ENSMUSG</code> (e.g., <code>ENSMUSG00000017146</code>)."))),
        column(4, tags$h6(HTML("Format for <i>Rattus norvegicus</i> (Rat)"), style = "margin-top:10px;"),
               HTML("<pre style='background-color:#f8f9fa; padding:10px; border-radius:5px;'>Brca1\n24169\nENSRNOG00000020708</pre>"),
               tags$p(HTML("<b>• SYMBOL:</b> First letter uppercase, rest lowercase (e.g., <code>Brca1</code>).")),
               tags$p(HTML("<b>• ENTREZID:</b> Numerical only, e.g., <code>24169</code>.")),
               tags$p(HTML("<b>• ENSEMBL:</b> Starts with <code>ENSRNOG</code> (e.g., <code>ENSRNOG00000020708</code>).")))),
      
      footer = modalButton("Close")))
    return(FALSE)}
  return(TRUE)}


#' Validate the uploaded files
#' @param input Shiny input object
#' @param session Shiny session object
#' @return TRUE if valid; FALSE (and shows modal) if any file had the wrong format

validate_file_format <- function(input, session) {
  bad_files <- list()
  
  if (input$selectedPackage == "Rankprodpack") {
    files <- if (input$inputMethod == "files") input$files$datapath else NULL
    file_names <- if (input$inputMethod == "files") input$files$name else NULL
  } else if (input$selectedPackage == "RRApack") {
    files <- if (input$inputMethod == "files") input$files_rra$datapath else NULL
    file_names <- if (input$inputMethod == "files") input$files_rra$name else NULL
  } else {
    return(TRUE)}
  
  if (is.null(files)) return(TRUE)
  
  for (i in seq_along(files)) {
    lines <- readLines(files[i], warn = FALSE, encoding = "UTF-8")
    lines <- trimws(lines[nzchar(lines)])

    if (any(grepl("#", lines))) {
      bad_files[[file_names[i]]] <- "Contains invalid characters like '#'."
      next}
    
    if (input$selectedPackage == "RRApack") {
      if (any(grepl("[,\t]", lines))) {
        bad_files[[file_names[i]]] <- "RRA files must contain only one column of gene symbols (no commas or tabs)."}
    } else if (input$selectedPackage == "Rankprodpack") {
      sep_counts <- sapply(strsplit(lines, "[,\t]"), length)
      
      if (!any(sep_counts == 2)) {
        bad_files[[file_names[i]]] <- "Weighted files must contain two columns separated by comma or tab."
      } else if (any(sep_counts > 2)) {
        bad_files[[file_names[i]]] <- paste(
          "Weighted files must not have more than two columns (more than one ‘tab’ or ‘comma’ is detected per line).",
          "This may be due to:",
          " Extra separators (e.g., multiple tabs or commas) in a single line.",
          " Misuse of comma as a decimal separator (use dot '.' instead, e.g., 1.23 not 1,23).",
          sep = "\n")}}}
  
  if (length(bad_files) > 0) {
    showModal(modalDialog(
      title = "⚠️ Format Errors Detected",
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close"),
      tags$ul(
        lapply(names(bad_files), function(f) {
          tags$li(HTML(paste("<b>", f, "</b>: ", bad_files[[f]])))})),
      tags$hr(),
      fluidRow(
        column(12,
               tags$h4("Expected Format"),
               tags$p("Each uploaded file must contain a plain list of gene identifiers, one per line, with no headers or separators."),
               tags$h5("Valid Examples by ID Type and Organism"),
               fluidRow(
                 column(4, wellPanel(
                   tags$strong("SYMBOL / HUMAN"),
                   HTML("<pre>BRCA1\nTP53\nMYC</pre>"))),
                 column(4, wellPanel(
                   tags$strong("SYMBOL / MOUSE"),
                   HTML("<pre>Trp53\nMyc\nEgfr</pre>"))),
                 column(4, wellPanel(
                   tags$strong("SYMBOL / RAT"),
                   HTML("<pre>Tp53\nMyc\nEgfr</pre>")))),
               fluidRow(
                 column(4, wellPanel(
                   tags$strong("ENTREZID / HUMAN"),
                   HTML("<pre>672\n7157\n4609</pre>"))),
                 column(4, wellPanel(
                   tags$strong("ENTREZID / MOUSE"),
                   HTML("<pre>22059\n17869\n13649</pre>"))),
                 column(4, wellPanel(
                   tags$strong("ENTREZID / RAT"),
                   HTML("<pre>24842\n116538\n25499</pre>")))),
               fluidRow(
                 column(4, wellPanel(
                   tags$strong("ENSEMBL / HUMAN"),
                   HTML("<pre>ENSG00000012048\nENSG00000141510\nENSG00000136997</pre>"))),
                 column(4, wellPanel(
                   tags$strong("ENSEMBL / MOUSE"),
                   HTML("<pre>ENSMUSG00000059552\nENSMUSG00000059552\nENSMUSG00000020122</pre>"))),
                 column(4, wellPanel(
                   tags$strong("ENSEMBL / RAT"),
                   HTML("<pre>ENSRNOG00000016647\nENSRNOG00000019902\nENSRNOG00000000007</pre>"))))))))
    return(FALSE)}
  return(TRUE)}


###############################################################################

#----------------------- AUXILIARY FUNCTIONS ---------------------------------#

###############################################################################


# Clean and normalize a vector of gene names
clean_gene_list <- function(lines) {
  lines <- as.character(lines)          # Forzar a character
  lines <- lines[lines != "" & !is.na(lines)]
  lines <- trimws(lines)
  lines <- gsub("///.*|//.*", "", lines)
  unique(lines[lines != ""])}

# Read RRApack gene lists from text files
read_rra_files <- function(paths) {
  lists <- lapply(paths, function(p) {
    lines <- readLines(p, warn = FALSE)
    clean_gene_list(as.character(lines))})
  names(lists) <- basename(paths)
  lists}

# Parse pasted RRApack lists separated by '###'
paste_rra_lists <- function(raw_text) {
  blocks <- unlist(strsplit(raw_text, "###", fixed = TRUE))
  blocks <- blocks[nzchar(trimws(blocks))]
  
  lists <- lapply(blocks, function(b) {
    lines <- unlist(strsplit(b, "\n"))
    clean_gene_list(as.character(lines))})
  
  lists <- Filter(function(x) length(x) > 0, lists)
  
  if (length(lists) == 0) {
    stop("No valid lists found after parsing. Check your ### separators.")}
  
  names(lists) <- paste0("List_", seq_along(lists))
  lists}


#' Read a single RankProdpack file (CSV/TSV/TXT)
#' @param path File path
#' @return Data frame with 'Gene' and 'Stat.data'

read_rankprod_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  df <- switch(ext,
               csv = read.csv(path, header = FALSE, stringsAsFactors = FALSE),
               tsv = read.delim(path, header = FALSE, stringsAsFactors = FALSE, sep = "\t"),
               txt = read.delim(path, header = FALSE, stringsAsFactors = FALSE, sep = "\t"),
               stop("Unsupported file format: ", path))
  if (ncol(df) < 2) {
    warning("Unexpected format: ", path)
    return(NULL)}
  colnames(df)[1:2] <- c("Gene", "Stat.data")
  df <- df %>%
    filter(!is.na(Gene), Gene != "") %>%
    mutate(Gene = sapply(strsplit(Gene, "//"), `[`, 1), Gene = (Gene)) %>%
    distinct(Gene, .keep_all = TRUE)
  buscar_redondeos_v2(df)}


#' Read multiple RankProdpack files
#' @param paths Vector of file paths
#' @return Named list of data frames

read_rankprod_files <- function(paths) {
  dfs <- lapply(paths, read_rankprod_file)
  dfs <- Filter(Negate(is.null), dfs)
  names(dfs) <- basename(paths)
  dfs}


#' Parse pasted RankProdpack blocks (CSV or TSV) separated by '###'
#' @param raw_text Single string
#' @param format 'csv' or 'tsv'
#' @return Named list of data frames

paste_rankprod_lists <- function(raw_text, format = c("csv", "tsv")) {
  fmt <- match.arg(format)
  blocks <- strsplit(raw_text, "###")[[1]]
  
  dfs <- lapply(blocks, function(b) {
    con <- textConnection(b)
    df <- if (fmt == "csv") read.csv(con, header = TRUE, stringsAsFactors = FALSE)
    else read.delim(con, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    close(con)
    
    df <- df %>%
      filter(!is.na(Gene), Gene != "") %>%
      mutate(
        Gene = as.character(Gene),
        Gene = ifelse(grepl("//", Gene), sapply(strsplit(Gene, "//"), `[`, 1), Gene)
      ) %>%
      distinct(Gene, .keep_all = TRUE)
    
    buscar_redondeos_v2(df)})
  
  names(dfs) <- paste0("List_", seq_along(dfs))
  return(dfs)}



#' Flag trailing “.00” numeric values as NA
#'
#' Scans each numeric column in the input data frame and replaces any values
#' that are whole numbers formatted with two or more trailing zeros (e.g. “5.00”,
#' “3.000”) with NA. Intended to catch and remove potential rounding artifacts.
#'
#' @param tabla A data.frame that may include a 'Gene' column and one or more numeric columns.
#' @return The same data.frame, with flagged numeric entries set to NA.

buscar_redondeos_v2 <- function(tabla) {
  if ("Gene" %in% names(tabla)) {
    tabla$Gene <- (tabla$Gene)}
  for (col in names(tabla)) {
    if (is.numeric(tabla[[col]])) {
      valores_formateados <- format(tabla[[col]], scientific = FALSE)
      indices <- which(
        tabla[[col]] == floor(tabla[[col]]) &
          grepl("\\.0{2,}", valores_formateados))
      if (length(indices) > 0) {
        tabla[[col]][indices] <- NA}}}
  return(tabla)}


#' Summarize gene occurrences across multiple ranked lists
#'
#' Builds a summary data frame listing each gene, the files in which it appears,
#' and the total count of distinct files per gene.
#'
#' @param gene_data A list of data.frames, each containing at least a 'Gene' column.
#' @param file_names A character vector of the same length as gene_data, with names for each list.
#' @return A data.frame with columns:
#'   - Gene: gene identifier  
#'   - File: comma-separated list of source files  
#'   - Count: number of distinct files in which the gene appears

count_gene_appearance_with_files <- function(gene_data, file_names) {
  if (length(gene_data) != length(file_names)) {
    stop("Length of 'gene_data' and 'file_names' must match.")}
  
  gene_file_tracker <- lapply(seq_along(gene_data), function(i) {
    data.frame(
      Gene = gene_data[[i]]$Gene,
      File = file_names[i],
      stringsAsFactors = FALSE)})
  
  combined_data <- do.call(rbind, gene_file_tracker)
  appearance_summary <- aggregate(
    File ~ Gene,
    data = combined_data,
    FUN = function(x) paste(unique(x), collapse = ", "))
  
  appearance_summary$Count <- sapply(
    strsplit(appearance_summary$File, ", "),
    length)
  
  return(appearance_summary)}


#' Count gene appearances for RRA-style ranked lists
#'
#' Given a list of ranked gene vectors (no scores, just gene symbols) and corresponding file names,
#' this function summarizes in which files each gene appeared and how many distinct files contained it.
#'
#' @param gene_data A list of character vectors, each representing one RRA input list of gene symbols.
#' @param file_names A character vector of the same length as \code{gene_data}, naming each list.
#' @return A data.frame with columns:
#'   \describe{
#'     \item{Gene}{The gene symbol.}
#'     \item{File}{Comma-separated list of input names where the gene was found.}
#'     \item{Count}{Number of distinct lists in which the gene appears.}

count_gene_appearance_rra <- function(gene_data, file_names) {
  if (length(gene_data) != length(file_names)) {
    stop("Length of 'gene_data' and 'file_names' must be the same.")}

  valid_gene_data <- gene_data[sapply(gene_data, length) > 0]
  valid_file_names <- file_names[sapply(gene_data, length) > 0]
  
  if (length(valid_gene_data) == 0) {
    return(data.frame(Gene = character(0), File = character(0), Count = integer(0), stringsAsFactors = FALSE))}
  
  gene_file_tracker <- do.call(rbind, lapply(seq_along(valid_gene_data), function(i) {
    data.frame(
      Gene = valid_gene_data[[i]],
      File = valid_file_names[i],
      stringsAsFactors = FALSE)}))
  
  appearance_summary <- aggregate(
    File ~ Gene,
    data = gene_file_tracker,
    FUN = function(x) paste(unique(x), collapse = ", "))
  
  appearance_summary$Count <- sapply(
    strsplit(appearance_summary$File, ", "),
    length)
  
  return(appearance_summary)}



#' Generate a statistics matrix for meta-analysis
#'
#' Constructs a matrix (data.frame) of statistic values for each gene across multiple datasets,
#' aligning genes by name and handling missing values according to the chosen strategy.
#'
#' @param gene_data A list of data.frames, each containing at least:
#'   \describe{
#'     \item{Gene}{Character vector of gene symbols.}
#'     \item{Stat.data}{Numeric vector of statistic values (e.g., p-values, fold changes).}
#'   }
#' @param na_management Character string indicating how to treat \code{NA} values:
#'   \itemize{
#'     \item \code{"ignore"}: Leave \code{NA}s as-is (downstream methods may impute or ignore).
#'     \item \code{"penalize"}: Replace \code{NA}s with the worst observed value (max or min) across all datasets.
#'   }
#' @param ranking_direction Character string, either \code{"ascending"} (lower values are better)
#'   or \code{"descending"} (higher values are better). Used when computing the penalty value.
#'
#' @return A data.frame where:
#'   \describe{
#'     \item{rows}{Unique genes (rownames).}
#'     \item{columns}{One column per input dataset, containing the statistic values.}

generate_statdata_matrix <- function(gene_data, na_management = "ignore", ranking_direction = "ascending") {
  
  all_genes <- unique(unlist(lapply(gene_data, function(dataset) dataset$Gene)))
  statdata_matrix <- sapply(gene_data, function(dataset) {
    stat_values <- setNames(as.numeric(dataset$Stat.data), dataset$Gene)
    stat_values[all_genes]})
  
  statdata_matrix <- as.data.frame(statdata_matrix)
  rownames(statdata_matrix) <- all_genes
  
  if (na_management == "penalize") {
    penalty_value <- if (ranking_direction == "ascending") {
      max(statdata_matrix, na.rm = TRUE)
    } else if (ranking_direction == "descending") {
      min(statdata_matrix, na.rm = TRUE)
    } else {
      stop("Unrecognized ranking_direction; use 'ascending' or 'descending'.")}
    statdata_matrix[is.na(statdata_matrix)] <- penalty_value}
  return(statdata_matrix)}


#' Perform basic Rank Products meta-analysis
#'
#' Wraps \code{RankProducts} to compute gene rankings across multiple datasets,
#' optionally tracking progress for a Shiny app.
#'
#' @param logFC_matrix A numeric matrix or data.frame where rows are genes and columns are datasets (values typically log-fold changes or scores).
#' @param cl Numeric vector assigning class labels (e.g., all 1s for one class).
#' @param logged Logical; if \code{TRUE}, data are already log-transformed.
#' @param na.rm Logical; if \code{TRUE}, remove \code{NA} pairs before computing.
#' @param plot Logical; if \code{TRUE}, generate diagnostic plots.
#' @param ranking_direction Character, either \code{"ascending"} (lower values = better) or \code{"descending"}.
#' @param progress Optional function to report progress (e.g., Shiny’s \code{incProgress}).
#' @param MinNumOfValidPairs Numeric; minimum number of valid sample pairs to include a gene.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{Gene}{Gene symbol.}
#'     \item{Rank}{Rank position (1 = best).}
#'     \item{RP_stat}{Rank Products statistic.}
#'     \item{PFP}{Percentage of false positives.}
#'     \item{pval}{P-value for the rank product.}

perform_RP <- function(logFC_matrix, cl, logged = TRUE, na.rm = TRUE, plot = FALSE,
                       ranking_direction = "ascending", progress = NULL, MinNumOfValidPairs = NA) {
  if (!is.null(progress)) progress(0.1, detail = "Starting RP computation...")
  
  result <- RankProducts(
    logFC_matrix,
    cl = cl,
    na.rm = na.rm,
    logged = logged,
    gene.names = rownames(logFC_matrix),
    plot = plot,
    MinNumOfValidPairs = MinNumOfValidPairs)
  
  if (!is.null(progress)) progress(0.3, detail = "RP computation finished. Processing results...")
  rank_col <- ifelse(ranking_direction == "ascending", 1, 2)
  
  ranked_genes <- data.frame(
    Gene    = rownames(logFC_matrix),
    Rank    = result$RPrank[, rank_col],
    RP_stat = result$RPs[, rank_col],
    PFP     = result$pfp[, rank_col],
    pval    = result$pval[, rank_col])
  
  ranked_genes$p.adjust <- p.adjust(ranked_genes$pval, method = "BH")
  
  ranked_genes <- ranked_genes[order(ranked_genes$Rank), ]
  if (!is.null(progress)) progress(0.8, detail = "Results processed.")
  return(ranked_genes)}


#' Perform advanced Rank Products meta-analysis with origin grouping
#'
#' Wraps \code{RP.advance} to compute gene rankings, allowing specification of group origins
#' (e.g., lab, platform) and tracking progress for a Shiny app.
#'
#' @param logFC_matrix A numeric matrix or data.frame where rows are genes and columns are datasets.
#' @param cl Numeric vector assigning class labels.
#' @param origin Numeric vector indicating the group/origin of each dataset column.
#' @param logged Logical; if \code{TRUE}, data are already log-transformed.
#' @param na.rm Logical; if \code{TRUE}, remove \code{NA} pairs before computing.
#' @param plot Logical; if \code{TRUE}, generate diagnostic plots.
#' @param ranking_direction Character, either \code{"ascending"} or \code{"descending"}.
#' @param progress Optional function to report progress.
#'
#' @return A data.frame with the same structure as \code{perform_RP}.

perform_RP_advance <- function(logFC_matrix, cl, origin, logged = TRUE, na.rm = TRUE, plot = FALSE,
                               ranking_direction = "ascending", progress = NULL) {
  if (!is.null(progress)) progress(0.1, detail = "Starting advanced RP computation...")
  
  result <- RP.advance(
    logFC_matrix,
    cl = cl,
    origin = origin,
    na.rm = na.rm,
    logged = logged,
    gene.names = rownames(logFC_matrix),
    plot = plot)
  
  if (!is.null(progress)) progress(0.6, detail = "Advanced RP computation finished. Processing results...")
  rank_col <- ifelse(ranking_direction == "ascending", 1, 2)
  
  ranked_genes <- data.frame(
    Gene    = rownames(logFC_matrix),
    Rank    = result$RPrank[, rank_col],
    RP_stat = result$RPs[, rank_col],
    PFP     = result$pfp[, rank_col],
    pval    = result$pval[, rank_col])
  
  ranked_genes$p.adjust <- p.adjust(ranked_genes$pval, method = "BH")
  
  ranked_genes <- ranked_genes[order(ranked_genes$Rank), ]
  if (!is.null(progress)) progress(1, detail = "Results processed.")
  return(ranked_genes)}


#' Run Robust Rank Aggregation (RRA) meta-analysis
#'
#' Performs rank aggregation across multiple gene lists using the RRA algorithm,
#' then merges with appearance counts for each gene.
#'
#' @param gene_data A named list of character vectors, each a ranked gene list.
#' @param method Character; aggregation method to pass to \code{aggregateRanks} (e.g., "RRA", "min", "mean").
#' @param full Logical; if \code{TRUE}, return the full ranking table, otherwise only top hits.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{Gene}{Gene symbol.}
#'     \item{Ranking}{Rank position assigned by RRA (1 = best).}
#'     \item{Score}{Aggregation score from RRA (lower = more significant).}
#'     \item{Appereance}{Number of input lists in which the gene appears.}
#'     \item{Files}{Comma-separated names of input lists containing the gene.}

run_rra_analysis <- function(gene_data, method = "RRA", full = TRUE) {
  if (!is.list(gene_data) || length(gene_data) == 0) {
    stop("gene_data must be a non-empty list of gene vectors.")}
  
  rra_result <- RobustRankAggreg::aggregateRanks(
    glist   = gene_data,
    N       = NA,
    method  = method,
    full    = full,
    exact   = FALSE,
    topCutoff = NA)
  
  rra_result <- as.data.frame(rra_result)
  rra_result$Ranking <- seq_len(nrow(rra_result))
  colnames(rra_result)[1:2] <- c("Gene", "Score")
  
  # Ajuste BH para Score
  rra_result$p.adjust <- p.adjust(rra_result$Score, method = "BH")
  
  appearance_info <- count_gene_appearance_rra(gene_data, names(gene_data))
  
  final_table <- merge(rra_result, appearance_info, by = "Gene", all.x = TRUE)
  colnames(final_table)[colnames(final_table) == "Count"] <- "Appereance"
  colnames(final_table)[colnames(final_table) == "File"]  <- "Files"
  
  final_table <- final_table[order(final_table$Ranking), 
                             c("Gene", "Ranking", "Score", "p.adjust", "Appereance", "Files")]
  
  colnames(final_table) <- c("GeneID", "Rank", "Score", "p.adjust", "FileCount", "FileNames")
  
  return(final_table)}


#' Create a combined ranking table (basic)
#'
#' Merges the output of the basic RankProd analysis with gene appearance counts.
#'
#' @param ranked_genes Data frame with columns Gene, Rank, RP_stat, PFP, pval.
#' @param gene_appearance Data frame with columns Gene, File, Count indicating in which files each gene appears.
#' @return A data frame combining rankings and appearance metadata, with columns Gene, Rank, Count, File, RP_stat, PFP, pval.

create_combined_rank_table <- function(ranked_genes, gene_appearance) {
  if (!"Gene" %in% colnames(ranked_genes) ||
      !"Gene" %in% colnames(gene_appearance)) {
    stop("Both tables must contain a 'Gene' column.")}
  
  gene_appearance$File <- sapply(strsplit(gene_appearance$File, " "), function(files) {
    paste(basename(files), collapse = " ")})
  merged_data <- merge(ranked_genes, gene_appearance, by = "Gene", all.x = TRUE, sort = FALSE)
  
  merged_data$Rank <- ranked_genes$Rank
  merged_data <- merged_data[, c("Gene","Rank","Count","File","RP_stat","PFP","pval","p.adjust")]
  
  colnames(merged_data) <- c("GeneID", "Rank", "FileCount", "FileNames", 
                             "RP_stat", "PFP", "pvalue", "p.adjust")
  
  return(merged_data)}


#' Create a combined ranking table (with penalization)
#'
#' Applies an extra penalty to genes that appear in fewer lists, adjusting their rank.
#'
#' @param ranked_genes Data frame with columns Gene, Rank, RP_stat, PFP, pval.
#' @param gene_appearance Data frame with columns Gene, File, Count indicating in which files each gene appears.
#' @return A data frame with penalized ranks, columns Gene, Rank, Count, File, RP_stat, PFP, pval.

create_combined_rank_table2 <- function(ranked_genes, gene_appearance) {
  if (!"Gene" %in% colnames(ranked_genes) ||
      !"Gene" %in% colnames(gene_appearance)) {
    stop("Both tables must contain a 'Gene' column.")}
  
  total_lists <- max(gene_appearance$Count, na.rm = TRUE)
  max_rank    <- max(ranked_genes$Rank, na.rm = TRUE)
  
  merged_data <- merge(ranked_genes, gene_appearance, by = "Gene", all.x = TRUE, sort = FALSE)
  merged_data$AdjustedRank <- merged_data$Rank +
    ((total_lists - merged_data$Count) * (max_rank / total_lists))
  merged_data <- merged_data[order(merged_data$AdjustedRank), ]
  merged_data$Rank <- seq_len(nrow(merged_data))
  
  # Selección columnas incluyendo p.adjust
  merged_data <- merged_data[, c("Gene","Rank","Count","File","RP_stat","PFP","pval","p.adjust")]
  
  colnames(merged_data) <- c("GeneID", "Rank", "FileCount", "FileNames", 
                             "RP_stat", "PFP", "pvalue", "p.adjust")
  
  return(merged_data)}


#' Create an UpSet plot from gene lists
#'
#' @param gene_data A named list where each element is either:
#'   - A data.frame with the first column containing gene identifiers, or
#'   - A character vector of gene identifiers.
#' @param text_size Numeric; base font size for plot text (default: 16).
#' @param sets_color Character; hex color for the horizontal set bars (default: "#56B4E9").
#' @param intersection_color Character; hex color for the vertical intersection bars (default: "#666666").
#'
#' @return Invisibly returns the UpSet plot object (invisibly printed).

create_upset_plot <- function(gene_data,
                              text_size = 16,
                              sets_color = "#56B4E9",
                              intersection_color = "#666666") {
  par(bg = "transparent")
  
  gene_lists <- lapply(gene_data, function(entry) {
    if (is.data.frame(entry)) {
      col <- entry[[1]]  
      if (is.null(col)) return(NULL)
      return(unique(na.omit(as.character(col))))
    } else if (is.vector(entry)) {
      return(unique(na.omit(as.character(entry))))
    } else {
      return(NULL)}})
  
  gene_lists <- Filter(Negate(is.null), gene_lists)
  
  if (is.null(names(gene_lists)) || any(names(gene_lists) == "")) {
    names(gene_lists) <- paste0("List_", seq_along(gene_lists))}
  
  if (length(gene_lists) == 0) stop("No hay listas válidas de genes para el UpSet plot.")
  
  all_genes <- unique(unlist(gene_lists))
  binary_data <- lapply(gene_lists, function(x) as.numeric(all_genes %in% x))
  binary_df <- as.data.frame(do.call(cbind, binary_data))
  
  colnames(binary_df) <- names(gene_lists)
  rownames(binary_df) <- all_genes
  
  binary_df[] <- lapply(binary_df, function(col) {
    suppressWarnings(as.numeric(as.character(col)))})
  
  if (!all(sapply(binary_df, is.numeric))) {
    print(str(binary_df))
    stop("Error: La tabla binaria contiene columnas no numéricas.")}
  
  plot <- UpSetR::upset(
    binary_df,
    sets = names(gene_lists),
    sets.bar.color = sets_color,
    main.bar.color = intersection_color,
    order.by = "freq",
    empty.intersections = "on",
    text.scale = c(
      text_size / 16, text_size / 16,
      text_size / 14, text_size / 14,
      text_size / 12, text_size / 12))
  
  return(plot)}


#' Save an UpSet plot as a transparent PNG
#'
#' @param plot The UpSet plot object returned by \code{create_upset_plot()}.
#' @param output_file File path (with .png extension) to save the transparent image.
#' @param width Numeric; width in inches (default: 14.5).
#' @param height Numeric; height in inches (default: 8.33).
#' @param dpi Numeric; resolution in dots per inch (default: 96).
#'
#' @return Invisibly returns \code{NULL} and writes the PNG file to disk.

save_upset_plot_transparent <- function(plot, output_file, width = 14.5, height = 8.33, dpi = 96) {
  temp_file <- tempfile(fileext = ".png")
  png(temp_file, width = width * dpi, height = height * dpi, res = dpi, bg = "white")
  print(plot)
  dev.off()
  
  img <- image_read(temp_file)
  img_trans <- image_transparent(img, "white", fuzz = 10)
  image_write(img_trans, path = output_file, format = "png")
  message("Plot saved in: ", output_file)}


#' Generate an interactive heatmap of proportional intersections between gene lists
#'
#' @param gene_data A named list of gene sets. Each element can be either:
#'   - A data.frame containing a column named "Gene", or  
#'   - A character vector of gene identifiers.
#' @param text_size Numeric; base font size for axis tick labels (default: 12).
#' @param title_size Numeric; font size for plot title (default: 17).
#' @param axis_title_size Numeric; font size for axis titles (default: 18).
#' @param tick_size Numeric; font size for axis tick labels and colorbar title (default: 12).
#' @param colorscale Character or list; color scale name for the heatmap (e.g., "Viridis") (default: "Viridis").
#'
#' @return A plotly heatmap object showing the pairwise proportional intersections between each gene list.

generate_heatmap <- function(gene_data, 
                             text_size = 12, 
                             title_size = 17, 
                             axis_title_size = 18, 
                             tick_size = 12, 
                             colorscale = "Viridis") {
  
  gene_lists <- lapply(gene_data, function(entry) {
    if (is.data.frame(entry) && "Gene" %in% colnames(entry)) {
      unique(na.omit(entry$Gene))
    } else if (is.vector(entry)) {
      unique(na.omit(entry))
    } else {
      NULL}})
  
  gene_lists <- Filter(Negate(is.null), gene_lists)
  names(gene_lists) <- names(gene_data)
  
  if (length(gene_lists) < 2) {
    stop("At least two gene lists are required to generate the heatmap.")}
  
  n <- length(gene_lists)
  intersection_matrix <- matrix(0, nrow = n, ncol = n)
  for (i in seq_along(gene_lists)) {
    for (j in seq_along(gene_lists)) {
      inter <- intersect(gene_lists[[i]], gene_lists[[j]])
      uni   <- union(gene_lists[[i]], gene_lists[[j]])
      intersection_matrix[i, j] <- length(inter) / length(uni)}}
  
  rownames(intersection_matrix) <- names(gene_lists)
  colnames(intersection_matrix) <- names(gene_lists)
  
  heatmap_plot <- plot_ly(
    z = intersection_matrix,
    x = colnames(intersection_matrix),
    y = rownames(intersection_matrix),
    type = "heatmap",
    colorscale = colorscale,
    colorbar = list(
      title     = "Proportional<br>Intersection",
      titleside = "right",
      titlefont = list(size = tick_size))) %>%
    
    layout(
      title = list(
        text = "Intersection Proportions Between Gene Lists",
        font = list(size = title_size)),
      xaxis = list(
        title     = "",
        titlefont = list(size = axis_title_size),
        tickfont  = list(size = tick_size)),
      yaxis = list(
        title     = "",
        titlefont = list(size = axis_title_size),
        tickfont  = list(size = tick_size)),
      margin = list(t = 80))
  
  return(heatmap_plot)}


#' Generate a flat (ggplot2) heatmap of proportional intersections between gene lists
#'
#' @param gene_data A named list of gene sets. Each element may be:
#'   - a data.frame containing a column named "Gene", or  
#'   - a character vector of gene identifiers.
#' @param text_size Numeric; base font size for axis tick labels (default: 12).
#' @param title_size Numeric; font size for the plot title (default: 17).
#' @param axis_title_size Numeric; font size for axis titles (default: 18).
#' @param tick_size Numeric; font size for axis tick labels (default: 12).
#' @param colorscale Character; one of "Viridis", "Cividis", or "Portland" to select the color palette (default: "Viridis").
#' @param output_file Character; base name for saving the plot (not used in this function but kept for consistency) (default: "heatmap").
#'
#' @return A ggplot2 object representing the heatmap of pairwise proportional intersections.

generate_flat_heatmap <- function(gene_data, 
                                  text_size = 12, 
                                  title_size = 17, 
                                  axis_title_size = 18, 
                                  tick_size = 12, 
                                  colorscale = "Viridis", 
                                  output_file = "heatmap") {
  
  gene_lists <- lapply(gene_data, function(entry) {
    if (is.data.frame(entry) && "Gene" %in% colnames(entry)) {
      unique(na.omit(entry$Gene))
    } else if (is.vector(entry)) {
      unique(na.omit(entry))
    } else {
      NULL}})
  
  gene_lists <- Filter(Negate(is.null), gene_lists)
  names(gene_lists) <- names(gene_data)
  
  if (length(gene_lists) < 2) {
    stop("At least two gene lists are required to generate the heatmap.")}
  
  n <- length(gene_lists)
  intersection_matrix <- matrix(0, ncol = n, nrow = n)
  for (i in seq_along(gene_lists)) {
    for (j in seq_along(gene_lists)) {
      inter <- intersect(gene_lists[[i]], gene_lists[[j]])
      uni   <- union(gene_lists[[i]], gene_lists[[j]])
      intersection_matrix[i, j] <- length(inter) / length(uni)}}
  
  rownames(intersection_matrix) <- names(gene_lists)
  colnames(intersection_matrix) <- names(gene_lists)

  df_matrix <- as.data.frame(as.table(intersection_matrix))
  colnames(df_matrix) <- c("Gene_List_1", "Gene_List_2", "Intersection_Score")
  
  if (colorscale == "Viridis") {
    palette <- viridis::viridis(256)
  } else if (colorscale == "Cividis") {
    palette <- viridis::cividis(256)
  } else if (colorscale == "Portland") {
    palette <- RColorBrewer::brewer.pal(9, "YlOrRd")
  } else {
    palette <- viridis::viridis(256)}
  
  heatmap_plot <- ggplot2::ggplot(df_matrix, 
                                  ggplot2::aes(x = Gene_List_1, 
                                               y = Gene_List_2, 
                                               fill = Intersection_Score)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colors = palette) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(size = tick_size, angle = 45, hjust = 1),
      axis.text.y  = ggplot2::element_text(size = tick_size),
      axis.title.x = ggplot2::element_text(size = axis_title_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size),
      plot.title   = ggplot2::element_text(size = title_size, hjust = 0.5)) +
    ggplot2::labs(
      title = "Intersection Proportions Between Gene Lists",
      x     = NULL,
      y     = NULL,
      fill  = "Proportion")
  
  return(heatmap_plot)}


#' Extract top N genes from combined results
#'
#' @param combined_results Data frame returned by combined_results() containing a 'Gene' column
#' @param top_n Integer, number of top genes to select (must be > 0)
#' @return Character vector of the first \code{top_n} gene names

get_top_genes <- function(combined_results, top_n) {
  req(combined_results)  
  req(top_n > 0)         
  
  top_genes <- head(combined_results$GeneID, n = top_n)
  return(top_genes)}


#' Generate position vectors for all genes across all input lists
#' 
#' @param gene_data List of data frames (RankProd) or list of character vectors (RRA)
#' @return Data frame with 'Gene' and 'GenePositions' (list-column of integer vectors)

get_gene_positions <- function(gene_data_output) {
  is_rra <- is.character(gene_data_output[[1]]) 
  
  all_genes <- if (is_rra) {
    unique(unlist(gene_data_output))
  } else {
    unique(unlist(lapply(gene_data_output, function(df) df[[1]])))}
  
  gene_positions <- lapply(all_genes, function(gene) {
    sapply(gene_data_output, function(list_or_df) {
      if (is_rra) {
        match(gene, list_or_df)
      } else {
        match(gene, list_or_df[[1]])}})})
  
  result_df <- data.frame(GeneID = all_genes, stringsAsFactors = FALSE)
  result_df$GenePositions <- gene_positions
  
  result_df$GenePositions <- sapply(result_df$GenePositions, function(pos_vector) {
    paste(ifelse(is.na(pos_vector), "NA", pos_vector), collapse = ",")})
  return(result_df)}


###############################################################################

#-------------------------------- NOT USED -----------------------------------#

###############################################################################

Load_gene_file_tsv <- function(files) {
  gene_data <- lapply(files, function(file) {
    data <- read.delim(file, header = FALSE, stringsAsFactors = FALSE, sep = "\t", row.names = NULL)
    
    if (ncol(data) < 2) {
      warning(paste("El archivo", file, "no tiene el formato esperado."))
      return(NULL)}
    
    colnames(data) <- c("Gene", "Stat.data")
    data <- data %>%
      filter(!is.na(Gene) & Gene != "") %>%
      mutate(Gene = sapply(strsplit(Gene, "//"), `[`, 1)) %>%
      mutate(Gene = (Gene)) %>%
      distinct(Gene, .keep_all = TRUE)
    
    data <- buscar_redondeos_v2(data)
    
    return(data)})
  file_names <- basename(files)
  names(gene_data) <- file_names
  return(gene_data = gene_data)}


Load_gene_file_csv <- function(files) {
  gene_data <- lapply(files, function(file) {
    data <- read.csv(file, header = FALSE, stringsAsFactors = FALSE, row.names = NULL)
    
    if (ncol(data) < 2) {
      warning(paste("El archivo", file, "no tiene el formato esperado."))
      return(NULL)}
    
    colnames(data) <- c("Gene", "Stat.data")
    
    data <- data %>%
      filter(!is.na(Gene) & Gene != "") %>%
      mutate(Gene = sapply(strsplit(Gene, "//"), `[`, 1)) %>%
      mutate(Gene = (Gene)) %>%
      distinct(Gene, .keep_all = TRUE)
    
    data <- buscar_redondeos_v2(data)
    
    return(data)})
  
  file_names <- basename(files)
  names(gene_data) <- file_names
  return(list(gene_data = gene_data, file_names = file_names))}


Load_gene_paste_tsv <- function(raw_text) {
  lists <- strsplit(raw_text, split = "###")[[1]]
  gene_data <- lapply(lists, function(list) {
    con <- textConnection(list)
    data <- read.delim(con, header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names = NULL)
    close(con)
    data <- data %>%
      filter(!is.na(Gene) & Gene != "") %>%
      mutate(Gene = sapply(strsplit(Gene, "//"), `[`, 1)) %>%
      mutate(Gene = (Gene)) %>%
      distinct(Gene, .keep_all = TRUE)
    
    data <- buscar_redondeos_v2(data)
    
    return(data)})
  return(gene_data)}


Load_gene_paste_csv <- function(raw_text) {
  lists <- strsplit(raw_text, split = "###")[[1]]
  gene_data <- lapply(lists, function(list) {
    con <- textConnection(list)
    data <- read.csv(con, header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
    close(con)
    data <- data %>%
      filter(!is.na(Gene) & Gene != "") %>%
      mutate(Gene = sapply(strsplit(Gene, "//"), `[`, 1)) %>%
      mutate(Gene = (Gene)) %>%
      distinct(Gene, .keep_all = TRUE)
    
    data <- buscar_redondeos_v2(data)
    
    return(data)})
  return(gene_data)}


normalize_symbols <- function(genes, species) {
  switch(species,
         "Hsa" = toupper(genes),
         "Mmu" = sub("^(.)", "\\U\\1", tolower(genes), perl = TRUE),
         "Rno" = sub("^(.)", "\\U\\1", tolower(genes), perl = TRUE),
         genes)}


perform_RS <- function(logFC_matrix, cl, logged = TRUE, na.rm = TRUE, plot = FALSE) {
  result <- RankProducts(
    logFC_matrix,
    cl = cl,
    na.rm = na.rm,
    logged = logged,
    gene.names = rownames(logFC_matrix),
    plot = plot,
    calculateProduct = FALSE)
  ranked_genes <- data.frame(
    Gene = rownames(logFC_matrix),
    Rank = result$RPrank[, 2],
    RS_stat = result$RSs[, 1],
    PFP = result$pfp[, 1],
    pval = result$pval[, 1])
  ranked_genes <- ranked_genes[order(ranked_genes$Rank), ]
  return(ranked_genes)}


perform_RSadvance <- function(logFC_matrix, cl, origin, logged = TRUE, na.rm = TRUE, plot = FALSE) {
  result <- RSadvance(
    logFC_matrix,
    cl = cl,
    origin = origin,
    na.rm = na.rm,
    logged = logged,
    gene.names = rownames(logFC_matrix),
    plot = plot)
  ranked_genes <- data.frame(
    Gene = rownames(logFC_matrix),
    Rank = result$RPrank[, 2],
    RS_stat = result$RSs[, 1],
    PFP = result$pfp[, 1],
    pval = result$pval[, 1])
  ranked_genes <- ranked_genes[order(ranked_genes$Rank), ]
  return(ranked_genes)}


save_combined_rank_table <- function(merged_data, output_file) {
  tsv_file <- paste0(output_file, ".tsv")
  csv_file <- paste0(output_file, ".csv")
  
  write.table(
    merged_data, 
    file = tsv_file, 
    sep = "\t", 
    row.names = FALSE, 
    quote = FALSE)
  
  write.csv(
    merged_data, 
    file = csv_file, 
    row.names = FALSE)}


save_heatmap <- function(gene_data, output_file) {
  heatmap <- generate_heatmap(gene_data)
  htmlwidgets::saveWidget(heatmap, file = output_file)
  message("Heatmap saved in: ", output_file)}


perform_ora_enrichment <- function(gene_list, 
                                   ont = "BP", 
                                   pvalue_cutoff = 0.05, 
                                   qvalue_cutoff = 0.2, 
                                   OrgDb = org.Hs.eg.db) {
  if (length(gene_list) == 0) {
    stop("La lista de genes está vacía. Proporcione genes válidos para el análisis.")}
  
  enrich_result <- enrichGO(
    gene          = gene_list,
    OrgDb         = OrgDb,  
    keyType       = "SYMBOL",  
    ont           = ont,       
    pvalueCutoff  = pvalue_cutoff,
    qvalueCutoff  = qvalue_cutoff,
    readable      = TRUE)
  
  if (!is.null(enrich_result) && nrow(enrich_result@result) > 0) {
    enrich_result@result$pvalue <- round(enrich_result@result$pvalue, digits = 8)
    enrich_result@result$p.adjust <- round(enrich_result@result$p.adjust, digits = 8)
    enrich_result@result$qvalue <- round(enrich_result@result$qvalue, digits = 8)}
  
  return(enrich_result)}


ora_result_tsv <- function(ora_result, file_path) {
  ora_df <- as.data.frame(ora_result)
  write.table(ora_df, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)}


ora_result_csv <- function(ora_result, file_path) {
  ora_df <- as.data.frame(ora_result)
  write.csv(ora_df, file = file_path, row.names = FALSE)}

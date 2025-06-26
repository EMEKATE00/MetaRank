# Set the output directory for annotation files
output_dir <- "./database_annotations"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load required libraries
library(org.Hs.eg.db)     # Human gene annotation database
library(org.Mm.eg.db)     # Mouse gene annotation database
library(org.Rn.eg.db)     # Rat gene annotation database
library(org.Dm.eg.db)     # (Optional) Drosophila gene annotation database
library(GO.db)            # Gene Ontology database
library(AnnotationDbi)    # Tools for annotation databases
library(ReactomePA)       # Reactome pathway analysis
library(biomaRt)          # Interface to Ensembl BioMart
library(KEGGREST)         # Interface to KEGG REST API
library(reactome.db)      # Reactome database



################################################################################
#                  FUNCTIONS TO GENERATE GO ANNOTATIONS                        #
################################################################################


# Function to select the appropriate OrgDb object based on the organism code
#----------------------------------------------------------------------------
getOrgDb <- function(organism) {
  if (organism == "Hsa") {
    return(org.Hs.eg.db)
  } else if (organism == "Mmu") {
    return(org.Mm.eg.db)
  } else if (organism == "Rno") {
    return(org.Rn.eg.db)
  } else {
    stop("Unsupported organism code")}}



# Function to generate and save GO annotation files for a given organism and ID type
#------------------------------------------------------------------------------------
generate_GO_files <- function(organism, ID) {
  
  orgdb <- getOrgDb(organism)
  keys_list <- keys(orgdb, keytype = ID)                                          # Retrieve the list of keys (gene identifiers) for the given ID type
  go_data <- AnnotationDbi::select(orgdb,
                                   keys = keys_list,
                                   columns = c(ID, "GO", "ONTOLOGY"),             # Extract GO annotations: gene ID, GO term, and ontology (BP/MF/CC)
                                   keytype = ID)
  go_data <- na.exclude(go_data)                                                  # Remove rows without GO annotation
  
  annotations <- list(
    BP = unique(go_data[go_data$ONTOLOGY == "BP", c("GO", ID)]),
    MF = unique(go_data[go_data$ONTOLOGY == "MF", c("GO", ID)]),                  # Split GO annotations by ontology category and remove duplicates
    CC = unique(go_data[go_data$ONTOLOGY == "CC", c("GO", ID)]))
  
  filename_BP <- file.path(output_dir, paste0(organism, "_GO_BP_", ID, ".txt"))
  filename_MF <- file.path(output_dir, paste0(organism, "_GO_MF_", ID, ".txt"))   # Define output filenames for each ontology category
  filename_CC <- file.path(output_dir, paste0(organism, "_GO_CC_", ID, ".txt"))
  
  write.table(annotations$BP, filename_BP, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(annotations$MF, filename_MF, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)     # Write the annotation files
  write.table(annotations$CC, filename_CC, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  cat("Files generated for", organism, "with ID type", ID, "\n")}



# Function to generate and save the GO term descriptions (TERM2NAME file)
#------------------------------------------------------------------------
generate_term2name <- function() {
  go_terms <- Term(GOTERM)                                                        # Get GO term descriptions
  term2name <- data.frame("GO" = names(go_terms), "NAME" = go_terms)
  out_fn <- file.path(output_dir, "TERM2NAME_GO.txt")
  write.table(term2name, out_fn, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  cat("TERM2NAME_GO.txt file generated\n")}



#-----------------------------------------------------------------------------#
#                              GO SCRIPT EXECUTION                            #
#-----------------------------------------------------------------------------#


organisms <- c("Hsa", "Mmu", "Rno")
ID_types <- c("ENTREZID", "SYMBOL", "ENSEMBL")                                    # Define the list of organisms (by code) and the types of gene IDs to process

for (org in organisms) {
  orgdb <- getOrgDb(org)
  available_IDs <- keytypes(orgdb)                                                # Loop through each organism and each ID type to generate GO annotation files
  for (id in ID_types) {
    if (id %in% available_IDs) {
      generate_GO_files(org, id)
    } else {
      cat("ID type", id, "is not available for organism", org, "\n")}}}

generate_term2name()                                                               # Generate the TERM2NAME file for GO terms







################################################################################
#                FUNCTIONS TO GENERATE KEGG ANNOTATIONS                        #
################################################################################

# Function to select the appropriate OrgDb object based on the organism code
#----------------------------------------------------------------------------
getGeneNameDb <- function(organism) {
  if (organism == "Hsa") {
    return(org.Hs.egGENENAME)
  } else if (organism == "Mmu") {
    return(org.Mm.egGENENAME)
  } else if (organism == "Rno") {
    return(org.Rn.egGENENAME)
  } else {
    stop("Unsupported organism code")}}


# Function to generate and save KEGG pathway annotation files for a given organism and ID type
#----------------------------------------------------------------------------------------------
generate_KEGG_files <- function(organism, ID) {
  orgdb <- getOrgDb(organism)
  
  if (!(ID %in% keytypes(orgdb))) {
    cat("ID type", ID, "is not available for", organism, "\n")                    # Check if the requested ID type is available for the organism
    return(NULL)}

  gene_keys <- keys(orgdb, keytype = ID)                                          # Retrieve gene identifiers of the selected type
  kegg_data <- AnnotationDbi::select(orgdb,
                                     keys = gene_keys,
                                     columns = c(ID, "PATH"),                     # Extract KEGG pathway annotations (ID and PATH)
                                     keytype = ID)
  
  kegg_data <- kegg_data[, c("PATH", ID)]                                         # Keep only the columns of interest: pathway and gene ID
  kegg_data <- na.exclude(kegg_data)                                              # Remove entries with missing annotations
  filename <- file.path(output_dir, paste0(organism, "_KEGG_", ID, ".txt"))       # Define the output filename
  write.table(kegg_data, filename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)    # Write the KEGG annotation file
  
  cat("File generated:", filename, "\n")}


# Function to generate and save gene ID to name mapping files (ID2NAME)
#-----------------------------------------------------------------------
generate_geneName_file <- function(organism) {
  geneNameDb <- getGeneNameDb(organism)
  id2name <- toTable(geneNameDb)
  filename <- file.path(output_dir, paste0("ID2NAME_", organism, ".txt"))                          # Define the output filename (one file per organism)
  write.table(id2name, filename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)   # Write the mapping file
  cat("File generated:", filename, "\n")}


# Function to generate a KEGG pathway description file for a specific organism
#----------------------------------------------------------------------
generate_KEGG_pathway_desc <- function(organism) {

  code_kegg <- tolower(organism)
  pw_list <- KEGGREST::keggList("pathway", code_kegg)
  
  df <- data.frame(
    pathway_id   = sub("^path:", "", names(pw_list)),                                # Build a data frame with pathway ID and description
    description  = pw_list,                                                          # Pathway name/description
    stringsAsFactors = FALSE)
  
  df$description <- sub(" –.*$", "", df$description)                                 # Clean up the description to remove text after "–" (dash)
  filename <- file.path(output_dir, paste0(organism, "_KEGG_pathway_desc.txt"))      # Define output file path and write the data frame to a tab-separated file
  write.table(df, filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  cat("File generated:", filename, "\n")}


# Function to generate general KEGG pathway descriptions (map-level, not species-specific)
#----------------------------------------------------------------------
generate_general_Kegg_map <- function() {

  pw_map <- KEGGREST::keggList("pathway", "map")
  df <- data.frame(                                                                    # Build a data frame with cleaned pathway IDs and descriptions
    pathway_id = sub("^map", "",                                                       # Remove "map" prefix
                     sub("^path:", "", names(pw_map))),                                # Remove "path:" prefix
    description = sub(" –.*$", "", pw_map),                                            # Clean description
    stringsAsFactors = FALSE)
  
  write.table(df,
              file.path(output_dir, "TERM2NAME_KEGG.txt"),                             # Write the table to a file as TERM2NAME_KEGG.txt
              sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  
  cat("File generado:", file.path(output_dir, "TERM2NAME_KEGG.txt"), "\n")}



#-----------------------------------------------------------------------------#
#                            KEGG SCRIPT EXECUTION                            #
#-----------------------------------------------------------------------------#

organisms <- c("Hsa", "Mmu", "Rno")
ID_types <- c("ENTREZID", "SYMBOL", "ENSEMBL")

for (org in organisms) {
  for (id in ID_types) {                                                            # Loop through each organism and ID type to generate KEGG annotation files
    generate_KEGG_files(org, id)}

  generate_geneName_file(org)
  generate_KEGG_pathway_desc(org)}                                                  # Generate gene name mapping file

generate_general_Kegg_map()







################################################################################
#                FUNCTIONS TO GENERATE REACTOME ANNOTATIONS                    #
################################################################################

species_mart <- list(
  Hsa = "hsapiens_gene_ensembl",
  Mmu = "mmusculus_gene_ensembl",                                                             # Mapping of organism codes to Ensembl dataset names
  Rno = "rnorvegicus_gene_ensembl")

ID_types <- c("ENTREZID", "SYMBOL", "ENSEMBL")
attr_map <- list(
  ENTREZID = "entrezgene_id",                                                                 # ID types to process and their corresponding biomaRt attributes
  SYMBOL   = "external_gene_name",
  ENSEMBL  = "ensembl_gene_id")


# Function to generate Reactome annotation files for a given organism and ID type
#---------------------------------------------------------------------------------
generate_Reactome_files <- function(org_code, mart_dataset, ID_type) {
  message("Processing ", org_code, " - ", ID_type, " ...")

  mart <- useEnsembl(biomart = "genes", dataset = mart_dataset)                              # Connect to Ensembl BioMart
  df <- getBM(
    attributes = c(attr_map[[ID_type]], "reactome"),                                         # 1) Retrieve gene-to-pathway associations
    mart       = mart)
  df <- df[!is.na(df[,1]) & !is.na(df[,2]), , drop=FALSE]                                    # Remove rows with missing gene ID or Reactome pathway ID
  
  df2 <- data.frame(
    term = df[, 2],
    gene = df[, 1],                                                                          # 2) Reformat columns: 'term' = pathway ID, 'gene' = gene ID
    stringsAsFactors = FALSE)
  
  out_fn <- file.path(output_dir, paste0(org_code, "_Reactome_", ID_type, ".txt"))           # 3) Write to file with appropriate headers
  write.table(df2, out_fn, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  message("  -> File generated: ", out_fn)}



# Function to generate TERM2NAME_REACTOME.txt: mapping from Reactome pathway ID to name
#---------------------------------------------------------------------------------------
generate_term2name_reactome <- function() {
  message("Generating TERM2NAME_REACTOME.txt ...")
  
  if (!requireNamespace("reactome.db", quietly = TRUE))                                          # Ensure the reactome.db package is installed
    BiocManager::install("reactome.db")

  all_ids <- keys(reactome.db, keytype = "PATHID")                                               # Retrieve all Reactome pathway IDs
  t2n <- AnnotationDbi::select(
    reactome.db,
    keys    = all_ids,
    columns = c("PATHID", "PATHNAME"),                                                           # Get pathway names corresponding to each ID
    keytype = "PATHID")
  colnames(t2n) <- c("term", "name")                                                             # Rename columns to match expected format
  
  out_fn <- file.path(output_dir, "TERM2NAME_REACTOME.txt")
  write.table(t2n, out_fn, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  message("  -> File generated: ", out_fn)}


#-----------------------------------------------------------------------------#
#                        REACTOME SCRIPT EXECUTION                            #
#-----------------------------------------------------------------------------#

for (org in names(species_mart)) {
  ds <- species_mart[[org]]
  for (id in ID_types) {                                                                         # Loop through all organisms and ID types to generate Reactome annotation files
    tryCatch(generate_Reactome_files(org, ds, id),
      error = function(e) message("  ! Error with ", org, "/", id, ": ", e$message))}}

generate_term2name_reactome()                                                                    # Generate Reactome term-to-name mapping file


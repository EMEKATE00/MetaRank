library(data.table)
library(clusterProfiler)

# Documentation for ORA enrichment functions with clusterProfiler
#
# These helper functions perform Over-Representation Analysis (ORA) for 
# Gene Ontology (GO), KEGG pathways, and Reactome pathways using custom 
# annotation files and clusterProfiler::enricher.
#
# General Notes:
# – All functions assume a global variable `annotation_dir` points to the folder 
#   containing annotation files.
# – Missing annotation files produce an error and halt execution.
# – Functions return an empty data.frame (rather than NULL) when no significant 
#   terms are identified.

annotation_dir <- file.path(getwd(), "database_annotations")

#------------------------------------------------------------------------------
# Helper to normalize the organism
#
#    • Input:
#        – organism: character string identifying the species 
#                    ("Hsa"/"human", "Mmu"/"mouse", "Rno"/"rat").
#    • Operation:
#        – Maps supported organism names to a standardized 3-letter code.
#    • Output:
#        – Returns one of "Hsa", "Mmu", or "Rno", or stops with an error 
#          if the organism is not supported.
#
#------------------------------------------------------------------------------

getOrgCode <- function(organism) {
  switch(organism,
         Hsa   = "Hsa", human = "Hsa",
         Mmu   = "Mmu", mouse = "Mmu",
         Rno   = "Rno", rat   = "Rno",
         stop("Organismo no soportado: ", organism))}


#------------------------------------------------------------------------------
# ORA GO
#
#    • Inputs:
#        – gene_list: character vector of gene IDs to test.
#        – organism: species name (see getOrgCode).
#        – ontology: GO sub-ontology ("BP", "CC", or "MF").
#        – pvalue: numeric p‑value cutoff for enrichment.
#        – ID_type: gene identifier type (e.g. "ENTREZID", "ENSEMBL").
#    • Operation:
#        1. Translate organism to code via getOrgCode().
#        2. Build path to TERM2GENE file for GO (org_code, ontology, ID_type).
#        3. Load term-to-gene mappings (TERM2GENE).
#        4. Load term-to-name mappings from TERM2NAME_GO.txt.
#        5. Call clusterProfiler::enricher with TERM2GENE, TERM2NAME, and 
#           pvalueCutoff.
#    • Output:
#        – A data.frame of enrichment results (columns: ID, Description, 
#          GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count), 
#          or an empty data.frame if no significant terms are found.
#
#------------------------------------------------------------------------------

ORA_GO <- function(gene_list, organism, ontology, pvalue, ID_type) {
  org_code <- getOrgCode(organism)
  
  file_gene <- file.path(
    annotation_dir,
    paste0(org_code, "_GO_", ontology, "_", ID_type, ".txt"))
  
  if (!file.exists(file_gene)) stop("Falta archivo GO: ", file_gene)
  term2gene <- fread(file_gene, col.names = c("term","gene"))
  
  file_name <- file.path(annotation_dir, "TERM2NAME_GO.txt")
  if (!file.exists(file_name)) stop("Falta TERM2NAME_GO.txt")
  term2name <- fread(file_name, col.names = c("term","name"))
  
  enr <- clusterProfiler::enricher(gene         = gene_list,
                                   TERM2GENE    = term2gene,
                                   TERM2NAME    = term2name,
                                   pvalueCutoff = pvalue)
  
  if (is.null(enr) || nrow(enr@result)==0) return(data.frame())
  as.data.frame(enr@result)}


#------------------------------------------------------------------------------
# ORA KEGG
#
#    • Inputs:
#        – gene_list: character vector of gene IDs.
#        – organism: species name (see getOrgCode).
#        – pvalue: numeric p‑value cutoff.
#        – ID_type: gene identifier type.
#    • Operation:
#        1. Translate organism to code.
#        2. Load KEGG TERM2GENE for the organism and ID_type.
#        3. Load general TERM2NAME_KEGG.txt for pathway names.
#        4. Strip the first three-letter prefix from KEGG IDs in TERM2GENE.
#        5. Run clusterProfiler::enricher with pvalueCutoff.
#        6. Format result IDs as zero-padded 5-digit strings.
#    • Output:
#        – A data.frame of KEGG enrichment results with formatted IDs, or 
#          an empty data.frame if no pathways are enriched.
#
#------------------------------------------------------------------------------

ORA_KEGG <- function(gene_list, organism, pvalue, ID_type) {
  org_code <- getOrgCode(organism)
  
  file_t2g <- file.path(annotation_dir, paste0(org_code, "_KEGG_", 
                                                            ID_type, ".txt"))
  
  if (!file.exists(file_t2g)) stop("Missing KEGG annotation file: ", file_t2g)
  term2gene <- data.table::fread(file_t2g,
                                 col.names = c("term","gene"),
                                 colClasses = c(term="character", 
                                                            gene="character"))

  term2gene[, term := as.character(term)]
  file_t2n <- file.path(annotation_dir, "TERM2NAME_KEGG.txt")
  if (!file.exists(file_t2n)) stop("Missing general TERM2NAME file: ", file_t2n)
  term2name <- data.table::fread(file_t2n,
                                 col.names = c("term","name"),
                                 colClasses = c(term="character", 
                                                            name="character"))
  
  term2name[, term := as.character(term)]
  term2gene[, term := sub("^[a-z]{3}", "", term)]
  
  enr <- clusterProfiler::enricher(gene         = gene_list,
                                   TERM2GENE    = term2gene,
                                   TERM2NAME    = term2name,
                                   pvalueCutoff = pvalue)
  
  if (is.null(enr) || nrow(enr@result)==0) return(data.frame())
  res <- as.data.frame(enr@result, stringsAsFactors = FALSE)
  res$ID <- sprintf("%05d", as.integer(res$ID))
  res}


#------------------------------------------------------------------------------
# ORA Reactome
#
#    • Inputs:
#        – gene_list: character vector of gene IDs.
#        – organism: species name (see getOrgCode).
#        – pvalue: numeric p‑value cutoff.
#        – ID_type: gene identifier type.
#    • Operation:
#        1. Translate organism to code.
#        2. Load Reactome TERM2GENE for the organism and ID_type.
#        3. Load TERM2NAME_REACTOME.txt for pathway names.
#        4. Remove any prefix before “: ” in pathway names.
#        5. Execute clusterProfiler::enricher with pvalueCutoff.
#    • Output:
#        – A data.frame of Reactome enrichment results, or an empty data.frame 
#          if there are no significant pathways.
#
#------------------------------------------------------------------------------

ORA_REACTOME <- function(gene_list, organism, pvalue, ID_type) {
  org_code <- getOrgCode(organism)
  
  file_gene <- file.path(
    annotation_dir,
    paste0(org_code, "_Reactome_", ID_type, ".txt"))
  
  if (!file.exists(file_gene)) stop("Falta archivo Reactome genes: ", file_gene)
  term2gene <- fread(file_gene, col.names = c("term","gene"))
  
  file_name <- file.path(annotation_dir, "TERM2NAME_REACTOME.txt")
  if (!file.exists(file_name)) stop("Falta TERM2NAME_REACTOME.txt")
  term2name <- fread(file_name, col.names = c("term","name"))
  term2name$name <- sub("^[^:]+: ", "", term2name$name)
  
  enr <- clusterProfiler::enricher(gene         = gene_list,
                                   TERM2GENE    = term2gene,
                                   TERM2NAME    = term2name,
                                   pvalueCutoff = pvalue)
  
  if (is.null(enr) || nrow(enr@result)==0) return(data.frame())
  as.data.frame(enr@result)}

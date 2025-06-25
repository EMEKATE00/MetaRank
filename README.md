---
title: "MetaRank Usage Instructions"
format:
  html:
    toc: true
    toc-depth: 4
    number-sections: true
    theme: cosmo
    self-contained: true
---

```{=html}
<style>
body {
  text-align: justify;
  font-size: 0.95em;
}
</style>
```

------------------------------------------------------------------------

MetaRank is a Shiny-based application designed for non-parametric meta-analysis of ranked gene lists. It allows users to combine multiple pre-ranked gene lists into a single consensus ranked list using robust statistical methods (RankProd and RobustRankAggreg). In addition, MetaRank allows functional interpretation by performing overrepresentation analysis (ORA) on the top ranked genes (up to the top 100) of the consensus list.

This comprehensive tutorial provides step-by-step instructions on how to load ranked lists, customise ranking parameters, and visualise or interpret biologically enriched terms from the resulting consensus ranking.

1.  **Overview**
2.  **Meta-analysis Options** (RankProd and RobustRankAggreg)
    -   Rank Product (RankProd)
    -   Robust Rank Aggregation (RRA)
    -   Summary table
3.  **RankProd Workflow**
    -   Inputs
    -   Parameters
    -   Outputs
4.  **RobustRankAggreg Workflow**
    -   Inputs
    -   Parameters
    -   Outputs
5.  **Shared Elements**
    -   Shared Plots
    -   Shared Enrichment Analysis
    -   Data Visualization Tab
6.  **Background Pipeline**
7.  **Best Practices**
8.  **Troubleshooting**

<br>

## Overview

MetaRank is a user-friendly Shiny application designed to unify ranked gene lists from multiple studies and extract a consensus ranking of the most consistently relevant genes. It provides two distinct analytical workflows, `Rank Product (RP)` and `Robust Rank Aggregation (RRA)`, allowing users to choose the method that best fits their data structure and research goals. The app includes a rich set of features for data input, configuration, enrichment, and visualization. Its key features include:

-   **Choice of meta-analysis algorithm**:

    -   `Rank Product (RP):` A **weighted** approach that incorporates both gene rankings and associated scores as pvalues or fold changes.
    -   `Robust Rank Aggregation (RRA):` A **non-weighted** method based on rank positions only, suitable for plain ordered gene lists.

-   **Flexible input options and parameter settings tailored to each method**:

    -   **RP mode**:
        -   Supports both file upload and text paste, using `"###"` to separate lists.
        -   Accepts `.TXT`, `.CSV`, and `.TSV` files containing genes with numerical values.
        -   Includes example datasets for quick testing and file downloads.
        -   Offers a choice between *basic* and *advanced* Rank Product functions.
        -   Customizable settings for handling `NA` values and filtering genes by minimum list appearance.
    -   **RRA mode**:
        -   Accepts plain-text files or pasted input with one gene per line, using `"###"` to separate lists.
        -   Only `.txt` format is supported to avoid structure conflicts.
        -   Provides example data and downloadable templates.
        -   Includes several aggregation options such as "*RRA*", *geometric mean*, *median* or *minimum rank.*
        -   Also allows configuration of `NA` handling and list inclusion thresholds.

-   **Post-ranking functional enrichment**:

    -   Enables *Over-Representation Analysis (ORA)* on the top-ranked genes (selectable from 10 to 100).
    -   Compatible with *Gene Ontology (GO)*, *KEGG*, and *Reactome* databases.
    -   Supports multiple organisms: *Homo sapiens*, *Mus musculus*, and *Rattus norvegicus*.
    -   Accepts gene identifiers in `SYMBOL`, `ENTREZID`, and `ENSEMBL` formats.

-   **Interactive visualization and result export**:

    -   Explore input overlap with **heatmaps** and **UpSet plots**.
    -   View enriched terms using interactive **bar plots** and **dot plots**.
    -   All result tables are interactive and downloadable in `.CSV` or `.TSV` formats.
    -   Customise each graph and table in the `Data Visualisation` tab.

<br>

![](www/images/metarank_overview.png){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: -7px;">
Figure 1: MetaRank overview
</p>

## Meta-analysis Options

MetaRank provides two robust statistical methods for integrating ranked gene lists from multiple studies: `Rank Product (RankProd)` and `Robust Rank Aggregation (RRA)`. These methods are designed to identify genes that consistently appear at the top of ranked lists, thereby highlighting potential candidates for further biological investigation.

### Rank Product (RankProd)

The Rank Product method is a non-parametric statistical approach for identifying differentially expressed genes based on the consistency of gene rankings across multiple datasets. It is particularly suitable for meta-analyses that combine results from different studies, as it does not rely on data normality and is relatively robust to outliers.

Key features:

-   **Non-parametric analysis**: Does not assume any particular data distribution, allowing application to heterogeneous datasets.
-   **Geometric mean aggregation**: Ranks are combined using the geometric mean, giving greater weight to genes consistently ranked at the top.
-   **False discovery rate estimation**: Provides an estimate of the proportion of false positives (pfp) to evaluate the statistical significance of the results.
-   **Cross-platform applicability**: Designed to integrate data from diverse experimental conditions or technologies.

Recommended use cases:

-   Works best with **complete and balanced gene lists**, where most genes are consistently represented across all datasets.
-   Suitable when datasets have similar quality and measurement platforms, minimizing unwanted variability.
-   Performs optimally with a **moderate number of lists** (approximately 5 to 20), ensuring a balance between sensitivity and computational cost.
-   Less effective in the presence of **high proportions of missing values** or when gene representation is inconsistent, although these limitations can be mitigated using filtering based on minimum gene appearance and applying penalization strategies.
-   Execution **time is relatively high**, especially when using permutation-based significance testing on large datasets.
-   May be moderately influenced by **outliers**, particularly in smaller datasets with high variability.

The Rank Product method is implemented in the Bioconductor package RankProd, which provides functions for performing the analysis and visualizing the results.



### Robust Rank Aggregation (RobustRrankAggreg)

Robust Rank Aggregation (RRA) is a probabilistic method designed to identify genes that are consistently ranked higher than expected by chance across multiple input lists. It is particularly effective in scenarios involving noisy data, incomplete lists, or substantial variability among datasets.

Key features:

-   **Probabilistic modeling**: Computes p-values by modeling the probability of observing gene rankings under a null model.
-   **Robustness to noise and variability**: Maintains performance in the presence of random noise, outliers, or inconsistencies between lists.
-   **Adaptable to varying list lengths**: Can accommodate lists of different sizes without requiring imputation or alignment.
-   **No parameter tuning required**: Offers a straightforward implementation without the need for user-defined parameters.

Recommended use cases:

-   Appropriate for **heterogeneous datasets** obtained from different experimental conditions, platforms, or studies.
-   Particularly suitable when **gene lists are incomplete** or vary significantly in content and length.
-   Scales efficiently with a large number of input lists (more than 20), taking advantage of increased data diversity to improve robustness.
-   Demonstrates high **resistance to noise**, performing reliably even if some input lists contain irrelevant or partially random data.
-   Does not consider the magnitude of expression differences, focusing solely on rank order.
-   Assumes **independence among ranked lists**, which may not always be valid in certain experimental designs.
-   The interpretation of significance scores may be less intuitive due to the probabilistic nature of the method.

The RRA method is implemented in the CRAN package RobustRankAggreg, which offers functions for list aggregation and significance estimation.

### Summary table


| Feature                         | Rank Product                                         | Robust Rank Aggregation                        |
| ------------------------------- | ---------------------------------------------------- | ---------------------------------------------- |
| **Data completeness**           | Requires complete gene presence across lists         | Supports partial and incomplete lists          |
| **List consistency**            | Performs best with uniform list lengths              | Handles varying list lengths and contents      |
| **Handling of missing data**    | Limited unless filtered or penalized                 | Naturally tolerant to missing genes            |
| **Number of input lists**       | Optimal with 5–20 lists                              | Scales well with more than 20 lists            |
| **Noise resistance**            | Moderate                                             | High                                           |
| **Execution time**              | Higher due to permutation testing                    | Lower, computationally efficient               |
| **Quantitative interpretation** | Considers expression magnitude indirectly            | Considers only rank order                      |
| **Recommended applications**    | Datasets with consistent platforms and full coverage | Integration of diverse and incomplete datasets |

<p style="text-align: center; font-style: italic; color: #666; margin-top: 15px; margin-bottom: 10px;">
Table 1: Comparison of meta-analysis methods
</p>


## RankProd Workflow

### Input Methods

MetaRank allows users to input ranked gene lists for Rank Product analysis in three flexible ways:

1.  **Upload Files:**\
    When the “Upload Files” mode is selected, users can upload one or more files in `.TXT`, `.TSV`, or `.CSV` format. Each file represents a ranked gene list from a separate study. The system expects each file to contain at least two columns: a **gene identifier** (e.g., `TP53`) and a **numeric score** representing the expression level, or any other ranking criterion. It is recommended to follow this instructions:

    -   Supported encodings: UTF-8.
    -   Supported delimiters: comma (`,`) or tab (`\t`) (automatically detected).
    -   Do not include headers: remove the corresponding headers for either column names or row names.
    -   Scores are mandatory: if no score are detected, an error will be displayed
    -   Complex gene entries (e.g., `HBA2///HBA1`) are parsed, and only the first gene is retained. (`HBA2`)
    -   NA values, blank lines, and duplicate genes are cleaned automatically.

    Clicking the ℹ️ icon opens a modal window showing the expected file structure and format. There is no strict limit to the number of gene lists that can be uploaded, as it depends on the size of each list. For example, when lists contain approximately 20,000 genes, up to 12 have been successfully processed. In contrast, for smaller lists (ranging from 100 to 500 genes), the system has handled up to 50 lists without issue.


::: callout-note
Many interface elements include contextual tooltips activated by hovering. These tooltips explain each input option, accepted formats, and internal validation steps. For example, hovering over the *Use Example Data* toggle reveals the origin of these datasets, while hovering over the text input area shows how to format pasted genes properly.
:::


2.  **Paste Genes:**\
    When the “Paste Genes” mode is enabled, users can manually paste ranked gene lists into a large text area. This mode supports both `.CSV` and `.TSV` formatting, selectable from a dropdown. Each list must be separated by the string `###`, and within each block, one gene per line is expected. The score must follow the gene, separated by a tab or comma. Unlike the file upload mode, the pasted input must include a header in each block with the exact column names: `Gene` and `Stat.data`:

    ```         
                            Example format (TSV)     Technical format (TSV)
                            Gene    Stat.data        Gene\tStat.data\n
                            TP53    0.95             TP53\t0.95\n
                            BRCA1   0.91             BRCA1\t0.91\n
                            EGFR    0.85             EGFR\t0.85\n
                            ###                      ###
                            Gene    Stat.data        Gene\tStat.data\n
                            BRCA1    0.95            BRCA1\t0.95\n
                            EGFR   0.91              EGFR\t0.91\n
                            TP53    0.85             TP53\t0.85

                            Example format (CSV)     Technical format (CSV)
                            Gene,Stat.data           Gene,Stat.data\n
                            MYC,0.93                 MYC,0.93\n
                            CDK2,0.88                CDK2,0.88\n
                            FOXO1,0.80               FOXO1,0.80n
                            ###                      ###
                            Gene,Stat.data           Gene,Stat.data\n
                            CDK2,0.93                CDK2,0.93\n
                            FOXO1,0.88               FOXO1,0.88\n
                            MYC,0.80                 MYC,0.80
    ```

    Similar to file upload, pasted inputs are automatically cleaned of duplicates and malformed entries. The placeholder text in the input box provides a working example for guidance.

::: {style="display: flex; justify-content: center; gap: 20px; margin-bottom: 1em;"}
<img src="www/images/metarank_rankprod_input.png" alt="MetaRank RankProd input methods" width="265"/> <img src="www/images/metarank_rankprod_info.png" alt="MetaRank RankProd info modal" width="415"/>
:::

<p style="text-align: center; font-style: italic; color: #666; margin-top: -7px;">
Figure 2: a) MetaRank input methods for Rank Product and b) details shown in the ℹ️ info window regarding the required file structure.
</p>

<br>

3.  **Use Example Data:**\
    Enabling the “*Use Example Data*” switch loads four datasets for demonstration purposes. These examples simulate real analysis scenarios with pre-ranked gene lists across multiple studies, allowing users to explore the workflow without providing their own data.

    The table above shows four gene lists used in our example analysis. These lists come from four independent studies related to lung cancer and associated with the following identifiers: *GSE10072, GSE19188, GSE63459, GSE75037*. Each two columns corresponds to a list, containing 22283, 54675, 24526 and 48803 gene identifiers (*SYMBOL*) respectively, including duplicate or missing entries. Each gene has its associated statistical value in the second column (in this case, pvalue). This arrangement allows direct comparison of the size and composition of the lists across studies, highlighting the diverse scope of each dataset prior to subsequent meta-analysis. If the user wishes to study this data in depth, it is possible to download these datasets, as well as view their distribution in the UpsetPlot (*Section 5.1.1.*) and Heatmap (*Section 5.1.2.*). 

<br>

```{r, echo=FALSE, warning=FALSE, message=FALSE}
library(DT)

example_files <- c(
  "./example_data/data1.tsv",
  "./example_data/data2.tsv",
  "./example_data/data3.tsv",
  "./example_data/data4.tsv")

gene_tables <- lapply(example_files, function(path) {
  read.delim(path, header = FALSE, stringsAsFactors = FALSE, col.names = c("Gene", "Score"))})
max_rows <- max(sapply(gene_tables, nrow))
pad_df <- function(df, max_rows) {
  n <- nrow(df)
  if (n < max_rows) {
    pad <- data.frame(Gene = rep(NA, max_rows - n), Score = rep(NA, max_rows - n), stringsAsFactors = FALSE)
    df <- rbind(df, pad)} 
  df}

gene_tables_padded <- lapply(gene_tables, pad_df, max_rows = max_rows)
result_df <- data.frame(matrix(NA, nrow = max_rows, ncol = length(gene_tables_padded)*2))
colnames(result_df) <- unlist(lapply(seq_along(gene_tables_padded), function(i) {
  c(paste0("Gene File", i), paste0("Score File", i))}))
for (i in seq_along(gene_tables_padded)) {
  df <- gene_tables_padded[[i]]
  result_df[[paste0("Gene File", i)]] <- df$Gene
  result_df[[paste0("Score File", i)]] <- df$Score}

DT::datatable(result_df, rownames = FALSE, options = list(pageLength = 10))
```


<p style="text-align: center; font-style: italic; color: #666; margin-top: 20px;">
Table 2: Sample data content.
</p>

<br>

### Parameters

Once the gene lists are loaded, six configuration parameters become available to customize the RankProd analysis. These options provide full control over how genes are filtered and ranked. You can fine-tune aspects such as ranking direction, handling of missing values, penalization of genes with low recurrence, and the minimum number of lists a gene must appear in to be considered. This flexibility ensures the meta-analysis is aligned with your experimental design and data quality.

![](www/images/metarank_rankprod_parameters.png){fig-align="center" width="225"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: -7px;">
Figure 3: The parameter panel offered by MetaRank to adjust the meta-analysis.
</p>

<br>

#### Rank-based Method

A rank-based meta-analysis combines gene rankings across multiple studies or conditions instead of directly comparing raw values. This approach is especially useful when datasets are heterogeneous or measured on different scales, allowing robust integration based on gene order rather than absolute expression.

There are two available modes:

-   **Basic**: Uses `RankProd::RP`. It assumes that each gene list comes from a unique origin (i.e., no shared batches). This is ideal when the true origin of your data is unknown or when you prefer not to group them explicitly. Suitable for datasets with unknown or homogeneous background (e.g., mixed public data without batch labels).

-   **Advanced**: Uses `RankProd::RP.Advanced`. This mode allows specifying a vector of origins (or batches) for each list via the *Origin* field (explained in the section 3.2.2.). Recommended when your gene lists come from distinct experimental setups, platforms, conditions, or time points. It adjusts the ranking by grouping lists with the same origin, improving robustness in multi-batch scenarios.

#### Origin (*Advanced only*)

The *Origin* field is required when using the `Advanced` mode. It should be a comma-separated vector of integers (e.g., `1,1,2,2`), where each number indicates the batch or origin of the corresponding input file.

This field:

-   Must match the number of input gene lists.
-   Allows grouping lists from the same source.
-   Each batch must have replicas, i.e. at least two datasets from the same source.
-   Is validated with custom error messages if the format is incorrect or inconsistent.
-   Is accompanied by an ℹ️ info button with a usage example for user guidance.

![](www/images/metarank_rankprod_origin.png){fig-align="center" width="456"}


<p style="text-align: center; font-style: italic; color: #666; margin-top: -7px;">
Figure 4: Details displayed in the information window ℹ️ about the structure of the Origin variable 
</p>

<br>


#### Minimum Number of Datasets

This slider sets the minimum fraction of input lists in which a gene must appear to be included in the analysis. Genes present in fewer lists will be excluded before the ranking process. For example: If you set it to 4 and there are 5 input lists, only genes appearing in at least 4 lists will be considered (4 and 5).

This filter helps reduce statistical noise caused by infrequent genes that may distort the consensus ranking. Genes that appear only once are always shown separately in the "excluded genes" table, since they do not allow robust comparison and may bias the analysis if included.

#### Ranking Direction

This option indicates whether lower or higher values should be considered better rankings. This depends on the type of metric:

- `Ascending`: lower values are better (e.g., pvalues).
- `Descending`: higher values are better (e.g., logFC, z-scores, relevance scores).

It is important to choose the correct direction to ensure proper interpretation of the results.

#### NA Management

This option determines how to handle missing values (genes not present in some of the lists):

| Option        | Description                                                                                                       |
|---------------|-------------------------------------------------------------------------------------------------------------------|
| Impute NA     | Replaces NA with the median rank of the list. Allows applying an extra penalty based on the number of appearances. Useful when preserving all genes and reducing the impact of missing values, for example in exploratory analyses. |
| Ignore NA     | Uses only the available values, omitting NAs. Also allows extra penalization based on the number of appearances. Useful when preserving all genes and reducing the impact of missing values, for example in exploratory analyses.   |
| Penalize NA   | Assigns the worst possible rank to missing values, depending on whether the direction is ascending or descending. Recommended when missing values should be heavily penalized to increase robustness.  |


<p style="text-align: center; font-style: italic; color: #666; margin-top: 20px;">
Table 3: Missing values managment options.
</p>


::: callout-note

If you apply a `Minimum Number of Datasets` filter that requires genes to appear in all lists, there will be no missing values and this setting will have no effect. Note that the impact of missing value management depends on how strict or tolerant the user wants the analysis to be. More relaxed settings retain more genes but may introduce noise, while stricter settings increase reliability but may discard potentially relevant genes.

:::


#### Extra Penalization (*Impute & Ignore only*)

When this option is enabled, an additional penalty is applied to each gene depending on the number of lists in which it appears (calculated before the analysis, but applied after it). The fewer times a gene appears, the worse its adjusted ranking will be, even if it initially ranked well. An adjusted rank is calculated by adding a penalty proportional to the number of lists where the gene is missing.

Conceptual formula:

```
          AdjustedRank = Rank + ((TotalLists - Count) * (MaxRank / TotalLists))
```

Where:

- `Rank`: the original consensus rank of the gene.
- `TotalLists`: the total number of gene lists loaded.
- `Count`: the number of lists in which the gene appears.
- `MaxRank`: the worst (highest) rank in the current ranking.

This adjustment is especially useful when using the Impute or Ignore NA options, to ensure that genes with limited support across datasets are penalized accordingly and do not dominate the top of the consensus ranking. This helps to prioritize genes that are consistently present and reduce the impact of rare, potentially spurious genes.


### Outputs

#### Results Table (RankProd)

The main output of the RankProd analysis is a table with the following columns and their meanings:

| Column Name    | Description                                                                                                                                              |
|----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------|
| **GeneID**     | Unique gene identifier, which can be a HUGO symbol (e.g., TP53), an Entrez ID (e.g., 7157), or an Ensembl ID (e.g., ENSG00000141510), depending on input. |
| **Rank**       | Consensus ranking of the gene across all input lists; lower values indicate higher overall relevance or consistency among the lists.                      |
| **FileCount**  | Number of input gene lists in which this gene appears; a higher count suggests greater consistency across datasets.                                       |
| **FileNames**  | Names of the input files where the gene was found, separated by spaces; useful for identifying the sources supporting the gene's relevance.              |
| **GenePositions** | Rank positions of the gene in each individual input list; provides insight into the gene's performance across different datasets.                     |
| **RP_stat**    | Rank Product statistic calculated to assess the significance of the gene's ranking across multiple lists; lower values suggest higher significance.       |
| **PFP**        | Estimated Proportion of False Positives, analogous to False Discovery Rate (FDR); lower values indicate more reliable findings.                           |
| **pvalue**     | Raw p-value from the meta-analysis, indicating the probability of observing the gene's ranking by chance; lower values suggest higher significance.        |
| **p.adjust**   | Adjusted p-value accounting for multiple hypothesis testing using the Benjamini-Hochberg method; helps control the FDR.                                  |

<p style="text-align: center; font-style: italic; color: #666; margin-top: 20px;">
Table 4: The names of the columns and their meanings in the main table generated by RankProd in MetaRank.
</p>


The table can be downloaded in `.TSV` and `.CSV` formats. Each column has a tooltip (question mark icon) that shows this information when hovered over. The table is interactive: columns can be filtered by value ranges or keywords, and their order can be customized.

::: callout-tip

If data filtering is applied, either by value or by selecting specific columns (see Section 5.3.1), the downloaded file will reflect only the currently displayed data. This includes both the filtered rows and the visible columns selected in the user interface.

:::


#### Excluded Genes

A secondary table is generated and accessible via the eye icon button, also downloadable as `.TSV`. It always contains genes excluded by the `Minimum Number of Datasets` filter and those appearing only once. If, for example, a filter of 4/4 is applied, genes appearing in 1, 2, or 3 lists are moved to this excluded table, while only genes appearing in all 4 lists remain in the main table. This table allows tracking of excluded genes and understanding of filtering effects.


| Column Name   | Description                                                                                                   |
|---------------|---------------------------------------------------------------------------------------------------------------|
| **GeneID**    | Unique identifier of the excluded gene, in the same format as the input.                                      |
| **FileCount** | Number of input gene lists in which this gene appears.                                                        |
| **FileNames** | Names of the input files where the gene was found.                                                            |

<p style="text-align: center; font-style: italic; color: #666; margin-top: 20px;">
Table 5: The names of the columns and their meanings in the excluded table generated by RankProd in MetaRank.
</p>

<br>


## RobustRankAggreg Workflow

### Input Methods

1.  **Upload Files:** When the “Upload Files” mode is enabled, it is possible to select one or more text files (`.TXT`) via the file upload control. The system recognizes each file as a list of genes (one identifier per line, without a header), automatically removes duplicates and missing values, and correctly handles both Unix (`\n`) and Windows (`\r\n`) line endings. If more than one gene is provided on a single line separated by delimiters (e.g., `BRCA1///BRCA2`), only the first entry (`BRCA1`) is retained. 

    Clicking the ℹ️ icon opens a modal showing a sample file structure, and example datasets can be downloaded for in-depth study and reference. There is no strict limit to the number of gene lists that can be uploaded, as it depends on the size of each list. For example, when lists contain approximately 20,000 genes, up to 12 have been successfully processed. In contrast, for smaller lists (ranging from 100 to 500 genes), the system has handled up to 50 lists without issue.

2.  **Paste Genes:** When the “Paste Genes” mode is enabled, gene lists can be entered directly into a text area. Each list is delimited by `###`, and within each section the system expects one gene per line, with no header row. 

    ```         
                              Example format    Technical format
                              TP53              TP53\n
                              BRCA1             BRCA1\n
                              EGFR              EGFR\n
                              ###               ###\n
                              BRCA1             BRCA1\n
                              EGFR              EGFR\n
                              TP53              TP53
    ```

    Duplicate entries and blank lines are cleaned up automatically, and if a line contains multiple gene identifiers (e.g., `BRCA1///BRCA2`), only the first is used. The placeholder text illustrates this formatting.

::: {style="display: flex; justify-content: center; gap: 20px; margin-bottom: 1em;"}
<img src="www/images/metarank_rra_input.png" alt="MetaRank RankProd input methods" width="265"/> <img src="www/images/metarank_rra_info.png" alt="MetaRank RankProd info modal" width="415"/>
:::

<p style="text-align: center; font-style: italic; color: #666; margin-top: -7px;">
Figure 5: a) MetaRank input methods for Robust Rank Aggreg and b) details shown in the ℹ️ info window regarding the required file structure.
</p>


3.  **Use Example Data**: Enabling the “*Use Example Data*” switch loads predefined files that represent various analysis scenarios, allowing users to explore the workflow without providing their own data.

    The table above shows four gene lists used in our example analysis. These lists come from four independent studies related to lung cancer and associated with the following identifiers: GSE10072, GSE19188, GSE63459, GSE75037. Each column corresponds to a list, containing 19417, 21752, 17509 and 13099 gene identifiers (*SYMBOL*) respectively, including duplicate or missing entries. This arrangement allows direct comparison of the size and composition of the lists across studies, highlighting the diverse scope of each dataset prior to subsequent meta-analysis.

```{r, echo=FALSE, warning=FALSE, message=FALSE}
library(DT)
example_files <- c(
  "./example_data/gene1.txt",
  "./example_data/gene2.txt",
  "./example_data/gene3.txt",
  "./example_data/gene4.txt")

gene_lists <- lapply(example_files, function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines[lines != ""]})
max_length <- max(sapply(gene_lists, length))
gene_lists_padded <- lapply(gene_lists, function(vec) {
  length(vec) <- max_length
  vec})
df_gene_lists <- as.data.frame(
  setNames(gene_lists_padded, paste0("List_", seq_along(gene_lists_padded))),
  stringsAsFactors = FALSE)

datatable(df_gene_lists, rownames = FALSE)
```

<p style="text-align: center; font-style: italic; color: #666; margin-top: 20px;">
Table 6: The content of the 4 example files together in one table.
</p>

<br>

### Parameters

Once the gene lists are loaded, three configuration parameters become available to customize the `RobustRankAggreg` analysis. These options provide some control over how genes are filtered and ranked. You can fine-tune aspects such as selecting the aggregation method, handling of missing values, or even the minimum number of lists a gene must appear in to be considered. This flexibility ensures the meta-analysis is aligned with your experimental design and data quality.

![](www/images/metarank_rra_parameters.png){fig-align="center" width="225"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: -7px;">
Figure 6: The parameter panel offered by RobustRankAggreg to adjust the meta-analysis. 
</p>

#### Aggregation Method

The RobustRankAggreg package offers five aggregation methods to combine rankings across multiple gene lists, ranging from simple statistical approaches to the more sophisticated probabilistic scoring native to the package:

-   **RRA**: Uses a probabilistic model to assign p-values to ranks, based on the minimum probability across all lists. It evaluates how surprising a gene's ranking is across the datasets using a beta-uniform mixture model.
-   **Median**: Takes the median rank of each gene across all lists. It is robust to outliers and provides a central tendency measure.
-   **Stuart**: A method based on order statistics. It combines ranks using a meta-analysis approach, particularly suitable for independent rankings.
-   **Geometric Mean**: Computes the geometric mean of ranks across lists, giving more weight to consistently low ranks.
-   **Arithmetic Mean**: Averages the rank values directly. This method is sensitive to outliers but intuitive and easy to interpret.

All of these methods rely on the position of genes in the individual rankings to compute a consensus order. The `RRA method` is unique in that it transforms rankings into *p-values* and evaluates their statistical significance, accounting for both the number of lists and the positions within each.

#### Minimum Number of Datasets

This slider sets the minimum fraction of input lists in which a gene must appear to be included in the analysis. Genes present in fewer lists will be excluded before the ranking process. For example: If you set it to 4 and there are 5 input lists, only genes appearing in at least 4 lists will be considered.

This filter helps reduce statistical noise caused by infrequent genes that may distort the consensus ranking. Genes that appear only once are always shown separately in the "excluded genes" table, since they do not allow robust comparison and may bias the analysis if included.

#### NA Management

This option controls how to handle missing values (i.e., when a gene does not appear in a list). Two strategies are provided:

-   **Ignore NA**: Exclude missing values from the analysis.
-   **Penalize NA**: Assign worst rank for missing entries

In this context, additional penalization is not required beyond what the algorithm already incorporates. The score, also known as rho, is a significance measure used in RobustRankAggreg to reflect how strongly a gene is supported across the rankings. It is based on the minimum p-value method:

-   Each gene’s position in a list is converted to a probability. For example, if a gene ranks 5th in a list of 1000, we calculate the probability of randomly selecting a gene ranked 5th or better.
-   The lowest (best) of these probabilities across all lists is taken.
-   The final score is calculated using a beta-uniform distribution, estimating the likelihood of observing such a good ranking by chance, given how many lists exist and in how many the gene appears.

If a gene is absent from some lists, the method does not assign an artificially bad rank. Instead, the score inherently adjusts for the fact that a gene appeared in fewer lists. This naturally penalizes low-frequency genes unless they show extremely strong evidence in the lists they do appear in. A low score means the gene’s strong ranks are unlikely to be due to chance and that it is consistently important across studies.

### Outputs

#### Results Table (RRA)

The output table from the RobustRankAggreg (RRA) workflow differs slightly from the one used in RankProd. The available columns are:

| Column Name        | Description                                                                 |
|---------------|-----------------------------------------------------------------------------|
| **GeneID**    | Unique gene identifier, which can be a HUGO symbol (e.g., TP53), an Entrez ID (e.g., 7157), or an Ensembl ID (e.g., ENSG00000141510), depending on input.                          |
| **Rank**      | Consensus ranking of the gene across all input lists; lower values indicate higher overall relevance or consistency among the lists.       |
| **Score**     | Also called *rho*, this is the probabilistic score assigned by RRA, reflecting the significance of the observed ranks (lower values indicate stronger evidence). |
| **p.adjust**  | Adjusted p-value (multiple testing correction) for the Score using the Benjamini-Hochberg method; helps control the FDR.              |
| **FileCount** | Number of input gene lists in which this gene appears; a higher count suggests greater consistency across datasets.                          |
| **FileNames** | Names of the input files where the gene was found, separated by spaces; useful for identifying the sources supporting the gene’s relevance.                      |
| **GenePositions** | Rank positions of the gene in each individual input list; provides insight into the gene’s performance across different datasets.             |

<p style="text-align: center; font-style: italic; color: #666; margin-top: 20px;">
Table 7: The names of the columns and their meanings in the main table generated by RobustRankAggreg in MetaRank.
</p>

The table can be downloaded in TSV and CSV formats. Each column has a tooltip (question mark icon) that shows information when hovered over. The table is interactive: columns can be filtered by intervals or keywords and reordered.

#### Excluded Genes

A secondary table is generated and accessible via the eye icon button, also downloadable as TSV. It always contains genes excluded by the `Minimum Number of Datasets` filter and those appearing only once. If, for example, a filter of 4/4 is applied, genes appearing in 1, 2, or 3 lists are moved to this excluded table, while only genes appearing in all 4 lists remain in the main table.


| Column Name   | Description                                                                                                   |
|---------------|---------------------------------------------------------------------------------------------------------------|
| **GeneID**    | Unique identifier of the excluded gene, in the same format as the input.                                      |
| **FileCount** | Number of input gene lists in which this gene appears.                                                        |
| **FileNames** | Names of the input files where the gene was found.                                                            |

<p style="text-align: center; font-style: italic; color: #666; margin-top: 20px;">
Table 8: The names of the columns and their meanings in the excluded table generated by RobustRankAggreg in MetaRank.
</p>


This table allows tracking of excluded genes and understanding of filtering effects.

<br>



## Shared Elements

While each method, `RankProd` or `RobustRankAggreg (RRA)`, has its own specific parameters and analysis pipeline, the app also includes a set of shared components that remain available regardless of the selected method. Among the shared functionalities, the app includes:  

- A set of raw data visualizations (UpSet plot and interactive Heatmap) to examine overlaps and differences between gene lists.  
- An enrichment analysis module based on the top 100 consensus genes obtained from the ranking process, with independent outputs.  
- A visualization settings panel, allowing customization of the appearance of plots and tables (such as number of terms shown, font sizes, and colors).

These features provide robust tools for evaluating gene list consistency, biological relevance, and presentation quality of the results.

### Shared Plots

Regardless of the selected method `RankProd` or `RobustRankAggreg (RRA)` the output always includes two additional visualizations: an **UpSet plot** and a **Heatmap**. These plots display the distribution of the raw input data, helping users to understand the relationships among the gene lists before any ranking or aggregation is performed. They allow checking for common genes across lists, unique genes in each list, and their proportions.

#### UpSet Plot

The UpSet plot visualizes intersections between the input gene lists. The top panel shows the size of each intersection (i.e., how many genes are shared between specific combinations of lists), while the left panel shows the size of each individual list. This plot is particularly useful when dealing with multiple sets where traditional Venn diagrams become difficult to interpret. The features of this plot are:

-   The UpSet plot can be downloaded as PNG or JPG.
-   Users can customize the colors of both horizontal and vertical bars.
-   Text size (including axis labels, titles, and legends) can also be adjusted.
-   For more details, see section 5.2 Data Visualization Tab.

**Interpretation:**  
- Tall bars in the top panel indicate large overlaps between specific sets of gene lists.  
- The connected dots below the bars specify which lists are involved in each intersection.  
- This allows quick identification of genes common to many lists or unique to one or more specific lists.

![](www/images/metarank_upsetplot.jpg){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: -7px;">
Figure 7: Upset Plot generated with the example data.
</p>


#### Heatmap

The heatmap visualizes the pairwise similarity between all input gene lists. Each cell represents the proportion of shared genes between two lists, calculated as the ratio of common genes to the total number of genes in the corresponding reference list.

-   The diagonal cells always show a value of 1.00, as each list is identical to itself.
-   Off-diagonal cells indicate the degree of overlap between different lists (e.g., a value of 0.59 between List 1 and List 3 means that 59% of genes in List 1 are also present in List 3).
-   The comparison is based on raw gene content, without considering rank or associated statistics.

**Features**: 

-   The heatmap can be downloaded in PNG, JPG, or HTML formats.
-   Users can adjust the title size, axis label size, and select among different color scales (e.g., Viridis, Cividis, Portland) to enhance readability and presentation (for more details, refer to section 5.2 Data Visualization Tab).

**Interpretation:**  

- High similarity values indicate strong agreement in gene presence between lists.
- Lower values suggest variability or dataset-specific gene composition.  
- This visualization helps identify outlier datasets or assess overall consistency among inputs.

<iframe src="/www/images/metarank_heatmap.html" style="width:100%; height:500px;" frameborder="0">

</iframe>

<p style="text-align: center; font-style: italic; color: #666; margin-top: -7px;">
Figure 8: Heatmap generated with the example data.
</p>

::: callout-note  

This heatmap is fully interactive. Users can zoom in on specific regions, hover over individual cells to see detailed information (such as gene ID, list name, and rank/value), and explore patterns in greater detail. This interactivity enhances the ability to identify key trends and outliers within the data.

:::

<br>

### Shared Enrichment Analysis

Once the consensus ranking has been generated (whether by `RankProd` or `RobustRankAggreg`) MetaRank offers the option to perform *Over-Representation Analysis (ORA)* using the top-ranked genes. This analysis helps identify biological terms or pathways that are significantly associated with the consensus gene list, providing biological context and insight into the aggregated results.

Up to the top 100 genes from the consensus list can be used for enrichment, offering a flexible way to explore functional relevance across methods.

![](www/images/metarank_enrich.png){fig-align="center" width="200"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: -7px;">
Figure 9: Enrichment analysis parameters at MetaRank.
</p>


#### Number of genes

A slider allows users to choose the number of top genes used in the enrichment analysis, ranging from 10 to 100 genes, increasing by steps of 1. For example, selecting exactly 41 genes is possible.

::: callout-tip

If the enrichment result shows no significant terms (e.g., the message “no biological information found” appears), consider increasing the number of genes selected. A larger input list increases the chance of capturing enriched pathways or categories.

:::

#### Database

-   **Gene Ontology (GO)**: A structured and controlled vocabulary used to describe the functions of genes and their products in a consistent and standardized way. Its content is divided into three sub-ontologies:

    -   **Biological Process (BP)**: Pathways and larger processes (e.g., *cell cycle, signal transduction*).
    -   **Molecular Function (MF)**: Biochemical activities (e.g., *ATP binding, kinase activity*).
    -   **Cellular Component (CC)**: Subcellular locations (e.g., *nucleus, ribosome*).

-   **KEGG (Kyoto Encyclopedia of Genes and Genomes)**: Database resource for understanding high-level functions and utilities of the biological system, such as the cell, the organism and the ecosystem, from molecular-level information, especially large-scale molecular datasets generated by genome sequencing and other high-throughput experimental technologies.

-   **Reactome**: A curated database of human biological pathways and reactions, including signaling, metabolism, and immune system processes. It also supports several model organisms via ortholog mapping.

All three resources serve as the reference background in the `Over-Representation Analysis (ORA)` step, where the frequency of your input genes in each term or pathway is statistically compared against a genomic background to identify the most significantly enriched biological categories.

#### Onology (*GO only*)

When *Gene Ontology (GO)* is selected as the database, an Ontology choice must be specified. Each GO branch provides a distinct view of gene function:

-   **Biological Process (BP)** Describes high-level biological objectives accomplished by ordered assemblies of molecular functions—such as “*cell cycle*,” “*signal transduction*,” or “*immune response*.” Use BP to discover which overarching pathways or processes your genes collectively influence.
-   **Molecular Function (MF)** Captures the elemental activities of proteins or gene products at the biochemical level—examples include “*ATP binding*,” “*kinase activity*,” or “*transcription factor binding*.” MF is ideal for pinpointing the specific enzymatic or binding roles enriched in your gene set.
-   **Cellular Component (CC)** Defines where gene products exert their function within the cell, such as “*nucleus*,” “*mitochondrion*,” or “*ribosome*.” CC helps reveal the subcellular localization patterns common to your genes, indicating, for instance, whether they cluster in particular organelles.

Selecting the appropriate ontology refines the enrichment analysis, focusing it on either broader process-level insights (BP), detailed activity-level functions (MF), or spatial context within the cell (CC).

#### Organism

The following organisms are supported, each identified by an internal code and NCBI Taxonomy ID:

-   **Homo sapiens** (`Hsa`; Taxonomy ID: 9606) – Human gene symbols follow the HGNC standard (e.g., `TP53`).
-   **Mus musculus** (`Mmu`; Taxonomy ID: 10090) – Mouse gene symbols use the MGI nomenclature (e.g., `Trp53`).
-   **Rattus norvegicus** (`Rno`; Taxonomy ID: 10116) – Rat gene symbols use RGD notation (e.g., `Rps6kb1`).

#### GeneID

Specify the identifier system for your gene lists:

-   **SYMBOL**: Common gene names or symbols (e.g., `BRCA1`).
-   **ENTREZID**: Unique numerical IDs assigned by NCBI (e.g., `672`).
-   **ENSEMBL**: Stable gene IDs from the Ensembl database (e.g., `ENSG00000012048`).

::: callout-warning
It is crucial to match your gene list format to the selected ID type. Selecting the wrong identifier system will lead to failed mappings and inaccurate enrichment results.
:::

#### Outputs

##### Enrichment Table

The enrichment table displays the results of the Over-Representation Analysis (ORA) performed on the consensus gene list. It contains the following columns and can be downloaded as CSV or TSV files for further use:

| Column          | Description                                                                                                                                                                        |
| --------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **ID**          | The unique identifier of the enrichment term. Examples include *GO:0006915* for “apoptotic process” (Gene Ontology) or *05200* for “Pathways in cancer” (KEGG).                    |
| **Description** | A concise, human-readable name or description of the biological term or pathway, providing clear context.                                                                          |
| **GeneRatio**   | The ratio of input genes associated with the term over the total number of genes selected for enrichment (e.g., 5/100).                                                            |
| **BgRatio**     | The ratio of all genes annotated with the term in the background/reference genome (e.g., 200/20000).                                                                               |
| **pvalue**      | Raw p-value from the enrichment test (e.g., Fisher’s exact test), representing the probability of observing the enrichment by chance. Lower values indicate stronger significance. |
| **p.adjust**    | Adjusted p-value after multiple testing correction using the Benjamini–Hochberg False Discovery Rate (FDR) method. Lower values (< 0.05) suggest reliable significance.            |
| **GeneCount**   | Number of genes from the input list associated with the enrichment term.                                                                                                           |
| **GeneID**      | List of gene identifiers (e.g., SYMBOL, ENTREZID, or ENSEMBL) from the input that contributed to the enrichment, helping to identify specific genes driving the signal.            |


<p style="text-align: center; font-style: italic; color: #666; margin-top: 20px;">
Table 9: The names of the columns and their meanings in the enrichment table generated by MetaRank.
</p>


::: callout-tip
Each column in the results table includes a small question mark icon that provides additional information when hovered over. These tooltips offer concise explanations of the column’s purpose and how to interpret its values, allowing users to quickly understand the data without needing to refer back to the documentation. Simply place your cursor over the icon to view the tooltip — no need to click.
:::




##### Enrichment Plot

The data resulting from the meta-analysis and summarized in the main results table is also represented visually through interactive plots. These plots provide a complementary overview of enriched terms and support intuitive interpretation of the results with this following features:

-   Two plot type are available: *Dot Plot* and *Bar Plot*, both depicting the number of genes associated with each term (gene count) and their corresponding adjusted p-values (`p.adjust`). The color scale can be customized to reflect statistical significance.
-   Tooltips appear on hover, displaying key information such as the full term or pathway name, GeneCount, and adjusted p-value, enhancing interpretability.
-   Download Options include `PNG`, `JPG`, or `HTML` (interactive) formats, allowing users to export plots for presentations or further analysis.

::: callout-note
Table-Plot reactivity ensures consistency. For example, filtering the results table (e.g., by the word "cell") dynamically updates the plot to show only the matching terms.
:::

::: panel-tabset
###### Dotplot

![](www/images/metarank_dot.png){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: -7px;">
Figure 10: Dot plot generated to represent the enrichment analysis result.
</p>

###### Barplot

![](www/images/metarank_plot.png){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: 12px;">
Figure 11: Bar plot generated to represent the enrichment analysis result.
</p>

###### Interactive

<iframe src="/www/images/metarank_bar_html.html" style="width:100%; height:500px;" frameborder="0">

</iframe>

<p style="text-align: center; font-style: italic; color: #666; margin-top: -7px;">
Figure 12: Interactive bar plot generated to represent the enrichment analysis result (with only term ID selected to improve the visualisation).
</p>

:::

<br>

### Data Visualization Tab

The Data Visualization tab, located in the sidebar panel, provides several customization options to tailor both the appearance and content of the output tables and plots. These settings are useful for enhancing presentation clarity, adjusting for accessibility needs, or focusing on specific aspects of the analysis.

#### Meta-analysis Table Settings

Users can toggle the visibility of specific columns in the consensus ranking table. The available columns differ depending on the selected meta-analysis method:

-   RankProd includes columns like: `GeneID`, `Rank`, `FileCount`, `FileNames`, `GenePositions`, `RP_stat`, `PFP`, `pvalue`, and `p.adjust`.
-   RRA includes: `GeneID`, `Rank`, `Score`, `p.adjust`, `FileCount`, `FileNames`, and `GenePositions`.

Only selected columns are shown in the table and included in the downloaded file. This allows exporting customized summaries focused on the user's needs.

![](www/images/metarank_data_1.png){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: 12px;">
Figure 13: Main table column selection.
</p>

#### Upset plot

Customization options include:

-   **Text size**: Affects axis labels, title, and legend. Recommended value: 27px.
-   **Color customization**: Separate color selectors for vertical and horizontal bars. Default colors are green (`#2ba915`) and blue (`#0838a0`), respectively.

These options help adapt the plot for presentations or specific visual preferences.

![](www/images/metarank_data_2.png){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: 12px;">
Figure 14: UpsetPlot customisation.
</p>


#### Heatmap plot

This section provides options to:

-   Adjust axis text size (recommended: 16px).
-   Adjust title size (recommended: 20px).
-   Select a color scale: Choose among *Viridis, Cividis*, or *Portland* palettes to match different visual or accessibility requirements.

![](www/images/metarank_data_3.png){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: 12px;">
Figure 15: Heatmap customisation.
</p>


#### Enrichment Table Settings

The interface allows selecting and deselecting specific columns from the enrichment result table. This customization enables users to tailor the displayed information to their needs. When using the download buttons, only the currently visible columns will be included in the exported file. For instance, if only the `ID` and `Description` columns are selected out of eight possible ones, the downloaded table will contain just those two.

![](www/images/metarank_data_4.png){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: 12px;">
Figure 16: Enrichment table column selection.
</p>


::: callout-note
-   Filters applied directly within the table interface (e.g., keyword searches such as "cell") are reflected in the downloaded file too, but it works after running the analysis
-   The table listing excluded genes or terms is not customizable—its structure and contents remain fixed for both display and export (full table).
:::

#### Enrichment Plot Settings

This section allows customization of the appearance of enrichment plots to better suit presentation or analysis needs. Several visual parameters can be adjusted:

-   **Number of terms to show**: Defines how many of the top-ranking enriched terms will be displayed in the plot. This helps focus the visualization on the most relevant results.
-   **Y-axis**: Users can choose what appears on the Y-axis—either the `Term ID`, the `Description`, or a combination of both. When selecting both, the term and its description are shown together, separated by a hyphen (“`-`”), providing more context for each entry.
-   **Plot Type**: Two types of plots are supported:
    -   *Dot plot*, which represents each term as a point, usually with size or color indicating significance or gene count.
    -   *Bar plot*, where each term is displayed as a bar, useful for comparing absolute or relative values.
-   **Color Scale**: A gradient color scale is applied based on the adjusted p-value (`p.adjust`). Users can customize both ends of the color gradient (low and high values) to match their preferred visual style or color scheme. This makes interpretation easier, especially when using consistent color themes across multiple plots.
-   **Text Size**: Allows control over the font size used in the plot. Increasing or decreasing this value can help adapt the plot for screens, print, or accessibility preferences.

![](www/images/metarank_data_5.png){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: 12px;">
Figure 17: Enrichment plot customisation.
</p>

<br>

## Background Pipeline

MetaRank performs consensus-based gene ranking using two complementary strategies: `RankProd (RP)` and `RobustRankAggreg (RRA)`. The workflow processes input gene lists through a series of well-defined steps to produce both tabular and graphical enrichment results:

### Analysis Method Selection

Users can select one of the following meta-ranking algorithms:

- **RankProd (RP)**:
  - A non-parametric method for identifying genes that are consistently up- or down-regulated across experiments.
  - Relies on ranking statistics typically derived from fold-change values.
  - Its functionality has been adapted to improve robustness in consensus ranking, including enhanced handling of `NA` values during computation.

- **Robust Rank Aggregation (RRA)**:
  - Detects genes that consistently appear at the top of ranked lists more often than expected by chance.
  - Operates solely on the order of ranks, without requiring fold-change values or p-values.
  - Supports multiple ranking strategies, from simple methods like the median rank to more sophisticated ones such as the original RRA algorithm.


### Data Ingestion and Preprocessing

Users can provide gene data by uploading multiple files (`.CSV`, `.TSV`, `.TXT`) or pasting lists directly into a text box (supports multiline input). Additionally, it is possible to work exclusively with the example datasets provided by the app, which are also available for download.

The required input format depends on the selected analysis package:

- **RankProd**:
  - Accepts `.CSV` (comma-separated) or `.TSV` (tab-separated) files.
  - Each file or list must contain **two columns**:
    - One with gene identifiers (e.g., `Gene`, `EntrezID`, or `Ensembl`).
    - One numeric column representing a ranking metric (e.g., `pvalue`, `logFC`, etc.).
  - Uploaded files **must NOT** contain headers. In contrast, when using the paste mode, each list **must include headers**, specifically `Gene` and `Stat.data`.
  - In paste mode, multiple lists must be separated using the delimiter `###`.
  - Make sure to check the **info (ℹ) button** for a detailed explanation of the correct input format.
  
- **RRA**:
  - Accepts `.TXT` files or pasted plain text.
  - Each file or list must consist of **a single column** containing gene identifiers (e.g., `Gene`, `EntrezID`, or `Ensembl`), listed one per line in descending order of significance.
  - In paste mode, multiple lists must be separated using the delimiter `###`.
  - Make sure to check the **info (ℹ) button** for a detailed explanation of the correct input format.

- **Validation**:
  - The app validates uploaded content to ensure consistent formatting, presence of required columns, and proper delimiters.
  - Automatic preprocessing includes:
    - Removing blank rows.
    - Trimming leading/trailing whitespace.
    - Attempting to detect the gene identifier format (SYMBOL, ENSEMBL, or ENTREZID).
  - If invalid input is detected, the app shows informative modals and provides example formats to guide the user.


### Gene Appearance Counting and Filtering

- For every gene across **all input lists**, the system:
  - **Counts the total number of appearances** (i.e., in how many lists the gene is found),
  - **Records the names** of the lists (or files) in which the gene appears,
  - **Stores the rank positions** of the gene in each list where it is present. For example, if a gene appears at position 45 in list 1, position 98053 in list 2, and position 1 in list 3, the position vector would be: `45, 98053, 1`.

- This information is compiled into a detailed **appearance table** for each gene, enabling complete traceability and data auditing.

- A **user-defined minimum appearance threshold** is applied:
  - Genes must appear in a minimum number of input lists to be included in the final meta-analysis.
  - This filtering step removes **low-frequency or list-specific genes**, which helps reduce background noise and increases the robustness of the consensus ranking.
  - The threshold is configurable by the user to balance inclusiveness and specificity.

- **Outputs**:
  - **Included genes**: Genes that meet or exceed the appearance threshold. These are used in the meta-ranking process and included in the final results.
  - **Excluded genes**: Genes that do **not** meet the threshold. These are completely **excluded from both the consensus ranking and any enrichment analysis**, which also helps reduce computational time. They are still accessible for **review and optional download**.


### Consensus Ranking Computation

The selected meta-ranking method is applied to the filtered gene lists:

- **RankProd**:
  - Computes rankings for **upregulated** and **downregulated** genes separately.
  - Utilizes the `RP` or `RP.advance` functions from the `RankProd` Bioconductor package.
  - Includes advanced options such as:
    - **Handling of missing values (`NA`)** gracefully.
    - **Custom directionality** settings to rank by high or low values depending on the metric.
    - **Optional penalization** of genes with low appearance frequency to further refine the consensus.

- **RRA**:
  - Aggregates multiple ranked gene lists into a single consensus ranking.
  - Uses the `aggregateRanks` function from the `RobustRankAggreg` package.
  - Performs a permutation-based statistical analysis to calculate:
    - **P-values**, representing the likelihood of observing such high rankings by chance.
    - **Adjusted p-values**, corrected for multiple testing using standard methods (Benjamini-Hochberg).
    
    
### Final table creation

After computing the consensus ranking, a final result table is generated by merging the ranking outputs with the detailed gene appearance data.

- For **included genes**, the final table contains:
  - Consensus ranking metrics (e.g., rank, p-value, score depending on method).
  - The number of appearances across all input lists.
  - A list of input files in which the gene appears.
  - A position vector, indicating the gene’s position in each list where it is present.
  - This integration enables full traceability and biological interpretability of the ranking.

- For **excluded genes**, a separate table is created containing:
  - The GeneID,
  - The number of appearances,
  - The names of the input files in which the gene was detected.
  - This simplified table is made available for inspection and optional download, but these genes are **not** used in any part of the ranking or enrichment process.

This two-table approach ensures a transparent analysis pipeline while maintaining performance and interpretability.


### Annotation Retrieval (Optional)

- A dedicated script located at `database_annotations/get_annotations.R` is used to generate local annotation files for **Gene Ontology (GO)**, **KEGG**, and **Reactome**. These files include `TERM2GENE` and `TERM2NAME` mappings required for enrichment analysis.
- Instead of relying on online-access functions like `enrichGO()` or `enrichKEGG()`, the system utilizes the more general `enricher()` function from the `clusterProfiler` package. This approach:
  - Loads annotations into memory at runtime.
  - Significantly improves performance.
  - Prevents errors caused by lack of internet connectivity or remote service timeouts.
- All annotation files are stored in the `/database_annotations/` directory and are **automatically loaded** by the app when enrichment is requested.

### Over-Representation Analysis (Optional)

- Based on user settings, a subset of the **top-ranked genes** (from 10 up to a maximum of 100) is selected from the consensus list.
- This selected gene set is then used to perform **over-representation analysis** against the locally loaded annotation databases.
- The enrichment analysis output includes:
  - **Term ID** (e.g., GO:0008150, R-HSA-123456),
  - **Description** of the biological term or pathway,
  - **Raw p-values**, and
  - **Associated gene sets** involved in the enrichment.

### Result Presentation

After the consensus analysis and optional enrichment, results are presented through multiple interactive and downloadable formats:

- **Interactive Results Table**:
  - Displays the final list of ranked genes.
  - Features include column visibility toggling, dynamic filtering, and downloadable formats (`.csv` and `.tsv`).

- **Excluded Gene Table**:
  - Displays genes filtered out due to low appearance frequency.
  - Includes number of appearances and list of files in which each gene was found.
  - Downloadable as `.TSV` only.
  - This table is optional and can be toggled on/off for inspection.

- **Upset Plot**:
  - Visualizes intersections between input lists (i.e., which genes are shared across how many lists).
  - Fully interactive with customization options.
  - Downloadable as `.PNG` and `.JPG`.

- **Heatmap**:
  - Shows the relative rank position of each gene across input lists.
  - Provides customization options for clustering, color schemes, and font sizes.
  - Downloadable as `.PNG`, `.JPG`, and interactive `.HTML`.

- **Enrichment Results Table** (optional):
  - Displays functional terms or biological pathways enriched among the selected genes.
  - Includes term ID, description, p-values, and matching genes.
  - Can be exported and explored alongside plots.

- **Plotting Options for Gene Ranking**:
  - Choose between **dot plot** or **bar plot** representations.
  - Customizable settings:
    - Number of top-ranked genes to show.
    - Color by p-value, rank, or appearance count.
    - Axis labels (e.g., Gene Symbol, Rank, Score).
    - Text size and color scale.
  - Download options include `.PNG`, `.JPG`, and interactive `.HTML`.

<br>

![](www/images/metarank_workflow.png){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: 12px;">
Figure 18: MetaRank workflow.
</p>


<br>

## Best Practices

To ensure accurate and meaningful results, users are advised to follow these best practices when using the app:

- **Choose the appropriate method for your data**:

  - Evaluate your dataset characteristics before selecting the meta-ranking algorithm.
  - If your data includes statistical values or fold changes, `RankProd` may be more suitable.
  - For simpler inputs or when only the order of genes matters, `RRA` offers a robust, non-parametric option.
  - Consider the presence of missing data, the number of lists, and the desired level of analytical stringency.

- **Ensure correct input formatting**:

  - RRA input files should contain one gene identifier per line (no headers).
  - RankProd files must contain two columns: gene ID and statistic (no header).
  - Input pasted directly into the app may include headers.
  - File format buttons and info modals provide validated examples — users are encouraged to consult them before uploading.
  - Use **official nomenclature** for gene identifiers (`SYMBOL`, `ENTREZID`, `ENSEMBL`) to ensure correct mapping.

- **Match gene ID type and organism**:

  - Confirm that the gene identifier type corresponds to the selected organism (e.g., `SYMBOL` with *Homo sapiens*).
  - Incorrect matching can lead to failed enrichment or unmapped terms.

- **Use provided example datasets**:

  - Preloaded example files demonstrate the expected input structure.
  - Running these examples helps verify that preprocessing, method selection, and plotting work as intended.

- **Tune the appearance filter carefully**:

  - A higher minimum appearance threshold increases the reliability of results but may exclude meaningful genes present in fewer lists.
  - A lower threshold includes more data but may increase background noise or false positives.
  - Adjust this parameter based on dataset size and study goals.

- **Set the enrichment gene count wisely**:

  - When running enrichment analysis, the number of top-ranked genes selected (e.g., 10 to 100) strongly impacts both result relevance and runtime.
  - Choose this number based on the biological context and size of the ranked list.

- **Inspect excluded genes and terms**:

  - If expected results are missing from the output, consult the excluded genes/terms tables.
  - This can reveal whether important items were filtered out due to low appearance or other criteria.

- **Customize visualizations for clarity**:

  - Modify axis labels, font sizes, color schemes, and number of top genes or terms to improve readability.
  - Interactive `.HTML` exports are ideal for in-depth exploration, while `.PNG` and `.JPG` are publication-ready.

- **Monitor computational performance**:

  - Large-scale analyses (e.g., >30 input lists or >1000 terms) may increase memory and runtime demands.
  - To optimize performance, reduce the number of input lists, narrow enrichment filters, or raise the minimum appearance threshold.

<br>

## Troubleshooting

Below is a list of common issues users may encounter during analysis, along with suggested solutions and interpretations.

| Issue | Possible Solution |
|-------|-------------------|
| **File Format Error** | This applies to both RankProd and RRA input files. Ensure that files are in `.txt` or `.csv` format without headers. For RankProd, each file must contain two columns: a gene identifier and a numeric statistic (e.g., logFC or fold change), separated by tabs or spaces. For RRA, the file must contain a single column of gene identifiers (one per line). Avoid special characters such as commas, semicolons, or `#`. Use info modals to preview accepted formats. |
| **Paste Format Error (RRA)** | When pasting data for RRA, each list must be separated by a line containing only `###`, and each gene must be on a separate line. Do not include headers or special characters. Avoid using tabular formatting or lists copied from spreadsheets. |
| **Paste Format Error (RankProd)** | For RankProd, pasted input should be structured with two columns per list, separated by tab or space: gene ID and statistic. Ensure there are no headers and that all values in the second column are numeric. Lines separating lists must include `###`. |
| **Invalid origin Field** | This applies to RankProd input only. Several possible issues can arise: (1) empty `origin` values (no data after `###`), (2) wrong separator (e.g., commas or semicolons instead of tab or space), (3) using non-numeric values (e.g., letters instead of logFC), (4) mismatch between number of provided `origin` values and number of lists, (5) missing replicates. The app will inform users of the specific problem encountered. |
| **Invalid Organism** | If the selected organism does not match the gene nomenclature used in the input, mapping may fail. Ensure that gene identifiers follow expected conventions: human (`TP53`, `BRCA1`), mouse (`Trp53`, `Brca1`), rat, etc. Use the dropdown to match the correct species code (e.g., *Hsa*, *Mmu*). |
| **Invalid Gene Identifiers** | Ensure consistency between the selected gene ID type and the actual format of your genes. `SYMBOL` = gene names like `TP53`; `ENTREZID` = numeric-only IDs; `ENSEMBL` = IDs like `ENSG00000...`. Mixing types can lead to errors or dropped genes. |
| **No Enrichment Results** | If enrichment tables return empty, it is likely due to a very strict filter or insufficient number of shared genes across lists. Try lowering the "Minimum Number of Lists" threshold and verify that the selected genes have known annotations. |
| **Appearance Threshold Too High** | If no genes remain after filtering, the "Minimum Number of Lists" threshold may be too restrictive. Reduce this value to allow inclusion of genes present in fewer lists. |
| **No Repeated Genes Across Lists** | If each list contains unique genes with no overlap, no consensus ranking or enrichment will be possible. Consider whether the gene lists are comparable and whether overlaps exist. |
| **Slow Performance** | Performance issues may arise when processing more than 30 lists or more than 1000 terms. To improve speed: reduce the number of lists, use fewer genes per list, or increase the minimum dataset threshold. Also consider simplifying visualization parameters. |

<p style="text-align: center; font-style: italic; color: #666; margin-top: 20px;">
Table 10: The list of possible errors that can be experienced and their possible cause
</p>


In addition to the summarized issues and suggestions presented in the table above, the following section visually illustrates the main error messages that may appear during the use of the application. Each message corresponds to a common input or processing issue, offering a brief explanation of what went wrong, how it affects the workflow, and what can be done to resolve it. In some cases, specific files or pieces of information that caused the problem are clearly indicated to help the user correct the input efficiently and continue without interruption.

::: panel-tabset
### File Format Error

![](www/images/metarank_error_1.png){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: 5px;">
Table 19: Pop-up window displayed when an error related to the file format is detected.
</p>

### Paste Format Error (RRA)

![](www/images/metarank_error_2.png){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: 5px;">
Table 20: Pop-up window displayed when an error related to the paste format (RRA) is detected.
</p>

### Paste Format Error (RankProd)

![](www/images/metarank_error_3.png){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: 5px;">
Table 21: Pop-up window displayed when an error related to the paste format (RP) is detected.
</p>

### Invalid origin Field

![](www/images/metarank_error_4.png){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: 5px;">
Table 22: Pop-up window displayed when an error related to the origin format is detected.
</p>

### Invalid Organism

![](www/images/metarank_error_5.png){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: 5px;">
Table 23: Pop-up window displayed when an error related to the organism selection is detected.
</p>


### Invalid Gene Identifiers

![](www/images/metarank_error_6.png){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: 5px;">
Table 24: Pop-up window displayed when an error related to the Gene ID selection is detected.
</p>

### No Enrichment Results

![](www/images/metarank_error_7.png){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: 5px;">
Table 25: Pop-up window displayed when no enrichment results were found.
</p>

### Appearance Threshold Too High

![](www/images/metarank_error_8.png){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: 5px;">
Table 26: Pop-up window that appears when no common genes are found between the lists established by the Appearance Threshold. 
</p>

### No Repeated Genes Across Lists

![](www/images/metarank_error_9.png){fig-align="center"}

<p style="text-align: center; font-style: italic; color: #666; margin-top: 5px;">
Table 27: Pop-up window that appears when no common genes are found. 
</p>

:::








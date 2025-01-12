# Project Title - Cancer Genomics


# Overview
This README provides an overview of the R code and analyses conducted for the Module 11 Assignment using the Maftools package. The assignment involves reading and summarizing Mutation Annotation Format (MAF) files, visualizing mutation data, and conducting various analyses related to genomic mutations in cancer datasets 

# Table of Contents
Reading and Summarizing MAF Files
Visualization
Analysis
    Somatic Interactions
    Detecting Cancer Driver Genes
    Adding and Summarizing Pfam Domains
    Survival Analysis
    Predicting Gene Sets Associated with Survival
    Comparing Two Cohorts (MAFs)
    Forest Plots
    Co-Oncoplots
    Lollipop Plot-2
    Clinical Enrichment Analysis
    Drug-Gene Interactions
    Oncogenic Signaling Pathways



# 6 Reading and summarizing maf files
6.1 Required input files
## 6.2 Reading MAF files.
In this section, the code focuses on reading MAF files using the maftools library. The MAF file (TCGA LAML) is loaded from a specified path, and optional clinical information, including survival data and histology, can be included. The read.maf function is then used to read the MAF file and associated clinical data, storing the result in the laml object for further analysis.

```{r}
library(maftools)

#path to TCGA LAML MAF file
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 

#clinical information containing survival information and histology. This is optional
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 

laml = read.maf(maf = laml.maf, clinicalData = laml.clin)
```

## 6.3 MAF object
In this section, various functions from the maftools library are used to analyze and summarize the loaded MAF file (laml). The code includes commands to display a basic summary of the MAF file, generate sample and gene summaries, show clinical data associated with samples, list all fields present in the MAF, and write a summary of the MAF to an output file with the basename 'laml'.

```{r}
#Typing laml shows basic summary of MAF file.
laml

#Shows sample summry.
getSampleSummary(laml)

#Shows gene summary.
getGeneSummary(laml)

#shows clinical data associated with samples
getClinicalData(laml)

#Shows all fields in MAF
getFields(laml)

#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = laml, basename = 'laml')
```

# 7 Visualization

## 7.1 Plotting MAF summary.
This code section uses the plotmafSummary function from the maftools library to generate a graphical summary of the MAF file (laml). The parameters used include removing outliers (rmOutlier = TRUE), adding the median mutation rate (addStat = 'median'), displaying a dashboard (dashboard = TRUE), and excluding the raw transition/transversion (Ti/Tv) ratio (titvRaw = FALSE). The resulting plot provides insights into the mutation landscape of the analyzed MAF file.
```{r}
# Plot a summary of the MAF file, including mutation statistics

plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

```

## 7.2 Oncoplots
## 7.2.1 Drawing oncoplots
In this section, the code is provided to draw oncoplots using the oncoplot function from the maftools library. The analysis focuses on the top ten mutated genes in the LAML (Acute Myeloid Leukemia) MAF file. The resulting oncoplot visually represents the mutation landscape of the specified genes in the analyzed dataset.

```{r}
#oncoplot for top ten mutated genes.
oncoplot(maf = laml, top = 10)
```


## 7.3 Transition and Transversions.
In this section, the code is provided to perform transition and transversion analysis using the titv function from the maftools library on the LAML (Acute Myeloid Leukemia) MAF file. The plotTiTv function is then used to visualize the transition and transversion summary. This analysis helps understand the ratio of transitions to transversions, providing insights into the mutational processes at play in the dataset.
```{r}
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)
```

## 7.5 Rainfall plots

In this section, the code focuses on creating rainfall plots using the maftools library. The BRCA (Breast Cancer) MAF file is loaded from the specified path, and the read.maf function is used to read the MAF file while suppressing verbose output. Subsequently, the rainfallPlot function is employed to generate a rainfall plot for the BRCA MAF file. The parameters include detecting change points (detectChangePoints = TRUE) and setting the point size to 0.4 (pointSize = 0.4). The resulting plot visually represents the distribution of mutations along the genome.

```{r}
# Specify the path to the BRCA (Breast Cancer) MAF file
brca <- system.file("extdata", "brca.maf.gz", package = "maftools")

# Read the BRCA MAF file, suppressing verbose output
brca = read.maf(maf = brca, verbose = FALSE)

# Generate a rainfall plot for the BRCA MAF file
rainfallPlot(maf = brca, detectChangePoints = TRUE, pointSize = 0.4)

```


## 7.6 Compare mutation load against TCGA cohorts
In this section, the code is provided to compare mutation loads against TCGA cohorts using the tcgaCompare function from the maftools library. The LAML (Acute Myeloid Leukemia) MAF file is used for this comparison, and specific parameters are set, including the cohort name ('Example-LAML'), enabling log scaling (logscale = TRUE), and specifying a capture size of 50 (capture_size = 50). The result is stored in the laml.mutload object for further analysis.
```{r}
laml.mutload = tcgaCompare(maf = laml, cohortName = 'Example-LAML', logscale = TRUE, capture_size = 50)

```


## 7.7 Plotting VAF
In this section, the code is provided to plot Variant Allele Frequency (VAF) using the plotVaf function from the maftools library. The LAML (Acute Myeloid Leukemia) MAF file is used, and the VAF column 'i_TumorVAF_WU' is specified for the plot. This plot visualizes the distribution of VAF values in the analyzed dataset.

```{r}

plotVaf(maf = laml, vafCol = 'i_TumorVAF_WU')

```

# 9 Analysis

## 9.1 Somatic Interactions
In this section, the code is provided to conduct exclusive/co-occurrence event analysis using the somaticInteractions function from the maftools library. The analysis is performed on the top 25 mutated genes in the LAML MAF file, and significance thresholds for p-values are set at 0.05 and 0.1. This function helps explore potential interactions between mutations in different genes in the specified dataset.
```{r}
#exclusive/co-occurance event analysis on top 10 mutated genes. 
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))
```


## 9.2 Detecting cancer driver genes based on positional clustering
In this section, the code is provided to detect cancer driver genes based on positional clustering using the oncodrive function from the maftools library. The analysis is performed on the LAML (Acute Myeloid Leukemia) MAF file, considering the 'Protein_Change' column for amino acid changes, requiring a minimum of 5 mutations (minMut = 5), and using the z-score method for p-value calculation (pvalMethod = 'zscore'). The results are stored in the laml.sig object, and the first few rows are displayed. Additionally, an oncodrive plot is generated with a false discovery rate cutoff of 0.1, utilizing specific plotting parameters for visualization.

```{r}
laml.sig = oncodrive(maf = laml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
head(laml.sig)
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)

```

## 9.3 Adding and summarizing pfam domains
In this section, the code is provided to add and summarize Pfam domains in the LAML (Acute Myeloid Leukemia) MAF file using the pfamDomains function from the maftools library. The analysis considers the 'Protein_Change' column for amino acid changes and focuses on the top 10 mutated genes (top = 10). The resulting protein and domain summaries are displayed, showing the first few columns for each summary for display convenience.

```{r}
laml.pfam = pfamDomains(maf = laml, AACol = 'Protein_Change', top = 10)

#Protein summary (Printing first 7 columns for display convenience)
laml.pfam$proteinSummary[,1:7, with = FALSE]

#Domain summary (Printing first 3 columns for display convenience)
laml.pfam$domainSummary[,1:3, with = FALSE]
```

## 9.4 Survival analysis
## 9.4.1 Mutation in any given genes
In this section, the code is provided to perform survival analysis based on the grouping of DNMT3A mutation status in the LAML (Acute Myeloid Leukemia) MAF file. The mafSurvival function from the maftools library is used, specifying the gene of interest ('DNMT3A'), the survival time variable ('days_to_last_followup'), the survival status variable ('Overall_Survival_Status'), and indicating that the dataset is from TCGA (isTCGA = TRUE). The analysis aims to explore the impact of DNMT3A mutations on overall survival in the specified dataset.

```{r}
#Survival analysis based on grouping of DNMT3A mutation status
mafSurvival(maf = laml, genes = 'DNMT3A', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = TRUE)
```

## 9.4.2 Predict genesets associated with survival
In this section, the code is provided to predict gene sets associated with survival in the LAML (Acute Myeloid Leukemia) MAF file using the survGroup and mafSurvGroup functions from the maftools library. The analysis uses the top 20 mutated genes to identify a set of genes (of size 2) to predict poor prognostic groups. The identified gene set is displayed, and a subsequent survival analysis is performed based on this gene set ('DNMT3A' and 'FLT3').

```{r}
#Using top 20 mutated genes to identify a set of genes (of size 2) to predict poor prognostic groups
prog_geneset = survGroup(maf = laml, top = 20, geneSetSize = 2, time = "days_to_last_followup", Status = "Overall_Survival_Status", verbose = FALSE)

print(prog_geneset)

mafSurvGroup(maf = laml, geneSet = c("DNMT3A", "FLT3"), time = "days_to_last_followup", Status = "Overall_Survival_Status")

```


## 9.5 Comparing two cohorts (MAFs)
In this section, the code is provided to compare two cohorts (Primary APL and Relapse APL) based on Mutation Annotation Format (MAF) files using the mafCompare function from the maftools library. The Primary APL and Relapse APL MAF files are loaded, and the comparison is performed considering only genes mutated in at least 5 samples in one of the cohorts to avoid bias. The results of the comparison are displayed using the print function.

```{r}
#Primary APL MAF
primary.apl = system.file("extdata", "APL_primary.maf.gz", package = "maftools")
primary.apl = read.maf(maf = primary.apl)
#Relapse APL MAF
relapse.apl = system.file("extdata", "APL_relapse.maf.gz", package = "maftools")
relapse.apl = read.maf(maf = relapse.apl)

#Considering only genes which are mutated in at-least in 5 samples in one of the cohort to avoid bias due to genes mutated in single sample.
pt.vs.rt <- mafCompare(m1 = primary.apl, m2 = relapse.apl, m1Name = 'Primary', m2Name = 'Relapse', minMut = 5)
print(pt.vs.rt)
```

## 9.5.1 Forest plots
In this section, the code is provided to generate a forest plot based on the results of the cohort comparison using the forestPlot function from the maftools library. The forest plot visualizes the comparison results, highlighting significant differences between the two cohorts with a p-value threshold of 0.1.

```{r}
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.1)

```


## 9.5.2 Co-onco plots
In this section, the code is provided to generate a co-oncoplot based on the comparison of Primary APL and Relapse APL MAFs using the coOncoplot function from the maftools library. The analysis focuses on the specified list of genes, and non-mutated samples are removed (removeNonMutated = TRUE). The resulting co-oncoplot visualizes the co-occurrence and mutual exclusivity of mutations in the specified genes between the two cohorts.

```{r}
genes = c("PML", "RARA", "RUNX1", "ARID1B", "FLT3")
coOncoplot(m1 = primary.apl, m2 = relapse.apl, m1Name = 'PrimaryAPL', m2Name = 'RelapseAPL', genes = genes, removeNonMutated = TRUE)
```

## 9.5.4 Lollipop plot-2
In this section, the code is provided to generate a lollipop plot comparing amino acid changes in the gene "PML" between Primary APL and Relapse APL using the lollipopPlot2 function from the maftools library. The lollipop plot visually represents the mutations in the specified gene for both cohorts, providing insights into changes in amino acid sequences between the two groups.
```{r}
lollipopPlot2(m1 = primary.apl, m2 = relapse.apl, gene = "PML", AACol1 = "amino_acid_change", AACol2 = "amino_acid_change", m1_name = "Primary", m2_name = "Relapse")

```

## 9.6 Clinical enrichment analysis
In this section, the code is provided to perform clinical enrichment analysis based on FAB classification in the LAML (Acute Myeloid Leukemia) MAF file using the clinicalEnrichment and plotEnrichmentResults functions from the maftools library. The analysis results are stored in the fab.ce object, and significant associations (p-value < 0.05) are displayed. The enrichment results are then visualized with a p-value cutoff 

```{r}
fab.ce = clinicalEnrichment(maf = laml, clinicalFeature = 'FAB_classification')

#Results are returned as a list. Significant associations p-value < 0.05
fab.ce$groupwise_comparision[p_value < 0.05]

plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05, geneFontSize = 0.5, annoFontSize = 0.6)

```

## 9.7 Drug-Gene Interactions
In this section, the code is provided to perform drug-gene interaction analysis using the drugInteractions function from the maftools library. The analysis is conducted on the LAML (Acute Myeloid Leukemia) MAF file, and the results are stored in the dgi object. Additionally, a specific analysis for the gene 'DNMT3A' is performed, and selected columns are displayed for the drug-gene interactions related to DNMT3A.

```{r}
dgi = drugInteractions(maf = laml, fontSize = 0.75)
dnmt3a.dgi = drugInteractions(genes = "DNMT3A", drugs = TRUE)
#Printing selected columns.
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]

```

## 9.8 Oncogenic Signaling Pathways
In this section, the code is provided to analyze oncogenic signaling pathways in the LAML (Acute Myeloid Leukemia) MAF file using the OncogenicPathways function from the maftools library. The results provide insights into the activation status of various signaling pathways. Additionally, a specific oncogenic signaling pathway ("RTK-RAS") is plotted using the PlotOncogenicPathways function for visual representation.

```{r}
OncogenicPathways(maf = laml)
PlotOncogenicPathways(maf = laml, pathways = "RTK-RAS")

```

# References
Mayakonda, Anand, De-Chen Lin, Yassen Assenov, Christoph Plass, and H. Phillip Koeffler. 2018. “Maftools: Efficient and Comprehensive Analysis of Somatic Variants in Cancer.” *Genome Research* 28 (11): 1747–56..

Network, Cancer Genome Atlas Research, Timothy J. Ley, Christopher Miller, Li Ding, Benjamin J. Raphael, Andrew J. Mungall, A. Gordon Robertson, et al. 2013. “Genomic and Epigenomic Landscapes of Adult de Novo Acute Myeloid Leukemia.” *The New England Journal of Medicine* 368 (22): 2059–74..

**Acknowledgment**
Developed at Northeastern University for BINF6308 coursework.
Author: Renuka Athinarayanan



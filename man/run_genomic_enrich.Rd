% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downstream_analysis.R
\name{run_genomic_enrich}
\alias{run_genomic_enrich}
\title{Running genomic enrichment analysis}
\usage{
run_genomic_enrich(
  model_results,
  species = "hsa",
  pvalueCutoff = 0.1,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.1,
  ont = "all",
  BPPARAM = BiocParallel::SerialParam(),
  extracted_data = NULL,
  ...
)
}
\arguments{
\item{model_results}{Model integration results, typically a list containing
different types of genomic results}

\item{species}{Species to select for the enrichment analysis. Default is
'hsa' (Homo sapiens).}

\item{pvalueCutoff}{P-value cutoff for significant enrichment. Default is
0.1.}

\item{pAdjustMethod}{Method for adjusting p-values. Default is 'BH'
(Benjamini & Hochberg).}

\item{qvalueCutoff}{Q-value cutoff for significant enrichment. Default is
0.1.}

\item{ont}{Ontology to use for the enrichment analysis. Default is 'all'.}

\item{BPPARAM}{A BiocParallelParam object specifying parallelization
options. Default is BiocParallel::SerialParam().}

\item{extracted_data}{Pre-extracted data for enrichment analysis. If NULL,
function will extract relevant data from model_results.}

\item{...}{Additional arguments to be passed to the internal enrichment
function.}
}
\value{
A list containing enrichment results. If CNV and methylation data
are available, it returns a nested list with results for each data type.
}
\description{
Running genomic enrichment analysis
}
\examples{
# Example usage:
library(MultiAssayExperiment)
data(mmultiassay_ov)
tmp <- lapply(experiments(mmultiassay_ov), function(x) x[1:200,])
mmultiassay_ov <- MultiAssayExperiment(experiments = tmp)
#multiomics_integration <- run_multiomics(mmultiassay_ov)
#gen_enr <- run_genomic_enrich(multiomics_integration, qvalueCutoff = 1,
#pvalueCutoff = 0.05, pAdjustMethod = 'none')
}

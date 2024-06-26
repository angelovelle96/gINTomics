% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiomics.R
\name{run_genomic_integration}
\alias{run_genomic_integration}
\title{Integration of expression, Copy Number Variations and methylation data}
\usage{
run_genomic_integration(
  expression,
  cnv_data,
  methylation,
  sequencing_data = TRUE,
  normalize = TRUE,
  norm_method = "TMM",
  interactions = NULL,
  class = NULL,
  scale = TRUE,
  run_deg = TRUE,
  BPPARAM = SerialParam(),
  ...
)
}
\arguments{
\item{expression}{Matrix or data.frame containing the expression values
for each model. Rows represent samples, while each column represents
the different response variables of the models.}

\item{cnv_data}{Matrix or data.frame containing the Copy Number variation
status for the models. Rows represent samples, while columns represent
the different covariates. If \strong{interactions} are not provided, they will be
automatically generated and for each gene contained in \strong{expression}
the model will look for the same gene in \strong{cnv_data}}

\item{methylation}{Matrix or data.frame containing the methylation
values for the models. Rows represent samples, while columns represent
the different covariates. If \strong{interactions} are not provided, they will be
automatically generated and for each gene contained in \strong{expression}
the model will look for the same gene in \strong{methylation}}

\item{sequencing_data}{logical. Are expression data obtained from RNA
sequencing ? Default is set to TRUE}

\item{normalize}{logical.Should expression data be
normalized ? Default is set to TRUE}

\item{norm_method}{Normalization method to be used for
expression data. One of "TMM" (default), "TMMwsp", "RLE", "upperquartile",
"none".}

\item{interactions}{A list of character vectors containing the interactions
between response variable and covariates. The names of the list should
match the response variables while the character contained in each element
of the list should match the covariates. If NULL (default), the interactions
will be automatically defined according to response variable's colnames.}

\item{class}{Character vector specifying the classes for differential
expression analysis.}

\item{scale}{Logical. Should the data be scaled? Default is set to TRUE.}

\item{run_deg}{Logical. Should differential expression analysis be performed?
Default is set to TRUE.}

\item{BPPARAM}{A BiocParallelParam object specifying the parallel backend to
be used.}

\item{...}{Additional arguments to be passed to internal functions.}
}
\value{
A list or a \linkS4class{MultiClass} object if \strong{class} is provided
containing the results of the Genomic integration
}
\description{
This function will perform an integration of expression data and Copy Number
Variations data
}
\examples{
# Example usage_multi:
library(MultiAssayExperiment)
data("mmultiassay_ov")
tmp <- lapply(experiments(mmultiassay_ov), function(x) x[1:20,])
mmultiassay_ov <- MultiAssayExperiment(experiments = tmp)
meth_matrix <- t(as.matrix(assay(mmultiassay_ov[["methylation"]])))
gene_exp_matrix <- t(as.matrix(assay(mmultiassay_ov[["gene_exp"]])))
gene_cnv_matrix <- t(as.matrix(assay(mmultiassay_ov[["cnv_data"]])))
genomic_integration <- run_genomic_integration(
    expression = gene_exp_matrix,
    cnv_data = gene_cnv_matrix, methylation = meth_matrix
)
}

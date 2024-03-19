<img width="250" height="250" src="vignettes/figures/logo_gINTomics.png" class="center" />

  gINTomics is an R package for Multi-Omics data integration and visualization. gINTomics is designed to detect the association between the expression of a target and of its regulators, taking into account also their genomics modifications such as Copy Number Variations (CNV) and methylation. For RNA sequencing data, the counts will be fitted using a negative binomial model, while in the case of microarray or other types of data, a linear model will be applied. In some cases the number of regulators for a given target could be very high, in order to handle this eventuality, we provide a random forest selection that will automatically keep only the most important regulators. All the models are gene-specific, so each gene/miRNA will have its own model with its covariates. The package will automatically download information about TF-target couples (OmnipathR), miRNA-target couples (OmnipathR) and TF-miRNA couples (TransmiR). The couples will be used to define the covariates used in the integration models.


# Installation
To install this package:

``` r
devtools::install_github("angelovelle96/gINTomics")
```

# How to use gINTomics
## Structure of the package
gINTomics is designed to perform both two-Omics and Multi Omics integrations. There are four different functions for the two-omics integration. The `run_cnv_integration` function can be used to integrate Gene or miRNA Expression data with Copy Number Variation data. The `run_met_integration` is designed for the integration between Gene Expression data with methylation data. When both gene CNV and methylation data are available, the `run_genomic_integration` function can be used to integrate them together with gene expression. Finally, the `run_tf_integration` function can be used for the integration between the Expression of a target gene or miRNA and every kind of regulator, such as Transcription Factors or miRNAs expression.
The package also gives the possibility to perform Multi Omics integration. The `run_multiomics` function takes as input a MultiAssayExperiment, that can be generated with the help of the `create_multiassay` function, and will perform all the possible integration available for the omics provided by the user.
All the integration functions exploit different statistical models depending on the input data type, the available models are negative binomial and linear regression. When the number of covariates is too high, a random forest model will select only the most importants which will be included in the integration models.
In order to make the results more interpretable, gINTomics provides a comprehensive and interactive *shiny app* for the graphical representation of the results, the shiny app can be called with the `run_shiny` function, which takes as input the output of the `run_multiomics` function. Moreover additional functions are available to build single plots to visualize the results of the integration outside of the shiny app.
In the following sections, we use an example MultiAssayExperiment of ovarian cancer to show how to use `gINTomics` with the different available integrations. The object contains all the available input data types: Gene expression data, miRNA expression data, gene methylation data, gene Copy Number Variations and miRNA Copy Number Variations.

## Generate a MultiAssayExperiment
The package contains a pre built MultiAssayExperiment, anyway in this section we will divide it in signle omics matrices and see how to build a proper MultiAssayExperiment from zero.

``` r
# loading packages
library(gINTomics)
library(MultiAssayExperiment)
data("ov_test_tcga_omics")
``` 

Here we will extract the single omics from the MultiAssayExperiment
``` r
gene_exp_matrix <- as.matrix(assay(mmultiassay_ov[["gene_exp"]]))
miRNA_exp_matrix <- as.matrix(assay(mmultiassay_ov[["miRNA_exp"]]))
meth_matrix <- as.matrix(assay(mmultiassay_ov[["methylation"]]))
gene_cnv_matrix <- as.matrix(assay(mmultiassay_ov[["cnv_data"]]))
miRNA_cnv_matrix <- as.matrix(assay(mmultiassay_ov[["miRNA_cnv_data"]]))
``` 

Now let's build a proper MultiAssayExperiment starting from single matrices

``` r
new_multiassay <- create_multiassay(methylation = meth_matrix, 
                                    gene_exp = gene_exp_matrix,
                                    cnv_data = gene_cnv_matrix,
                                    miRNA_exp = miRNA_exp_matrix,
                                    miRNA_cnv_data = miRNA_cnv_matrix)
``` 


## Run Genomic integration

In this section we will see how to perform a gene_expression~CNV+met integration. The input data should be provided as matrices or data.frames. Rows represent samples, while each column represents the different response variables/covariates of the models. By default all the integration functions will assume that expression data are obtained through sequencing, you should set the `sequencing_data` argument to FALSE if they are not.  Expression data could be both normalized or not normalized, the function is able to normalize data by setting the `normalize` argument to TRUE (default). If you want you can specify custom interaction through the `interactions` argument otherwise  the function will automatically look for the genes of `expression` in `cnv_data`. 
NOTE: if you customize genomics interactions, linking more covariates to a single response, the output may
be not more compatible with the shiny visualization.

``` r
gene_genomic_integration <- run_genomic_integration(expression = t(gene_exp_matrix),
                                                    cnv_data = t(gene_cnv_matrix),
                                                    methylation = t(meth_matrix))
summary(gene_genomic_integration)
```


## Run CNV integration

In this section we will see how to perform an expression~CNV integration. The input data (for both expression and CNV) should be provided as matrices or data.frames. Rows represent samples, while each column represents the different response variables/covariates of the models. By default all the integration functions will assume that expression data are obtained through sequencing, you should set the `sequencing_data` argument to FALSE if they are not.  Expression data could be both normalized or not normalized, the function is able to normalize data by setting the `normalize` argument to TRUE (default). If you want you can specify custom interaction through the `interactions` argument otherwise  the function will automatically look for the genes of `expression` in `cnv_data`.

``` r
gene_cnv_integration <- run_cnv_integration(expression = t(gene_exp_matrix),
                                            cnv_data = t(gene_cnv_matrix))
summary(gene_cnv_integration)
```

## Run methylation integration

In this section we will see how to perform an expression~methylation integration. The input data should be the same of `run_cnv_integration`, but the covariates matrix will contain methylation data instead of Copy Number Variations data. Expression data could be both normalized or not normalized, the function is able to normalize data by setting the `normalize` argument to TRUE (default). If you want you can specify custom interaction through the `interactions` argument otherwise  the function will automatically look for the genes of `expression` in `methylation`.

``` r
gene_met_integration <- run_met_integration(expression = t(gene_exp_matrix),
                                            methylation = t(meth_matrix))
summary(gene_met_integration)
```

## Run TF-target integration

In this section we will see how to perform an expression~TF integration. In this case you can use as input data a single gene expression matrix if both your TF and targets are genes and are contained in the gene expression matrix. If you want the package to automatically download the interactions between TFs and targets you have to set the argument `type` to "tf". Otherwise you can also specify your custom interactions providing them with the `interactions` argument. You can handle data normalization as in the previous functions through the `normalize` argument. This function is designed to integrate TF expression data but it can handle every type of numerical data representing a gene expression regulator. So you can pass the regulators matrix to the `tf_expression` argument and specify your custom `interactions`. If expression data are not obtained through sequencing, remember to set `sequencing data` to FALSE.

``` r
tf_target_integration <- run_tf_integration(expression = t(gene_exp_matrix),
                                            type = "tf")
summary(tf_target_integration)
```


## Run miRNA-target integration

In this section we will see how to perform an expression~miRNA integration. Gene expression data will be provided through the `expression` argument while miRNA expression data through the `tf_expression` argument. Input matrices should follow the same rules of the previous functions. If you want the package to automatically download the interactions between miRNA and targets you have to set the argument `type` to "miRNA_target". Otherwise you can also specify your custom interactions providing them with the `interactions` argument. You can handle data normalization as in the previous functions through the `normalize` argument.

``` r
miRNA_target_integration <- run_tf_integration(expression = t(gene_exp_matrix),
                                               tf_expression = t(miRNA_exp_matrix),
                                               type = "miRNA_target")
summary(miRNA_target_integration)
```


## Run TF-miRNA integration

In this section we will see how to perform an miRNA~TF integration. miRNA expression data will be provided through the `expression` argument while gene expression data through the `tf_expression` argument. Input matrices should follow the same rules of the previous functions. If you want the package to automatically download the interactions between TF and miRNA you have to set the argument `type` to "tf_miRNA". Otherwise you can also specify your custom interactions providing them with the `interactions` argument. You can handle data normalization as in the previous functions through the `normalize` argument.

``` r
tf_miRNA_integration <- run_tf_integration(expression = t(miRNA_exp_matrix),
                                               tf_expression = t(gene_exp_matrix),
                                               type = "tf_miRNA")
summary(tf_miRNA_integration)
```


## Run complete Multi-Omics integration

Finally we will see how to perform a complete integration using all the available omics. In order to run this function you need to generate a MultiAssayExperiment as showed at the beginning of this vignette. The function will automatically use all the omics contained in the MultiAssayExperiment to perform all the possible integrations showed before. Moreover, if genomic data are available (CNV and/or Methylation), the first step will be the genomic integration and all the following integrations that contain gene expression data as response variable will use instead the residuals of the genomic model in order to filter out the effect of CNV and/or methylation. This framework is used also for miRNA, but in this case only CNV data are supported.

``` r
multiomics_integration <- run_multiomics(data = new_multiassay)
summary(multiomics_integration)
```


# Visualization

Now let's see how you can visualize the results of your integration models.

## Shiny app
gINTomics provides an interactive environment, based on Shiny, that allows the user to easily visualize output results from integration models and to save them for downstream tasks and reports creation.
Once multiomics integration has been performed users can provide the results to the `run_shiny` function.
NOTE: Only the output of the `run_multiomics` function is compatible with `run_shiny`
``` r
run_shiny(multiomics_integration)
```



## plots
### network plot
Network plot shows the significant interactions between transcriptional regulators (TFs and miRNAs, if present) and their targets genes/miRNA. Nodes and edges are selected ordering them by the most high coefficient values (absolute value) and by default the top 200 interactions are showed. You can change the number of interactions with the `num_interactions` argument.
Here you can see the code to generate a network from a multiomics integration
``` r
data_table <- extract_model_res(multiomics_integration)
data_table <- data_table[data_table$cov!="(Intercept)",]
plot_network(data_table, num_interactions = 200)

```

### Venn Diagram
The Venn Diagram is designed for the genomic integration. It can help to identify genes which are significantly regulated by both CNV and methylation.
Here you can see the code to generate a Venn Diagram from a multiomics integration
``` r
plot_venn(data_table)

```


### Volcano plot
The Volcano Plot shows the distribution of integration coefficients for every integration type associated with a genomic class (cnv, met, cnv-mirna). For each integration coefficient, on the y axis you have the -log10 of Pvalue/FDR and on the x axis the value of the coefficient.
Here you can see the code to generate a Volcano plot for CNV and one for methylation from a multiomics integration
``` r
plot_volcano(data_table, omics = "gene_genomic_res", cnv_met = "cnv")
plot_volcano(data_table, omics = "gene_genomic_res", cnv_met = "met")

```


### Ridgeline plot
The ridgeline plot is designed to compare different distributions, it has been integrated in the package with the aim to compare the distribution of significant and non significant coefficients returned by our integration models.
For each distribution, on the y axis you have the frequencies and on the x axis the values of the coefficients.
Here you can see the code to generate a Ridgeline plot for CNV and one for methylation from a multiomics integration
``` r
plot_ridge(data_table, omics = "gene_genomic_res", cnv_met = "cnv")
plot_ridge(data_table, omics = "gene_genomic_res", cnv_met = "met")

```


### Chromosome distribution plot
This barplot highlights the distribution of significant and non significant covariates among chromosomes. This kind of visualization will identify chromosomes in which the type of regulation under analysis is particularly active.
Here you can see the code to generate a chromosome distribution plot specif for the genomic integration, starting from the results of the multiomics integration.
``` r
plot_chr_distribution(data_table = data_table,
                     omics = "gene_genomic_res")

```

### TF distribution plot
This barplot highlights the number of significant tragets for each TFs.
Here you can see the code to generate a TF distribution plot starting from the results of the multiomics integration.
``` r
plot_tf_distribution(data_table = data_table)

```

### Enrichment plot
The enrichment plot shows the enrichment results obtained with enrichGO and enrichKEGG (clusterProfiler). The genomic enrichment is performed providing the list of genes significantly regulated by methylation or CNV, while the transcriptional one with the list of genes significantly regulated by each transcription factor (we run an enrichment for each TF that significantly regulates at least 12 targets)
Here you can see the code to run a genomic enrichment starting from the results of the multiomics integration and to visualize the top results with a DotPlot
``` r
gen_enr <- run_genomic_enrich(multiomics_integration, qvalueCutoff = 1, pvalueCutoff = 0.05, pAdjustMethod = "none")
dot_plotly(gen_enr$cnv[[1]]$go)

```


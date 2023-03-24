run_multiomics <- function(data=NULL,
                           use_cnv_res=T,
                           interactions_met=NULL,
                           interactions_miRNA=NULL,
                           interactions_tf=NULL,
                           ...){


#si lavora con multiassay experiment, va creato con la funzione create_multia..
#prima integro il copy number e l'espressione e do la possibilità
#di usare i residui per le analisi successive
#prima scelta, se l'espressione è un dato di conte usiamo edgeR altrimenti
#usiamo modello lineare

    if(class(data)!="MultiAssayExperiment"){
      stop(str_wrap("You need to provide a MultiassayExperiment as input
                    data. The names of the experiments included in the
                    MultiAssay should follow the rule provided in the
                    help section. You can use the function create_multiassay
                    to create proper input data"))
    }

    if(!is.null(data@ExperimentList$cnv_data) &
       !is.null(data@ExperimentList$gene_exp)){
        gene_cnv_res <- run_cnv_integration(
          expression = t(assay(data, i = "gene_exp")),
          cnv_data = t(assay(data, i = "cnv_data")),
          threads=2,
          sequencing_data = F)

    }

    if(!is.null(data@ExperimentList$miRNA_cnv_data) &
       !is.null(data@ExperimentList$miRNA_exp)){
      mirna_cnv_res <- run_cnv_integration(
        expression = t(assay(data, i = "miRNA_exp")),
        cnv_data = t(assay(data, i = "miRNA_cnv_data")),
        threads=8,
        sequencing_data = T)

    }

  if(!is.null(data@ExperimentList$methylation) &
     !is.null(data@ExperimentList$miRNA_exp)){
    met_res <- run_met_integration(
      expression = t(assay(data, i = "gene_exp")),
      methylation = t(assay(data, i = "methylation")),
      threads=16,
      sequencing_data = T)

  }







#proseguo con le altre integrazioni


}

run_multiomics <- function(cnv_data=NULL,
                           gene_expr=NULL,
                           methylation=NULL,
                           miRNA_expr=NULL,
                           interactions_met=NULL,
                           interactions_miRNA=NULL,
                           interactions_tf=NULL){



#prima integro il copy number e l'espressione e do la possibilità
#di usare i residui per le analisi successive
#prima scelta, se l'espressione è un dato di conte usiamo edgeR altrimenti
#usiamo modello lineare

    if(is.null(cnv_data)){
        run_edgeR_integration(response_var = gene_expr,
                                covariates = cnv_data)


    }




#proseguo con le altre integrazioni


}

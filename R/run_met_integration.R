run_met_integration <- function( expression,
                                methylation,
                                sequencing_data=T,
                                interactions=NULL,
                                ...){


  if(is.null(interactions)){
    interactions <- as.list(intersect(colnames(methylation),
                                      colnames(expression)))
    names(interactions) <- unlist(interactions)
  }


  if(sequencing_data==T){
    met_res <- run_edgeR_integration(response_var = expression,
                                    covariates = methylation,
                                    interactions = interactions,
                                    ...)
  }else{

    met_res <- run_lm_integration(response_var = expression,
                                     covariates = methylation,
                                     interactions = interactions,
                                     ...)

  }



  return(met_res)


}

run_met_integration <- function( expression,
                                methylation,
                                sequencing_data=T,
                                ...){


  if(sequencing_data==T){
    met_res <- run_edgeR_integration(response_var = expression,
                                    covariates = methylation,
                                    cnv_mode = T,
                                    ...)
  }



  return(met_res)


}

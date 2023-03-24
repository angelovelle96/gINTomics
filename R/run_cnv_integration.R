run_cnv_integration <- function(expression,
                                cnv_data,
                                sequencing_data=T,
                                ...){


    if(sequencing_data==T){
        cnv_res <- run_edgeR_integration(response_var = expression,
                                        covariates = cnv_data,
                                        ...)
    }else{

        cnv_res <- run_lm_integration(response_var = expression,
                                     covariates = cnv_data,
                                     ...)
    }



    return(cnv_res)
}

run_tf_integration <- function( expression,
                                tf_expression,
                                interactions,
                                sequencing_data=T,
                                ...){


    if(sequencing_data==T){
      tf_res <- run_edgeR_integration(response_var = expression,
                                       covariates = tf_expression,
                                       interactions = interactions,
                                       ...)
    }else{

      tf_res <- run_lm_integration(response_var = expression,
                                    covariates = tf_expression,
                                    interactions = interactions,
                                    ...)
    }



    return(tf_res)


}

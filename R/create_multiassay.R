create_multiassay <- function(methylation=NULL,
                              cnv_data=NULL,
                              gene_exp=NULL,
                              miRNA_exp=NULL,
                              miRNA_cnv_data=NULL,
                              ...){


  if((is.null(methylation)+is.null(cnv_data)+
      is.null(gene_exp)+is.null(miRNA_exp)+is.null(miRNA_cnv_data))>=4){
    stop(str_wrap("You need to provide at least two assays or
                  SummarizedExperiment"))
  }




  eexperiments <- Filter(Negate(is.null),
                    list(gene_exp=gene_exp,
                         methylation=methylation,
                         cnv_data=cnv_data,
                         miRNA_exp=miRNA_exp,
                         miRNA_cnv_data=miRNA_cnv_data))

  mmultiassay <- MultiAssayExperiment(experiments = eexperiments)
  mmultiassay <- intersectColumns(mmultiassay)



}

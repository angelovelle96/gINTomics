formatR::tidy_dir("~/Documenti/Rpackages/integrazione/R/allgene_edgeR_model.R",
                  indent = 4, recursive = T, width.cutoff=80)
styler::style_file("~/Documenti/Rpackages/integrazione/R/run_edgeR_integration.R",
                   transformers = tidyverse_style(indent_by = 4))



load("~/Documenti/Rpackages/integrazione/data/miRBase_conversion.rda")
load("~/Documenti/Rpackages/integrazione/data/miRNA_TF.rda")

save(mirna_hsa, tf_mirna,
     file = "~/Documenti/Rpackages/integrazione/R/sysdata.rda")



load("~/Documenti/Ovarian/rdata/parametri_integrazione_mirna_tf_singleseq_isoall_20220329.Rdata")
load("~/Documenti/Ovarian/rdata/parametri_integrazione_mirna_all_samples_ace20220311.Rdata")


interactions <- expanded_tf_mirna_couples
interactions <- lapply(interactions, function(x) x[,1])
response_var <-  run_edgeR_test2_input$mirna_exp_model
covariates <-  run_edgeR_test2_input$tf_expression_model

a4 <- run_lm_integration( response_var = response_var,
                    covariates = covariates,
                    interactions = interactions,
                    threads = 10,
                    step = T)



counts <- run_edgeR_test1_input$counts
data <- run_edgeR_test1_input$data
tmp <- data[, 1:(ncol(data)-3)]
tmp <- apply(tmp, 2, function(x) ifelse(x>2, "amp", ifelse(x==2, "norm", "del")))
tmp <- as.data.frame(tmp)
tmp[sapply(tmp, is.character)] <-  lapply(tmp[sapply(tmp, is.character)], as.factor)
data[, 1:(ncol(data)-3)] <- tmp
tmp <- sapply(data[,1:ncol(data)], function(x) length(levels(x)))
data <- data[, tmp!=1]
counts <- counts[intersect(rownames(counts), colnames(data)),]

expectedres <- run_edgeR_test1_output
design_all <- model.matrix(~W_1+W_2+W_3, data = data)
y_all_offset <- matrix(1, ncol(counts), nrow(counts))
y_gene_offset <- lapply(1:nrow(counts), function(x) {
  matrix(1, 1, ncol(counts))
})
prova <- run_edgeR_integration( response_var = t(counts),
                                covariates = data,
                                steady_covariates = c("W_1", "W_2", "W_3"),
                                design_mat_allgene = design_all,
                                offset_allgene = t(y_all_offset),
                                offset_singlegene = y_gene_offset,
                                threads = 16,
                                norm_method = "TMM",
                                cnv_mode = T)

a <- sapply(colnames(data), function(x) length(levels(data[,x])))
b <- data[,names(a)[a==2]]


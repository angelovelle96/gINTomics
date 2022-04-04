formatR::tidy_dir("~/Documenti/Rpackages/integrazione/R/allgene_edgeR_model.R",
                  indent = 4, recursive = T, width.cutoff=80)
styler::style_file("~/Documenti/Rpackages/integrazione/R/run_edgeR_integration.R",
                   transformers = tidyverse_style(indent_by = 4))

load("~/Documenti/Rpackages/integrazione/data/mirna_exp_model.rda")
load("~/Documenti/Rpackages/integrazione/data/tf_expression_model.rda")
load("~/Documenti/Rpackages/integrazione/data/tf_mirna_couples.Rdata")
load("~/Documenti/Ovarian/rdata/parametri_integrazione_mirna_tf_singleseq_isoall_20220329.Rdata")
load("~/Documenti/Ovarian/rdata/parametri_integrazione_mirna_all_samples_ace20220311")



# save(run_edgeR_test1_input,
# run_edgeR_test1_output,
# run_edgeR_test2_input,
# run_edgeR_test2_output,
# file =
# '~/Documenti/Rpackages/integrazione/R/sysdata.rda')


a <- generate_design(response_var = run_edgeR_test1_input$counts,
                     covariates = run_edgeR_test1_input$data,
                     steady_covariates = c("W_1", "W_2", "W_3"),
                     cnv_mode = T)


design_mat_allgene <- model.matrix(~1, data = covariates)
y_all <- DGEList(counts = response_var)
y_all <- calcNormFactors(y_all, method = norm_method)
if (!is.null(offset_allgene)) y_all$offset <- offset_allgene
y_all <- estimateGLMCommonDisp(y_all, design_mat_allgene)
y_all <- estimateGLMTagwiseDisp(y_all, design_mat_allgene)
fit <- glmFit(y_all, design_mat_allgene)



design_mat_allgene2 <- model.matrix(~AKT3, data = covariates)
y_all2 <- DGEList(counts = response_var)
y_all2 <- calcNormFactors(y_all2, method = norm_method)
y_all2 <- estimateGLMCommonDisp(y_all2, design_mat_allgene2)
y_all2 <- estimateGLMTagwiseDisp(y_all2, design_mat_allgene2)
fit2 <- glmFit(y_all2, design_mat_allgene2)












mirna_exp_model <- run_edgeR_test2_input$mirna_exp_model
tf_expression_model <- run_edgeR_test2_input$tf_expression_model
interactions <- run_edgeR_test2_input$expanded_tf_mirna_couples
interactions <- lapply(interactions, function(x) x[,1])
expectedres <- run_edgeR_test2_output
y_all_offset <- matrix(1, ncol(mirna_exp_model), nrow(mirna_exp_model))

myres <- run_edgeR_integration( response_var = mirna_exp_model,
                                covariates = tf_expression_model,
                                interactions = interactions,
                                threads = 16,
                                norm_method = "TMM")
# myres2 <- run_edgeR_integration(response_var = mirna_exp_model,
#                                 design_mat_allgene = design_all,
#                                 interactions = interactions,
#                                 covariates = tf_expression_model,
#                                 threads = 16,
#                                 norm_method = "TMM")

colnames(myres$coef_data) <- paste0(colnames(myres$coef_data), "_coef")
colnames(myres$pval_data) <- paste0(colnames(myres$pval_data), "_pvalue")
tmp <- cbind(myres$coef_data, myres$pval_data)
tmp <- tmp[rownames(expectedres), colnames(expectedres)]
# colnames(myres2$coef_data) <- paste0(colnames(myres2$coef_data), "_coef")
# colnames(myres2$pval_data) <- paste0(colnames(myres2$pval_data), "_pvalue")
# tmp2 <- cbind(myres2$coef_data, myres2$pval_data)
# tmp2 <- tmp2[rownames(expectedres), colnames(expectedres)]
identical(tmp, expectedres)
# expect_identical(tmp2, expectedres)
# expect_identical(myres, myres2)
# })

generate_design(response_var = run_edgeR_test1_input$counts,
                            covariates = run_edgeR_test1_input$data,
                            interactions = NULL,
                            steady_covariates = NULL,
                            cnv_mode = T)

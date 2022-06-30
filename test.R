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


interactions <- expanded_tf_mirna_couples
interactions <- lapply(interactions, function(x) x[,1])
response_var <-  run_edgeR_test2_input$mirna_exp_model
covariates <-  run_edgeR_test2_input$tf_expression_model

a <- run_lm_integration( response_var = response_var,
                    covariates = covariates,
                    interactions = interactions,
                    threads = 10,
                    step = T)





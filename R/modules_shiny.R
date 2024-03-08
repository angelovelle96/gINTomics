.server_histo_gen <- function(data_table,
                             input,
                             output,
                             session){

reactive_histo <- .prepare_reactive_histo(data_table,
                                          input = input,
                                          output = output,
                                          deg = FALSE)
output$plotly <- .render_histo(reactive_histo)
reactive_histo_table <- .prepare_reactive_histo_table(data_table,
                                                      input=input,
                                                      output=output,
                                                      deg=FALSE)
output$table <- .render_histo_table(reactive_histo_table)
output$download_csv <- download_csv(deg=FALSE,
                                    plotType="histo",
                                    input=input,
                                    output=output,
                                    data_table=data_table)

}


.server_histo_trans <- function(data_table,
                              input,
                              output,
                              session){

  reactive_histo_transcript <- .prepare_reactive_histo(data_table,
                                                       input=input,
                                                       output=output,
                                                       deg=FALSE)
  output$plotly <- .render_histo(reactive_histo_transcript)
  reactive_histo_table_transcript <- .prepare_reactive_histo_table(data_table,
                                                                   input=input,
                                                                   output=output,
                                                                   deg=FALSE)
  output$table <- .render_histo_table(reactive_histo_table_transcript)
  output$download_csv <- download_csv(deg=FALSE,
                                      plotType="histo",
                                      input=input,
                                      output=output,
                                      data_table=data_table)

  reactive_histo_tf_transcript <- .prepare_reactive_histo_tf(data_table,
                                                             input = input,
                                                             output = output,
                                                             deg = FALSE)
  output$plotly_tf <- .render_histo_TF(reactive_histo_tf_transcript)
  reactive_tf_table <- .prepare_reactive_histoTF_table(data_table,
                                                       input=input,
                                                       output=output,
                                                       deg=FALSE)
  output$download_csv_tf <- download_csv(deg=FALSE,
                                         tf=TRUE,
                                         plotType="histo",
                                         input=input,
                                         output=output,
                                         data_table=data_table)
  output$table_tf <- .render_histo_tf_table(reactive_tf_table)

}


.server_histo_deg <- function(data_table,
                                input,
                                output,
                                session){
  reactive_histo_deg <- .prepare_reactive_histo(data_table,
                                                input=input,
                                                output=output,
                                                deg=TRUE)
  output$plotly <- .render_histo(reactive_histo_deg)
  reactive_histo_table_deg <- .prepare_reactive_histo_table(data_table,
                                                            input=input,
                                                            output=output,
                                                            deg=TRUE)
  output$table <- .render_histo_table(reactive_histo_table_deg)
  output$download_csv <- download_csv(deg=TRUE,
                                      plotType="histo",
                                      input=input,
                                      output=output,
                                      data_table=data_table)
  reactive_histo_tf_deg <- .prepare_reactive_histo_tf(data_table,
                                                      input = input,
                                                      output = output,
                                                      deg = TRUE)
  output$plotly_tf <- .render_histo_TF(reactive_histo_tf_deg)
  reactive_tf_table <- .prepare_reactive_histoTF_table(data_table,
                                                      input=input,
                                                      output=output,
                                                      deg=FALSE)
  reactive_tf_table_deg <- .prepare_reactive_histoTF_table(data_table,
                                                          input=input,
                                                          output=output,
                                                          deg=TRUE)
  output$download_csv_tf <- download_csv(deg=TRUE,
                                        tf=TRUE,
                                        plotType="histo",
                                        input=input,
                                        output=output,
                                        data_table=data_table)
  output$table_tf <- .render_histo_tf_table(reactive_tf_table_deg)

}


# ns('IntegrationSelect')
# paste0("input.", ns("IntegrationSelect"),"=='gene_genomic_res'")
# ns("genomicTypeSelect")
# ns('ClassSelect')
# ns('ChrSelect')
# ns('SignificativityCriteria')
# paste0("input.", ns("SignificativityCriteria"),"=='pval'")
# ns("PvalRange")
# paste0("input.", ns("SignificativityCriteria"),"=='FDR'")
# ns("FdrRange")



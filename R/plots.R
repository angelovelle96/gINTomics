
.circos_preprocess <- function(data){

  library(GenomicRanges)
  dataframes <- lapply(unique(data$omics), function(x) {
    single_omic_df <- data[data$omics==x,]
    return(single_omic_df)
  })
  names(dataframes) <- paste0("df_", unique(data$omics))
  if("df_gene_genomic_res"%in%names(dataframes)){
    dataframes$df_cnv <- filter(dataframes$df_gene_genomic_res,
                                cnv_met == 'cnv')
    dataframes$df_met <- filter(dataframes$df_gene_genomic_res,
                                cnv_met == 'met')
    tmp <- lapply(which(names(dataframes)!="df_gene_genomic_res"), function(x){
      ans <- dataframes[[x]]
      ans <- ans[ans$cov!="(Intercept)",]
      return(ans)
    })
    names(tmp) <- names(dataframes)[
      which(names(dataframes)!="df_gene_genomic_res")]
    dataframes <- tmp
  }
  gr <- lapply(dataframes, function(x){
    ans <- makeGRangesFromDataFrame(x,
                                    seqnames.field = 'chr',
                                    start.field = 'start',
                                    end.field = 'end',
                                    keep.extra.columns = TRUE,
                                    na.rm = TRUE)
    return(ans)
  })
  return(gr)
}

############################################################################
#############################################################################
.create_single_track <- function(data,
                                 dataValue,
                                 x_axis,
                                 xe_axis,
                                 y_axis,
                                 colorField,
                                 colorDomain,
                                 colorRange,
                                 tooltipField1,
                                 tooltipTitle,
                                 tooltipAlt1,
                                 tooltipField2,
                                 tooltipAlt2,
                                 tooltipField3,
                                 tooltipAlt3,
                                 tooltipField4,
                                 tooltipAlt4,
                                 legend) {
  return(
    add_single_track(
      data = track_data_gr(data, chromosomeField = 'seqnames', genomicFields = c('start','end'), value = dataValue),
      mark = 'bar',
      x = visual_channel_x(field = 'start', type = 'genomic', axis = x_axis),
      xe = visual_channel_x(field = 'end', type = 'genomic', axis = xe_axis),
      y = visual_channel_y(field = dataValue, type = 'quantitative', axis = y_axis),
      color = visual_channel_color(field = colorField, type = 'nominal', domain = colorDomain, range = colorRange),
      tooltip = visual_channel_tooltips(
        visual_channel_tooltip(field = "start", type = "genomic", alt = 'Start Position:'),
        visual_channel_tooltip(field = "end", type = "genomic", alt = "End Position:"),
        visual_channel_tooltip(field = tooltipField1, title = tooltipTitle, type = "quantitative", alt = paste(tooltipTitle, "Value:"), format = "0.2"),
        visual_channel_tooltip(field = tooltipField2, type = 'nominal', alt = tooltipAlt2),
        visual_channel_tooltip(field = tooltipField3, type = 'nominal', alt = tooltipAlt3),
        visual_channel_tooltip(field = tooltipField4, type = 'nominal', alt = tooltipAlt4)
      ),
      size = list(value = 1),
      legend = legend
    )
  )
}

#######################################################################
########################################################################

.create_tracks <- function(data_table, gr){

  tracks <- list()

  track_cyto <- add_single_track(
    id = "track2",
    data = track_data(
      url = "https://raw.githubusercontent.com/sehilyi/gemini-datasets/master/data/UCSC.HG38.Human.CytoBandIdeogram.csv",
      type = "csv",
      chromosomeField = "Chromosome",
      genomicFields = c("chromStart",
                        "chromEnd")
    ),
    mark = "rect",
    x = visual_channel_x(field = "chromStart",
                         type = "genomic"),
    xe = visual_channel_x(field = "chromEnd",
                          type = "genomic"),
    color = visual_channel_color(
      field = "Stain",
      type = "nominal",
      domain = c(
        "gneg",
        "gpos25",
        "gpos50",
        "gpos75",
        "gpos100",
        "gvar",
        "acen"           # acen: centromeric region (UCSC band files)
      ),
      range = c(
        "white",
        "#D9D9D9",
        "#979797",
        "#636363",
        "black",
        "#F0F0F0",
        "red"
      )
    ),
    stroke = visual_channel_stroke(
      value = "lightgray"
    ),
    strokeWidth = visual_channel_stroke_width(
      value = 0.5
    ),
  )

  if ("gene_genomic_res" %in% unique(data_table$omics)){
    #####da mettere in un'altra funzione
    track_cnv <- .create_single_track(data=gr$df_cnv,
                                     dataValue='cov_value',
                                     x_axis="none",
                                     xe_axis="none",
                                     y_axis="none",
                                     colorField="direction_cov",
                                     colorDomain=c("positive","negative"),
                                     colorRange=c("red","blue"),
                                     tooltipField1="cov_value",
                                     tooltipTitle="cnv",
                                     tooltipAlt1="CNV Value:",
                                     tooltipField2="gene",
                                     tooltipAlt2="Gene Name:",
                                     tooltipField3="class",
                                     tooltipAlt3="Class:",
                                     tooltipField4="cnv_met",
                                     tooltipAlt4="Integration Type:",
                                     legend=FALSE)

    track_met <- .create_single_track(data=gr$df_met,
                                     dataValue='cov_value',
                                     x_axis="none",
                                     xe_axis="none",
                                     y_axis="none",
                                     colorField="direction_cov",
                                     colorDomain=c("positive","negative"),
                                     colorRange=c("red","blue"),
                                     tooltipField1="cov_value",
                                     tooltipTitle="met",
                                     tooltipAlt1="MET Value:",
                                     tooltipField2="gene",
                                     tooltipAlt2="Gene Name:",
                                     tooltipField3="class",
                                     tooltipAlt3="Class:",
                                     tooltipField4="cnv_met",
                                     tooltipAlt4="Integration Type:",
                                     legend=FALSE)

    track_expr <- .create_single_track(data=gr$df_cnv,
                                      dataValue='response_value',
                                      x_axis="none",
                                      xe_axis="none",
                                      y_axis="none",
                                      colorField="direction_cov",
                                      colorDomain=c("positive","negative"),
                                      colorRange=c("red","blue"),
                                      tooltipField1="response_value",
                                      tooltipTitle="expr",
                                      tooltipAlt1="Expression Value (log2):",
                                      tooltipField2="gene",
                                      tooltipAlt2="Gene Name:",
                                      tooltipField3="class",
                                      tooltipAlt3="Class:",
                                      tooltipField4="cnv_met",
                                      tooltipAlt4="Integration Type:",
                                      legend=FALSE)
    tracks <- c(tracks, list(track_cnv=track_cnv, track_met=track_met, track_expr=track_expr, track_cyto=track_cyto))
  }


  if ("cnv_gene_res" %in% unique(gr$omics)){
    track_cnv <- .create_single_track(data=data$df_cnv_gene_res,
                                     dataValue='cov_value',
                                     x_axis="none",
                                     xe_axis="none",
                                     y_axis="none",
                                     colorField="direction_cov",
                                     colorDomain=c("positive","negative"),
                                     colorRange=c("red","blue"),
                                     tooltipField1="cov_value",
                                     tooltipTitle="cnv",
                                     tooltipAlt1="CNV Value:",
                                     tooltipField2="gene",
                                     tooltipAlt2="Gene Name:",
                                     tooltipField3="class",
                                     tooltipAlt3="Class:",
                                     tooltipField4="cnv_met",
                                     tooltipAlt4="Integration Type:",
                                     legend=FALSE)

    track_expr <- .create_single_track(data=gr$df_cnv_gene_res,
                                      dataValue='response_value',
                                      x_axis="none",
                                      xe_axis="none",
                                      y_axis="none",
                                      colorField="direction_cov",
                                      colorDomain=c("positive","negative"),
                                      colorRange=c("red","blue"),
                                      tooltipField1="response_value",
                                      tooltipTitle="expr",
                                      tooltipAlt1="Expression Value (log2):",
                                      tooltipField2="gene",
                                      tooltipAlt2="Gene Name:",
                                      tooltipField3="class",
                                      tooltipAlt3="Class:",
                                      tooltipField4="cnv_met",
                                      tooltipAlt4="Integration Type:",
                                      legend=FALSE)

    tracks <- c(tracks,list(track_cnv=track_cnv, track_expr=track_expr, track_cyto=track_cyto))
  }

  if("met_gene_res"%in%unique(data_table$omics)){

    track_met <- .create_single_track(data=gr$df_met_gene_res,
                                     dataValue='cov_value',
                                     x_axis="none",
                                     xe_axis="none",
                                     y_axis="none",
                                     colorField="direction_cov",
                                     colorDomain=c("positive","negative"),
                                     colorRange=c("red","blue"),
                                     tooltipField1="cov_value",
                                     tooltipTitle="met",
                                     tooltipAlt1="MET Value:",
                                     tooltipField2="gene",
                                     tooltipAlt2="Gene Name:",
                                     tooltipField3="class",
                                     tooltipAlt3="Class:",
                                     tooltipField4="cnv_met",
                                     tooltipAlt4="Integration Type:",
                                     legend=FALSE)

    track_expr <- .create_single_track(data=gr$df_met_gene_res,
                                      dataValue='response_value',
                                      x_axis="none",
                                      xe_axis="none",
                                      y_axis="none",
                                      colorField="direction_cov",
                                      colorDomain=c("positive","negative"),
                                      colorRange=c("red","blue"),
                                      tooltipField1="response_value",
                                      tooltipTitle="expr",
                                      tooltipAlt1="Expression Value (log2):",
                                      tooltipField2="gene",
                                      tooltipAlt2="Gene Name:",
                                      tooltipField3="class",
                                      tooltipAlt3="Class:",
                                      tooltipField4="response",
                                      tooltipAlt4="Integration Type:",
                                      legend=FALSE)

    tracks <- c(tracks, list(track_met=track_met, track_expr=track_expr, track_cyto=track_cyto))
  }

  return(tracks)
}

#######################################################################
########################################################################

.create_composed_view <- function(tracks, layout="circular") {

  composed_views <- list()

  if(sum(c("track_expr", "track_cnv", "track_met")%in%names(tracks))==3){

    composed_view_circos_genomic <- compose_view(
      multi = TRUE,
      layout = layout,
      tracks = add_multi_tracks(tracks$track_expr,
                                tracks$track_cnv,
                                tracks$track_met,
                                tracks$track_cyto),
      alignment = 'stack',
      spacing = 0.01,
      linkingId = "detail"
    )
    composed_views <- c(composed_views, list(circos_genomic=composed_view_circos_genomic))
  }else{

    if(!("track_cnv"%in%names(tracks))){

      composed_view_circos_met_gene <- compose_view(
        multi = TRUE,
        layout = layout,
        tracks = add_multi_tracks(tracks$track_expr,
                                  tracks$track_met,
                                  tracks$track_cyto),
        alignment = 'stack',
        spacing = 0.01,
        linkingId = "detail"
      )
      composed_views <- c(composed_views, list(circos_met_gene=composed_view_circos_met_gene))
    }
    if(!("track_met"%in%names(tracks))){

      composed_view_circos_cnv_gene <- compose_view(
        multi = TRUE,
        layout = layout,
        tracks = add_multi_tracks(tracks$track_expr,
                                  tracks$track_cnv,
                                  tracks$track_cyto),
        alignment = 'stack',
        spacing = 0.01,
        linkingId = "detail"
      )
      composed_views <- c(composed_views, list(circos_cnv_gene=composed_view_circos_cnv_gene))
    }

  }


  return(composed_views)
}

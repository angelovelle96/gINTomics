
.circos_preprocess <- function(data){

  library(GenomicRanges)
  dataframes <- lapply(unique(data$omics), function(x) {
    single_omic_df <- data[data$omics==x,]
    if(max(single_omic_df$response_value, na.rm = T)>30){
      single_omic_df$response_value <- log2(single_omic_df$response_value+1)
    }
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
      ans$cov_value <- ans$cov_value-min(ans$cov_value, na.rm = T)
      ans$cov_value <- ans$cov_value/max(ans$cov_value, na.rm = T)
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
                                 dataValue=NULL,
                                 x_axis=NULL,
                                 xe_axis=NULL,
                                 y_axis=NULL,
                                 colorField=NULL,
                                 colorDomain=NULL,
                                 colorRange=NULL,
                                 tooltipField1=NULL,
                                 tooltipTitle=NULL,
                                 tooltipAlt1=NULL,
                                 tooltipField2=NULL,
                                 tooltipAlt2=NULL,
                                 tooltipField3=NULL,
                                 tooltipAlt3=NULL,
                                 tooltipField4=NULL,
                                 tooltipAlt4=NULL,
                                 legend=NULL,
                                 colorType=NULL) {
  return(
    add_single_track(
      data = track_data_gr(data, chromosomeField = 'seqnames', genomicFields = c('start','end'), value = dataValue),
      mark = 'bar',
      x = visual_channel_x(field = 'start', type = 'genomic', axis = x_axis),
      xe = visual_channel_x(field = 'end', type = 'genomic', axis = xe_axis),
      y = visual_channel_y(field = dataValue, type = 'quantitative', axis = y_axis),
      color = visual_channel_color(field = colorField, type = colorType, domain = colorDomain, range = colorRange, legend = legend),
      tooltip = visual_channel_tooltips(
        visual_channel_tooltip(field = "start", type = "genomic", alt = 'Start Position:'),
        visual_channel_tooltip(field = "end", type = "genomic", alt = "End Position:"),
        visual_channel_tooltip(field = tooltipField1, title = tooltipTitle, type = "quantitative", alt = paste(tooltipTitle, "Value:"), format = "0.2"),
        visual_channel_tooltip(field = tooltipField2, type = 'nominal', alt = tooltipAlt2),
        visual_channel_tooltip(field = tooltipField3, type = 'nominal', alt = tooltipAlt3),
        visual_channel_tooltip(field = tooltipField4, type = 'nominal', alt = tooltipAlt4)
      ),
      size = list(value = 1)
    )
  )
}

#######################################################################
########################################################################

.create_tracks <- function(data, gr){

  tracks <- list()
  track_cyto <- .create_cyto_track()
  if ("gene_genomic_res" %in% unique(data$omics)){
    #####da mettere in un'altra funzione
    track_cnv <- .create_cnv_track(gr$df_cnv)
    track_met <- .create_met_track(gr$df_met)
    track_expr <- .create_expr_track(gr$df_cnv)
    track_coef_cnv <- .create_coef_track(gr$df_cnv)
    track_coef_met <- .create_coef_track(gr$df_met)
    tracks <- c(tracks, list(track_cnv=track_cnv,
                             track_met=track_met,
                             track_expr=track_expr,
                             track_coef_cnv=track_coef_cnv,
                             track_coef_met=track_coef_met))
  }

  if ("df_gene_cnv_res" %in% unique(names(gr))){
    track_cnv <- .create_cnv_track(gr$df_gene_cnv_res)
    track_expr <- .create_expr_track(gr$df_gene_cnv_res)
    tracks <- c(tracks,list(track_cnv=track_cnv, track_expr=track_expr))
  }

  if("df_gene_met_res"%in%unique(names(gr))){
    track_met <- .create_met_track(gr$df_gene_met_res)
    track_expr <- .create_expr_track(gr$df_gene_met_res)
    tracks <- c(tracks, list(track_met=track_met, track_expr=track_expr))
  }

    if ("df_mirna_cnv_res" %in% unique(names(gr))){
      track_mirna_cnv <- .create_cnv_track(gr$df_mirna_cnv_res)
      track_mirna_expr <- .create_expr_track(gr$df_mirna_cnv_res)
      tracks <- c(tracks,list(track_mirna_cnv=track_mirna_cnv,
                              track_mirna_expr=track_mirna_expr))
    }
  tracks <- c(tracks,list(track_cyto=track_cyto))
  return(tracks)
}

#######################################################################
########################################################################

.create_cnv_track <- function(gr){

    gr$cov_value2 <- as.character(gr$cov_value)
    ii <- cut(unique(as.numeric(gr$cov_value2)),
              breaks = seq(min(gr$cov_value, na.rm = T),
                           max(gr$cov_value, na.rm = T),
                           len = 50),
              include.lowest = TRUE)
    ccol <- colorRampPalette(c("#f2e6e6", "red4"))(49)[ii]
    track_cnv <- .create_single_track(data=gr,
                                      dataValue='cov_value',
                                      x_axis="none",
                                      xe_axis="none",
                                      y_axis="none",
                                      colorField="cov_value2",
                                      colorDomain=unique(gr$cov_value2),
                                      colorRange=ccol,
                                      tooltipField1="cov_value",
                                      tooltipTitle="cnv",
                                      tooltipAlt1="CNV Value:",
                                      tooltipField2="gene",
                                      tooltipAlt2="Gene Name:",
                                      tooltipField3="class",
                                      tooltipAlt3="Class:",
                                      tooltipField4="cnv_met",
                                      tooltipAlt4="Integration Type:",
                                      legend=FALSE,
                                      colorType="nominal")
    return(track_cnv)
    }

#######################################################################
########################################################################

.create_met_track <- function(gr){
  gr$cov_value2 <- as.character(gr$cov_value)
  ii <- cut(unique(as.numeric(gr$cov_value2)),
            breaks = seq(min(gr$cov_value, na.rm = T),
                         max(gr$cov_value, na.rm = T),
                         len = 50),
            include.lowest = TRUE)
  ccol <- colorRampPalette(c("#d1d1e3", "navyblue"))(49)[ii]
  track_met <- .create_single_track(data=gr,
                                    dataValue='cov_value',
                                    x_axis="none",
                                    xe_axis="none",
                                    y_axis="none",
                                    colorField="cov_value2",
                                    colorDomain=unique(gr$cov_value2),
                                    colorRange=ccol,
                                    tooltipField1="cov_value",
                                    tooltipTitle="met",
                                    tooltipAlt1="MET Value:",
                                    tooltipField2="gene",
                                    tooltipAlt2="Gene Name:",
                                    tooltipField3="class",
                                    tooltipAlt3="Class:",
                                    tooltipField4="cnv_met",
                                    tooltipAlt4="Integration Type:",
                                    legend=FALSE,
                                    colorType="nominal")
  return(track_met)

}

#######################################################################
########################################################################

.create_expr_track <- function(gr){

  track_expr <- .create_single_track(data=gr,
                                     dataValue='response_value',
                                     x_axis="none",
                                     xe_axis="none",
                                     y_axis="none",
                                     colorField="response_value",
                                     colorDomain=NULL,
                                     colorRange=NULL,
                                     tooltipField1="response_value",
                                     tooltipTitle="expr",
                                     tooltipAlt1="Expression Value (log2):",
                                     tooltipField2="gene",
                                     tooltipAlt2="Gene Name:",
                                     tooltipField3="class",
                                     tooltipAlt3="Class:",
                                     tooltipField4="cnv_met",
                                     tooltipAlt4="Integration Type:",
                                     legend=T,
                                     colorType="quantitative")
  return(track_expr)
}

#######################################################################
########################################################################

.create_coef_track <- function(gr){
  track_coef <- .create_single_track(data=gr,
                                    dataValue='coef',
                                    x_axis="none",
                                    xe_axis="none",
                                    y_axis="none",
                                    colorField="direction_coef",
                                    colorDomain=c("negative", "positive"),
                                    colorRange=c("blue", "red"),
                                    tooltipField1="coef",
                                    tooltipTitle="coef",
                                    tooltipAlt1="Coef Value:",
                                    tooltipField2="gene",
                                    tooltipAlt2="Gene Name:",
                                    tooltipField3="class",
                                    tooltipAlt3="Class:",
                                    tooltipField4="cnv_met",
                                    tooltipAlt4="Integration Type:",
                                    legend=FALSE,
                                    colorType="nominal")
  return(track_coef)
}



#######################################################################
########################################################################

.create_cyto_track <- function(){
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
  return(track_cyto)
}


#######################################################################
########################################################################

.create_composed_view <- function(tracks,width, height) {

  composed_views <- list()

  if(sum(c("track_expr", "track_cnv", "track_met")%in%names(tracks))==3){

    composed_view_circos_genomic <- compose_view(
      multi = TRUE,
      tracks = add_multi_tracks(tracks$track_cyto,
                                tracks$track_expr,
                                tracks$track_cnv,
                                tracks$track_coef_cnv,
                                tracks$track_met,
                                tracks$track_coef_met),
      alignment = 'stack',
      spacing = 0.01,
      linkingId = "detail",
      width = width,
      height = height
    )
    composed_views <- c(composed_views, list(circos_genomic=composed_view_circos_genomic))
  }else{

    if("track_met"%in%names(tracks)){

      composed_view_circos_met_gene <- compose_view(
        multi = TRUE,
        tracks = add_multi_tracks(tracks$track_expr,
                                  tracks$track_met,
                                  tracks$track_cyto),
        alignment = 'stack',
        spacing = 0.01,
        linkingId = "detail"
      )
      composed_views <- c(composed_views, list(circos_met_gene=composed_view_circos_met_gene))
    }
    if("track_cnv"%in%names(tracks)){

      composed_view_circos_cnv_gene <- compose_view(
        multi = TRUE,
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
  if(sum(c("track_mirna_expr", "track_mirna_cnv")%in%names(tracks))==2){

    composed_view_circos_genomic_mirna <- compose_view(
      multi = TRUE,
      tracks = add_multi_tracks(tracks$track_cyto,
                                tracks$track_mirna_expr,
                                tracks$track_mirna_cnv),
      alignment = 'stack',
      spacing = 0.01,
      linkingId = "detail",
      width = width,
      height = height
    )
    composed_views <- c(composed_views,
                        list(circos_genomic_mirna=composed_view_circos_genomic_mirna))
  }
  return(composed_views)
}

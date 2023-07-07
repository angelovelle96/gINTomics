#'
#' @import multiMiR
#' @export


   download_mirna_target <- function(miRNAs,
                                     targets = NULL,
                                     species = 'hsa',
                                     table   = 'mirtarbase',
                                     ...){


     tmp <- get_multimir(org = species,
                         target = targets,
                         mirna   = miRNAs,
                         table   = table,
                         summary = F,
                         ...)
     tmp2 <- tmp@data[grep("Luciferase", tmp@data$experiment),]
     ans <- lapply(unique(tmp2$target_symbol), function(x)
       tmp2$mature_mirna_id[tmp2$target_symbol==x])

     names(ans) <- unique(tmp2$target_symbol)
     ans <- lapply(ans, unique)
     return(ans)
   }


#'
#' @import  dorothea
#' @export

   download_tf <- function(genes,
                           species="hsa",
                           include_predicted=F){

     if(!species%in%c("hsa", "mmu")) stop(str_wrap("Available species are hsa
                                                   and mmu for TF download"))

     if(species=="hsa"){
       data("dorothea_hs", package = "dorothea")
       data <- dorothea_hs
       if(!include_predicted) data <- data[data$confidence!="E",]
       data <- data[data$tf%in%genes,]
       data <- data[data$target%in%genes,]

      }

     if(species=="mmu"){
       data("dorothea_mm", package = "dorothea")
       data <- dorothea_mm
       if(!include_predicted) data <- data[data$confidence!="E",]
       data <- data[data$tf%in%genes,]
       data <- data[data$target%in%genes,]
     }

   ans <- lapply(unique(data$target), function(x)
     data$tf[data$target==x])

   names(ans) <- unique(data$target)
   return(ans)


}


#' @import biomaRt
#' @importFrom stringr str_wrap
#' @export

   download_gene_info <- function(genes,
                                  species = "hsa",
                                  filters = "hgnc_symbol"){

     tmp <- setNames(c("hsapiens_gene_ensembl",
                       "mmusculus_gene_ensembl",
                       "drerio_gene_ensembl",
                       "dmelanogaster_gene_ensembl",
                       "celegans_gene_ensembl",
                       "rnorvegicus_gene_ensembl"),
                     c("hsa", "mmu", "dre", "dme", "cel", "rno"))
     if(!species%in%names(tmp)) stop(str_wrap(paste("Supported species for gene
                                                    annotation are:",
                                                    paste(names(tmp),
                                                          sep = " ",
                                                          collapse = ", "))))

     dataset <- tmp[species]
     genes <- unique(genes)
     mart <- useMart(biomart = "ensembl",
                     dataset = dataset)
     ans <- getBM(attributes = c('hgnc_symbol',
                                 "ensembl_gene_id",
                                 "chromosome_name",
                                 "start_position",
                                 "end_position",
                                 "band",
                                 "description"),
                  mart = mart,
                  filters = filters,
                  values = genes)


       tmp <- c(1:100, "X", "Y", "MT")
       ans <- ans[ans$chromosome_name%in%tmp,]

       ans <- ans[!duplicated(ans$hgnc_symbol),]
       ans <- ans[ans$hgnc_symbol!="",]
       rownames(ans) <- ans$hgnc_symbol


       data("miRBase_conversion")
       genes <- cbind(genes, nop = gsub("-3p", "", genes))
       genes[,"nop"] <- gsub("-5p", "", genes[, "nop"])

       if(sum(mirna_hsa$Alias.symbols%in%genes[, "nop"])>0){

         mirna_hsa <- mirna_hsa[mirna_hsa$Alias.symbols%in%genes[, "nop"],]
         mirna_hsa <- mirna_hsa[!duplicated(mirna_hsa$Approved.symbol),]
         rownames(mirna_hsa) <- mirna_hsa$Approved.symbol
         tmp <- setNames(mirna_hsa$Alias.symbols, mirna_hsa$Approved.symbol)
         tmp <- strsplit(tmp, split = ",")
         tmp <- sapply(tmp, function(x) x[1])
         tmp <- tmp[!duplicated(tmp)]
         tmp <- tmp[!is.na(tmp)]
         tmp2 <- getBM(attributes = c('hgnc_symbol',
                                     "ensembl_gene_id",
                                     "chromosome_name",
                                     "start_position",
                                     "end_position",
                                     "band",
                                     "description"),
                      mart = mart,
                      filters = "hgnc_symbol",
                      values = names(tmp))
         tmp3 <- c(1:100, "X", "Y", "MT")
         tmp2 <- tmp2[tmp2$chromosome_name%in%tmp3,]
         tmp2 <- tmp2[!duplicated(tmp2$hgnc_symbol),]
         rownames(tmp2) <- tmp2$hgnc_symbol
         tmp <- tmp[intersect(rownames(tmp2), names(tmp))]
         tmp2 <- tmp2[names(tmp),]
         tmp2$miRBase_id <- tmp
         tmp2 <- tmp2[!duplicated(tmp2$miRBase_id),]
         tmp2$original_hgnc_symbol <- tmp2$hgnc_symbol
         tmp2$hgnc_symbol <- tmp2$miRBase_id
         tmp2 <- tmp2[, colnames(tmp2)!="miRBase_id"]
         rownames(tmp2) <- tmp2$hgnc_symbol
         tmp2 <- tmp2[genes[,"nop"],]
         rownames(tmp2) <- genes[, "genes"]
         tmp2 <- tmp2[!is.na(tmp2$hgnc_symbol),]
         ans$original_hgnc_symbol <- ans$hgnc_symbol
         ans <- ans[!ans$hgnc_symbol%in%tmp2$original_hgnc_symbol,]
         ans <- rbind(ans, tmp2)
       }


       return(ans)
   }






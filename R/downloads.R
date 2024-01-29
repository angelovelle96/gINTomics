#'
#' @import OmnipathR miRBaseConverter

   .download_mirna_target <- function(miRNAs,
                                     species = 'hsa',
                                     ...){
     if(!species%in%c("hsa", "mmu", "rno")) stop(str_wrap("Available species
                                                          are hsa, mmu and rno
                                                          for miRNA target
                                                          download"))
     org <- setNames(c(9606, 10090, 10116), c("hsa", "mmu", "rno"))
     mirna_target <-  import_mirnatarget_interactions(organism = org[species])
     tmp <- getAllMiRNAs(species = species)
     rownames(tmp) <- tmp$Accession
     mirna_target$mirbase <- tmp[mirna_target$source, "Name"]
     mirna_target$mirbase[is.na(mirna_target$mirbase)] <-
       mirna_target$source_genesymbol[is.na(mirna_target$mirbase)]
     names(miRNAs) <- tolower(miRNAs)
     mirna_target <- mirna_target[
       tolower(mirna_target$mirbase)%in%names(miRNAs),]
     mirna_target$miRNAs <- miRNAs[tolower(mirna_target$mirbase)]
     ans <- lapply(unique(mirna_target$target_genesymbol), function(x)
       mirna_target$miRNAs[mirna_target$target_genesymbol==x])
     names(ans) <- unique(mirna_target$target_genesymbol)
     ans <- lapply(ans, unique)
     return(ans)
   }


#'
#' @import OmnipathR miRBaseConverter

    .download_tf_mirna <- function(miRNAs,
                                      species = 'hsa',
                                      ...){

      if(!species%in%c("hsa", "mmu", "rno")) stop(str_wrap("Available species
                                                          are hsa, mmu and rno
                                                          for miRNA target
                                                          download"))
      org <- setNames(c(9606, 10090, 10116), c("hsa", "mmu", "rno"))
      tf_mirna <-  import_tf_mirna_interactions(organism = org[species], resources = "TransmiR")
      tmp <- getAllMiRNAs(species = species)
      rownames(tmp) <- tmp$Accession
      tf_mirna$mirbase <- tmp[tf_mirna$target, "Name"]
      tf_mirna$mirbase[is.na(tf_mirna$mirbase)] <-
        tf_mirna$target_genesymbol[is.na(tf_mirna$mirbase)]
      tmp <- tf_mirna
      tmp$mirbase <- gsub("-3p", "", tmp$mirbase)
      tmp$mirbase <- gsub("-5p", "", tmp$mirbase)
      tmp$mirbase <- gsub("\\*", "", tmp$mirbase)
      tf_mirna <- rbind(tf_mirna, tmp)
      tf_mirna$check <- paste(tf_mirna$source_genesymbol,
                              tf_mirna$mirbase,
                              sep = "_")
      tf_mirna <- tf_mirna[!duplicated(tf_mirna$check),]
      names(miRNAs) <- tolower(miRNAs)
      tf_mirna <- tf_mirna[
        tolower(tf_mirna$mirbase)%in%names(miRNAs),]
      tf_mirna$miRNAs <- miRNAs[tolower(tf_mirna$mirbase)]
      ans <- lapply(unique(tf_mirna$miRNAs), function(x)
        tf_mirna$source_genesymbol[tf_mirna$miRNAs==x])
      names(ans) <- unique(tf_mirna$miRNAs)
      ans <- lapply(ans, unique)
      return(ans)


    }




#'
#' @import  OmnipathR

   .download_tf <- function(genes,
                           species="hsa",
                           ...){

     if(!species%in%c("hsa", "mmu", "rno")) stop(str_wrap("Available species
                                                          are hsa, mmu and rno
                                                          for TF target
                                                          download"))
     org <- setNames(c(9606, 10090, 10116), c("hsa", "mmu", "rno"))
     tf_target <-  import_tf_target_interactions(organism = org[species])
     tf_target <- tf_target[tf_target$source_genesymbol%in%genes,]
     tf_target <- tf_target[tf_target$target_genesymbol%in%genes,]
     ans <- lapply(unique(tf_target$target_genesymbol), function(x)
       tf_target$source_genesymbol[tf_target$target_genesymbol==x])
     names(ans) <- unique(tf_target$target_genesymbol)
     ans <- lapply(ans, unique)
     return(ans)


}


#' @import biomaRt
#' @importFrom stringr str_wrap
#' @importFrom stats setNames
#' @importFrom plyr rbind.fill

   .download_gene_info_biomart <- function(genes,
                                          species = "hsa",
                                          filters = c("hgnc_symbol",
                                                      "ensembl_gene_id",
                                                      "entrezgene_id"),
                                          ...){

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
     ans <- lapply(filters, function(x){
       ans <- getBM(attributes = c('hgnc_symbol',
                                   "ensembl_gene_id",
                                   "entrezgene_id",
                                   "chromosome_name",
                                   "start_position",
                                   "end_position",
                                   "band"),
                    mart = mart,
                    filters = x,
                    values = genes,
                    ...)
       tmp <- c(1:100, "X", "Y", "MT")
       ans <- ans[ans$chromosome_name%in%tmp,]
       ans <- ans[!duplicated(ans[,x]),]
       ans <- ans[ans[,x]!="",]
       ans$names <- ans[,x]
       return(ans)
       })
     names(ans) <- filters
     ans <- rbind.fill(ans)
     rownames(ans) <- ans$names
     ans <- ans[, colnames(ans)!="names"]

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
                                    "entrezgene_id",
                                    "chromosome_name",
                                    "start_position",
                                    "end_position",
                                    "band"),
                    mart = mart,
                    filters = "hgnc_symbol",
                    values = names(tmp),
                    ...)
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


   #' @import Homo.sapiens Mus.musculus
   #' @importFrom stringr str_wrap
   #' @importFrom stats setNames
   #' @importFrom plyr rbind.fill
   #' @importFrom GenomicFeatures genes
   #' @importFrom dplyr mutate_all
   #' @importFrom BiocGenerics which.max

   .download_gene_info_org <- function(genes,
                                   species = "hsa",
                                   ...){

     tmp <- c("hsa", "mmu")
     if(!species%in%tmp) stop(str_wrap(paste("Supported species for gene
                                                    annotation are:",
                                                    paste(tmp,
                                                          sep = " ",
                                                          collapse = ", "))))
     if(species=="hsa") dataset <- Homo.sapiens
     if(species=="mmu") dataset <- Mus.musculus
     genes <- unique(genes)
     all_genes <- suppressMessages(genes(dataset,
                                         columns=c("SYMBOL",
                                                   "MAP",
                                                   "ENSEMBL")))
     all_genes <- as.data.frame(all_genes)
     all_genes <- mutate_all(all_genes, function(x)
       unlist(lapply(x, function(y) y[1])))
     all_genes$ENTREZ <- rownames(all_genes)
     pos <- apply(all_genes, 2, function(x) length(intersect(genes,x)))
     ans <- data.frame()
     if(sum(pos)>0){
       pos <- which(pos!=0)
       ans <- lapply(pos, function(x){
                 ans <- all_genes[!is.na(all_genes[,pos]),]
                 ans <- ans[!duplicated(ans[,pos]),]
                 rownames(ans) <- ans[,pos]
                 return(ans)
                 })
       ans <- rbind.fill(ans)
       }

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

       tmp2 <- all_genes
       tmp2 <- tmp2[!duplicated(tmp2$SYMBOL),]
       tmp2 <- tmp2[!is.na(tmp2$SYMBOL),]
       rownames(tmp2) <- tmp2$SYMBOL
       tmp <- tmp[intersect(rownames(tmp2), names(tmp))]
       tmp2 <- tmp2[names(tmp),]
       tmp2$miRBase_id <- tmp
       tmp2 <- tmp2[!duplicated(tmp2$miRBase_id),]
       rownames(tmp2) <- tmp2$miRBase_id
       tmp2 <- tmp2[, colnames(tmp2)!="miRBase_id"]
       tmp2 <- tmp2[genes[,"nop"],]
       rownames(tmp2) <- genes[, "genes"]
       tmp2 <- tmp2[!is.na(tmp2$SYMBOL),]
       ans <- rbind(ans, tmp2)
     }

     tmp <- intersect(rownames(ans), genes[, "genes"])
     ans <- ans[tmp,]
     colnames(ans) <- gsub("seqnames", "chromosome_name", colnames(ans))
     colnames(ans) <- gsub("start", "start_position", colnames(ans))
     colnames(ans) <- gsub("end", "end_position", colnames(ans))
     colnames(ans) <- gsub("ENSEMBL", "ensembl_gene_id", colnames(ans))
     colnames(ans) <- gsub("MAP", "band", colnames(ans))
     colnames(ans) <- gsub("SYMBOL", "hgnc_symbol", colnames(ans))
     colnames(ans) <- gsub("ENTREZ", "entrezgene_id", colnames(ans))
     ans <- ans[, c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id",
                    "chromosome_name", "start_position",
                    "end_position", "band")]
     ans$chromosome_name <- gsub("chr", "", ans$chromosome_name)
     ans$band <- gsub(".*q", "q", ans$band)
     ans$band <- gsub(".*p", "p", ans$band)
     return(ans)
   }



   .download_gene_info <- function(genes,
                                   species = "hsa",
                                   biomaRt = F,
                                   ...){

     if(biomaRt){
       ans <- .download_gene_info_biomart(genes = genes,
                                         species = species,
                                         ...)
     }else{
       ans <- .download_gene_info_org(genes = genes,
                                     species = species,
                                     ...)
     }
     return(ans)
   }


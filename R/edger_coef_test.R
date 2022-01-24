edger_coef_test <- function(fit, threads = 1) {

    top_list <- mclapply(fit, function(x) {
        coef <- colnames(x$coefficients)
        lrt <- lapply(coef, function(z) glmLRT(x,
            coef = z))
        names(lrt) <- coef
        top <- as.matrix(sapply(names(lrt),
            function(z) topTags(lrt[[z]])$table))
        coef <- sapply(names(lrt), function(z) lrt[[z]]$coefficients[,
            z])
        top <- rbind(top, coef)
        top <- as.data.frame(top)
        return(top)

    }, mc.cores = threads)

    return(top_list)
}

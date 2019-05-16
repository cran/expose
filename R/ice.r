#' Extract the information from the simulation data frame to analyse the individual conditional expectation
#'
#' @param allsim dataset with all simulations values
#' @param dataset dataset with all variables
#' @param dr a vector with dose response values
#' @param squem squeme of the values of the prediction values
#' @param remove_extrem boolean parameter to remove the extrem values
#' @return a data frame with interactions
#' @examples
#' data(expose_data)
#' data(simu)
#' data(gen)
#' delta=c(1,0)
#' seku <- seq(0,1,0.05)
#' Exposures<- c('Var1','Var2','Var3','Var4','Var5')
#' summary_table_lines <- gen[[2]]
#' ice_res <- ice(allsim = simu[[1]], dataset = expose_data, dr = seku, 
#' squem = summary_table_lines, remove_extrem = FALSE)
#' @export

ice <- function(allsim, dataset, dr = seq(0, 1, 0.1), squem, remove_extrem = FALSE) {
    
    ext <- squem[grep("DR_", squem$Group), ]
    
    hm_pol <- unique(ext$Case)
    
    n <- nrow(dataset)
    
    ext2 <- ext[ext$Case == hm_pol[1], ]
    
    df <- data.frame(allsim)
    
    ini <- ext2$From[1]
    fin <- ext2$To[nrow(ext2)]
    
    allpcb <- data.frame(df[ini:fin, 1])
    
    allpcb$id <- rep(1:nrow(dataset), length(dr))
    
    mm <- NULL
    len2 <- length(dr) - 1
    for (g in 0:len2) {
        mm <- c(mm, rep(g, n))
    }
    
    allpcb$X <- mm
    names(allpcb)[1] <- paste0(hm_pol[1], "_pred")
    
    if (length(hm_pol) > 1) {
        hm_pol2 <- hm_pol[-1]
        for (h in hm_pol2) {
            # extracting again for each pollutant
            ext2 <- ext[ext$Case == h, ]
            
            df <- data.frame(allsim)
            
            ini <- ext2$From[1]
            fin <- ext2$To[nrow(ext2)]
            
            allpcb2 <- data.frame(df[ini:fin, 1])
            
            allpcb2$id <- rep(1:nrow(dataset), length(dr))
            
            allpcb2$X <- mm
            names(allpcb2)[1] <- paste0(h, "_pred")
            
            allpcb <- merge(allpcb, allpcb2, by = c("id", "X"))
            
        }
    }
    
    ice <- merge(dataset, allpcb, by = "id")
    
    seku_breaks <- as.character(sort(as.numeric(as.character(unique(ice$X)))))
    
    ice$X <- factor(ice$X, levels = seku_breaks)
    
    if (remove_extrem) {
        sekut <- dr[c(1, length(dr))]
        ice <- ice[!ice$X %in% sekut, ]
    }
    
    return(ice)
}

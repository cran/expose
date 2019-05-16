#' Extract the information from the simulation data frame to analyse the dose response effects
#'
#' @param allsim dataset with all simulations values
#' @param dataset dataset with all variables
#' @param exposures a vector with exposures
#' @param dr a vector with dose response values
#' @param ic_dis choose between ic (interval confidences) and dis (distribution)
#' @param st summary table from general function
#' @return a data frame with dose response values
#' @examples
#' data(expose_data)
#' data(simu)
#' data(gen)
#' delta=c(1,0)
#' seku <- seq(0,1,0.05)
#' Exposures<- c('Var1','Var2','Var3','Var4','Var5')
#' summary_table_lines <- gen[[2]]
#' drr.grp <- dose_resp (allsim = simu[[1]], dataset = expose_data, st = summary_table_lines,
#'                       dr = seku, exposures = Exposures) 
#' @export


dose_resp <- function(allsim, dataset, exposures, dr, ic_dis = "IC", st) {
    st2 <- as.character(st[, 1])
    un_dr <- unique(grep("DR_", st2, value = TRUE))
    len_exp <- length(exposures)
    len_dr <- length(un_dr)
    df_ace <- data.frame(matrix(NA, len_exp * len_dr, 4))
    names(df_ace) <- c("Group", "Mean", "ICa", "ICb")
    b2 <- data.frame(matrix(NA, 1, 4))
    names(b2) <- c("Quantile", "Mean", "SE", "Exp")
    h <- 0
    for (ex in 1:len_exp) {
        h <- h + 1
        stp <- st[st$Case == exposures[ex], c(4, 5)]
        from <- as.numeric(stp[3, 1])
        to <- as.numeric(stp[dim(stp)[1], 2])
        mdata <- allsim[from:to, ]
        b <- dose_resp_ind(allsim = mdata, dataset = dataset, dr = dr)
        b$Exp <- exposures[ex]
        b2 <- rbind(b2, b)
    }
    b2 <- b2[-1, ]
    return(b2)
}

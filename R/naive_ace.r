#' Extract the information from the simulation data frame to analyse the naive causal effects
#'
#' @param allsim dataset with all simulations values
#' @param dataset dataset with all variables
#' @param exposures a vector with exposures
#' @param delta a vector with two values
#' @param ic_dis choose between ic (interval confidences) and dis (distribution)
#' @param st summary table from general function
#' @return a data frame with naive ace and confident intervals
#' @examples
#' data(expose_data)
#' data(simu)
#' data(gen)
#' delta=c(1,0)
#' Exposures<- c('Var1','Var2','Var3','Var4','Var5')
#' summary_table_lines <- gen[[2]]
#' ace.df.g <- naive_ace (allsim = simu[[1]], dataset = expose_data,
#' ic_dis = 'IC', st = summary_table_lines,
#' exposures = Exposures, delta = delta)
#' @export


naive_ace <- function(allsim, dataset, exposures, delta = c(0, 1), ic_dis = "IC", st) {
    len_exp <- length(exposures)
    df_ace <- data.frame(matrix(NA, len_exp, 4))
    names(df_ace) <- c("Group", "Mean", "ICa", "ICb")
    h <- 0
    for (ex in 1:len_exp) {
        h <- h + 1
        df_ace[h, "Group"] <- paste0(exposures[ex])
        stp <- st[st$Group == "ACE" & st$Case == exposures[ex], c(4, 5)]
        from <- as.numeric(stp[1, 1])
        to <- as.numeric(stp[2, 2])
        mdata <- allsim[from:to, ]
        b <- naive_ace_ind(allsim = mdata, dataset = dataset, ic_dis = "IC")
        df_ace[h, "Mean"] <- b[[1]]
        df_ace[h, "ICa"] <- b[[2]]
        df_ace[h, "ICb"] <- b[[3]]
    }
    return(df_ace)
}

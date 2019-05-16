#' Extract the information from the simulation data frame to analyse the dose response effects
#'
#' @param allsim dataset with all simulations values
#' @param dataset dataset with all variables
#' @param dr a vector with dose response values
#' @return a data frame with dose response values
#' @export

dose_resp_ind <- function(allsim, dataset, dr = seq(0, 1, 0.1)) {
    N <- dim(dataset)[1]
    sim <- dim(allsim)[2]
    size <- length(dr)
    m2 <- NA
    sum_dr <- data.frame(matrix(NA, size, 3))
    names(sum_dr) <- c("Quantile", "Mean", "SE")
    pos <- 1
    for (g in 1:size) {
        newpos <- pos + N - 1
        m3 <- as.matrix(allsim[pos:newpos, ])
        for (f in 1:sim) {
            m1 <- mean(as.matrix(allsim[pos:newpos, f]))
            m2 <- c(m2, m1)
        }
        m2 <- m2[-1]
        se <- stats::sd(m2)
        sum_dr[g, 2] <- mean(m3)
        sum_dr[g, 3] <- se
        sum_dr[g, 1] <- paste0("DR_", dr[g])
        pos <- newpos + 1
    }
    return(sum_dr)
}

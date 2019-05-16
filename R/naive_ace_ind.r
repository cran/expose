#' Extract the information from the simulation data frame to analyse the interaction effects
#'
#' @param allsim dataset with all simulations values
#' @param dataset dataset with all variables
#' @param ic_dis choose between ic (interval confidences) and dis (distribution)
#' @return a data frame with naive ace and confident intervals
#' @export

naive_ace_ind <- function(allsim, dataset, ic_dis = "IC") {
    N <- dim(dataset)[1]
    sim <- dim(allsim)[2]
    m2 <- NA
    mean.res <- mean(as.matrix(allsim[1:N, ])) - mean(as.matrix(allsim[(N + 1):(N * 2), ]))
    for (f in 1:sim) {
        m1 <- mean(as.matrix(allsim[1:N, f])) - mean(as.matrix(allsim[(N + 1):(N * 2), f]))
        m2 <- c(m2, m1)
    }
    m2 <- m2[-1]
    if (ic_dis == "IC") {
        sdd <- stats::sd(m2[-1])
        lci <- mean.res - sdd * 1.96
        uci <- mean.res + sdd * 1.96
    } else if (ic_dis == "dis") {
        lci <- stats::quantile(m2, 0.025)
        uci <- stats::quantile(m2, 0.975)
    }
    return(list(mean.res, lci, uci))
}

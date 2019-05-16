#' Extract the information from the simulation data frame to analyse the interaction effects
#'
#' @param allsim dataset with all simulations values
#' @param dataset dataset with all variables
#' @param exposures a vector with exposures
#' @param confounders a vector with confounders
#' @param squem squeme of the values of the prediction values
#' @return data frame with interaction values
#' @examples
#' data(expose_data)
#' data(simu)
#' data(gen)
#' delta=c(1,0)
#' seku <- seq(0,1,0.05)
#' Exposures<- c('Var1','Var2','Var3','Var4','Var5')
#' summary_table_lines <- gen[[2]]
#' it <- interact (allsim = simu[[1]], dataset = expose_data,exposures = Exposures,
#' confounders = c('sex'), squem = summary_table_lines)
#' @export

interact <- function(allsim, dataset, exposures, confounders, squem) {
    
    N <- dim(dataset)[1]
    sim <- dim(allsim)[2]
    
    m2 <- NA
    len_cof <- length(confounders)
    len_exp <- length(exposures)
    len_tot <- len_exp + len_cof
    var_tot <- c(exposures, confounders)
    
    ext2 <- squem[grep("ITE_", squem$Group), ]
    
    ini <- ext2$From[1]
    fin <- ext2$To[nrow(ext2)]
    
    allsim <- allsim[ini:fin, ]
    
    h <- 0
    n <- 1
    z <- len_exp + len_cof
    z2 <- (z * (z - 1)) / 2
    
    azken_taula <- data.frame(matrix(NA, z2, sim + 1))
    names(azken_taula)[1] <- c("Interaction")
    ec <- 1
    h <- 1
    for (ex in 1:len_tot) {
        ec <- ec + 1
        if (len_tot == ex) 
            break
        for (ex2 in ec:len_tot) {
            inter2 <- paste0(var_tot[ex], "-", var_tot[ex2])
            azken_taula[h, 1] <- inter2
            h <- h + 1
        }
    }
    ec <- 1
    h3 <- 1
    gehi <- 0
    for (ex in 1:len_tot) {
        ec <- ec + 1
        if (len_tot == ex) 
            break
        for (ex2 in ec:len_tot) {
            inter2 <- paste0(var_tot[ex], "-", var_tot[ex2])
            for (si in 1:sim) {
                n <- 1 + gehi
                n2 <- n + N
                a <- allsim[n:(n2 - 1), si]
                n <- n2
                n2 <- n + N
                b <- allsim[n:(n2 - 1), si]
                n <- n2
                n2 <- n + N
                
                c <- allsim[n:(n2 - 1), si]
                n <- n2
                n2 <- n + N
                
                d <- allsim[n:(n2 - 1), si]
                n <- n2
                azken_taula[h3, si + 1] <- mean(a - b - c + d)
            }
            gehi <- gehi + 4 * N
            h3 <- h3 + 1
        }
    }
    azken_taula2 <- data.frame(matrix(NA, z2, 3))
    m1 <- apply(azken_taula[, -1], 1, mean)
    m2 <- apply(azken_taula[, -1], 1, stats::sd)
    azken_taula2 <- data.frame(azken_taula$Interaction, m1, m2)
    names(azken_taula2) <- c("Interaction", "Mean", "SD")
    return(azken_taula2)
}

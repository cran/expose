#' General function create the basic structure for the analysis
#'
#' @param dataset dataset with all variables
#' @param exposures a vector with exposures
#' @param confounders a vector with confounders
#' @param outcomes outcome´s name
#' @param delta a vector with two values
#' @param dr a vector with dose response values
#' @return a list with 2 objects. One is the dataframe with all the values and the other is a summary of the groups and the corresponding rows in the first dataframe.
#' @examples
#' data(expose_data)
#' N <- dim(expose_data)[1]
#' Outcome='Y4'
#' seku <- seq(0,1,0.05)   #c(0,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,1)
#' our.num.sim <- 5
#' delta=c(1,0)
#' Exposures<- c('Var1','Var2','Var3','Var4','Var5')
#' Confounders<- c('sex')
#' Outcome <- c('Y4')
#' gen <- general_function (dataset = expose_data, exposures = Exposures, 
#'                          confounders = Confounders,
#'                          outcomes = Outcome[1], delta=delta, dr = seku)
#' @export
#' @details .

general_function <- function(dataset, exposures,
                             confounders, outcomes, delta, dr) {

    if (class(exposures) %in% c("numeric", "integer")) {
        exposures <- names(dataset[, exposures])
    }
    if (class(confounders) %in% c("numeric", "integer")) {
        confounders <- names(dataset[, confounders])
    }
    if (class(outcomes) %in% c("numeric", "integer")) {
        outcomes <- names(dataset[, outcomes])
    }

    len_exp <- length(exposures)
    len_cof <- length(confounders)
    var_tot <- c(exposures, confounders)
    len_tot <- len_exp + len_cof

    mid_data <- dataset[, c(exposures, confounders)]

    N <- dim(dataset)[1]

    len_delta <- length(delta)
    newdata <- rep(NA, len_tot)
    len_exp2 <- len_exp * len_delta


    for (bat in 1:len_exp2) {
        newdata <- rbind(newdata, mid_data)
    }
    newdata <- newdata[-1, ]

    h <- 0
    n <- 1
    laburpencc <- NA
    for (bat in 1:len_exp) {
        for (bi in 1:len_delta) {
            h <- h + 1
            n2 <- n + N
            newdata[n:(n2 - 1), bat] <- delta[bi]
            laburpencc <- c(laburpencc, "ACE", exposures[bat], 
                            delta[bi], n, (n2 - 1))
            n <- n2
        }
    }
    laburpencc <- laburpencc[-1]
    laburpencc2 <- (data.frame(matrix(laburpencc, ncol = 5, byrow = TRUE)))

    zenbatcc <- n2 - 1

    sekuentzia <- dr
    len_sek <- length(sekuentzia)

    quar <- paste0("newdata2")
    assign(quar, data.frame(matrix(NA, 1, len_tot)))
    quary <- "newdata_sek"
    assign(quary, data.frame(matrix(NA, 1, len_tot)))
    assign("tartey", get(quary))
    names(tartey) <- names(dataset[, var_tot])
    assign(quary, tartey)

    len_sek_exp <- len_sek * len_exp

    for (r in 1:len_sek) {
        tarteko <- get("newdata2")
        names(tarteko) <- names(dataset[, var_tot])
        assign(quar, rbind(tarteko, dataset[, var_tot]))
    }

    assign(quar, get(quar)[-1, ])

    per_value <- NA


    h <- 0
    for (s in 1:len_exp) {
        for (f in sekuentzia) {
            h <- h + 1
            valor <- dataset[, exposures[s]]
            per_value <- c(per_value, as.numeric(stats::quantile(valor, f)))
       
        }
        
    }
    
    
    per_value <- per_value[-1]
    
    h <- 0
    n <- 1
    vec <- n
    for (bat in 1:len_sek) {
        n2 <- n + N
        vec <- c(vec, (n2 - 1))
        n <- n2
        vec <- c(vec, (n2))
    }
    vec <- vec[-length(vec)]
    
    z2gehi <- vec[length(vec)]
    
    laburpendr <- NA
    h <- 0
    
    zgehi <- 0
    for (s in 1:len_exp) {
        quar1 <- paste0("newdata2_", exposures[s])
        tarteko <- get("newdata2")
        names(tarteko) <- names(dataset[, var_tot])
        assign(quar1, tarteko)
        
        g <- 0
        h <- 0
        
        for (f in 1:len_sek) {
            g <- g + 1
            h1 <- h + 1
            h2 <- h1 + 1
            
            tarteko2 <- get(quar1)
            tarteko2[vec[h1]:vec[h2], exposures[s]] <- per_value[g]
            laburpendr <- c(laburpendr, paste0("DR_", f), exposures[s], 
                            f, vec[h1] + zgehi, vec[h2] + zgehi)
            assign(quar1, tarteko2)
            
            h <- h2
            
        }
        
        zgehi <- zgehi + z2gehi
        
        assign(quary, rbind(get(quary), get(quar1)))
    }
    
    laburpendr <- laburpendr[-1]
    laburpendr2 <- data.frame(matrix(laburpendr, ncol = 5, byrow = TRUE))
    
    new_data3 <- get(quary)
    new_data3 <- new_data3[-1, ]
    
    zenbatdr <- n2 - 1
    
    laburpenite <- NA
    
    quar <- paste0("newdata_it")
    assign(quar, data.frame(matrix(NA, 1, len_tot)))
    names(newdata_it) <- names(dataset[, var_tot])
    
    z <- len_tot
    n_inter <- z * (z - 1) / 2
    len_4_int <- 4 * n_inter
    
    for (r in 1:len_4_int) {
        tarteko <- get("newdata_it")
        names(tarteko) <- names(dataset[, var_tot])
        assign(quar, rbind(tarteko, dataset[, var_tot]))
    }
    
    assign(quar, get(quar)[-1, ])
    
    ec <- 1
    h <- 0
    n <- 1
    
    for (ex in 1:len_tot) {
        
        ec <- ec + 1
        if (len_tot == ex) 
            break
        for (ex2 in ec:len_tot) {
            h <- h + 1
            n2 <- n + N
            newdata_it[n:(n2 - 1), ex] <- delta[1]
            newdata_it[n:(n2 - 1), ex2] <- delta[1]
            laburpenite <- c(laburpenite, paste0("ITE_", var_tot[ex]), 
                             var_tot[ex2], 0, n, (n2 - 1))
            n <- n2
            n2 <- n + N
            h <- h + 1
            newdata_it[n:(n2 - 1), ex] <- delta[1]
            newdata_it[n:(n2 - 1), ex2] <- delta[2]
            laburpenite <- c(laburpenite, paste0("ITE_", var_tot[ex]), 
                             var_tot[ex2], 1, n, (n2 - 1))
            n <- n2
            n2 <- n + N
            h <- h + 1
            newdata_it[n:(n2 - 1), ex] <- delta[2]
            newdata_it[n:(n2 - 1), ex2] <- delta[1]
            laburpenite <- c(laburpenite, paste0("ITE_", var_tot[ex]), 
                             var_tot[ex2], 10, n, (n2 - 1))
            n <- n2
            n2 <- n + N
            h <- h + 1
            newdata_it[n:(n2 - 1), ex] <- delta[2]
            newdata_it[n:(n2 - 1), ex2] <- delta[2]
            laburpenite <- c(laburpenite, paste0("ITE_", var_tot[ex]), 
                             var_tot[ex2], 11, n, (n2 - 1))
            n <- n2
            
        }
    }
    
    laburpenite <- laburpenite[-1]
    laburpenite2 <- data.frame(matrix(laburpenite, len_4_int, 5, byrow = TRUE))
    
    new_tot <- rbind(newdata, new_data3, newdata_it)
    
    
    names(laburpencc2) <- names(laburpendr2) <- names(laburpenite2) <- c("Group", "Case", "SubGroup", "From", "To")
    
    laburtot <- rbind(laburpencc2, laburpendr2, laburpenite2)
    
    laburtot$From <- as.numeric(as.character(laburtot$From))
    laburtot$To <- as.numeric(as.character(laburtot$To))
    
    lastlab <- dim(laburtot)[1]
    len3 <- len_delta * len_exp + 1
    len4 <- len3 + len_sek_exp
    
    laburtot$From[c(len3:lastlab)] <- laburtot$From[len3:lastlab] + zenbatcc
    laburtot$To[c(len3:lastlab)] <- laburtot$To[len3:lastlab] + zenbatcc
    laburtot$From[(len4:lastlab)] <- laburtot$From[len4:lastlab] + (zenbatdr * len_exp)
    laburtot$To[(len4:lastlab)] <- laburtot$To[len4:lastlab] + (zenbatdr * len_exp)
    
    
    return(list(new_tot, laburtot))
}

#' This function run all the simulations based on the structure of the basic data frame. 
#'
#' @param dataset dataset with all variables
#' @param exposures a vector with exposures
#' @param confounders a vector with confounders
#' @param libraries a vector of libraries to use
#' @param outcomes outcome´s name
#' @param num.sim number of simulations
#' @param delta a vector with two values
#' @param dr a vector with dose response values
#' @param newdata the dataframe with new values
#' @param show_times boolean parameter to see the time
#' @param show_num_sim boolean parameter to see the iteration of simulations
#' @param save_time boolean parameter to save the time in the result list
#' @param verbose boolean parameter to see the verbose of superlearner
#' @param family a character parameter to describe the family of the model
#' @param method a character parameter to choose the method in the superlearner
#' @return a list with 4 objects: a data frame with all simulations, risk and coefficients of the crossvalidation and the time of the proces.
#' @export
#' @details libraries could be SL if we don't select nothing or 'SL.glm', 'SL.glm.interaction','SL.glmnet', 'SL.gam','SL.xgboost','SL.polymars','SL.randomForest'


run_simulations <- function(dataset, exposures, confounders, libraries, outcomes, num.sim = 50, delta = c(0, 1), dr, 
    newdata, show_times = FALSE, show_num_sim = TRUE, save_time = FALSE, verbose = FALSE, family = "gaussian", method = "method.NNLS") {
    if (class(exposures) %in% c("numeric", "integer")) {
        exposures <- names(dataset[, exposures])
    }
    if (class(confounders) %in% c("numeric", "integer")) {
        confounders <- names(dataset[, confounders])
    }
    if (class(outcomes) %in% c("numeric", "integer")) {
        outcomes <- names(dataset[, outcomes])
    }
    
    
    len_libraries <- length(libraries)
    if (len_libraries == 0) {
        simple_libr <- c("")
    } else {
        simple_libr <- gsub("SL.", "", libraries)
        simple_libr <- substr(gsub("\\.", "", simple_libr), 1, 5)
        simple_libr <- c("", simple_libr)
    }
    len_libraries <- length(libraries)
    var_tot <- c(exposures, confounders)
    
    d <- dataset
    N <- dim(d)[1]
    len_new <- dim(newdata)[1]
    
    big_matrix <- data.frame(matrix(NA, len_new, num.sim))
    
    risk <- data.frame(libraries, matrix(NA, len_libraries, num.sim))
    coef <- data.frame(libraries, matrix(NA, len_libraries, num.sim))
    
    if (save_time == TRUE) 
        time3 <- NA
    t1 <- Sys.time()
    
    if (show_num_sim == TRUE) {
        for (i in 1:num.sim) {
            print(paste0("Run number ", i, "/", num.sim))
            
            idx <- sample(1:N, N, replace = T)
            dataset <- d[idx, ]
            
            modelo <- SuperLearner::SuperLearner(Y = dataset[, outcomes], X = dataset[, var_tot], SL.library = libraries, 
                newX = newdata[, var_tot], verbose = verbose, family = family, method = method)
            risk[, i + 1] <- modelo$cvRisk
            coef[, i + 1] <- modelo$coef
            big_matrix[, i] <- modelo$SL.predict[, 1]
        }
        
    } else {
        print(paste0("Start running - ", Sys.time()))
        for (i in 1:num.sim) {
            idx <- sample(1:N, N, replace = T)
            dataset <- d[idx, ]
            
            modelo <- SuperLearner::SuperLearner(Y = dataset[, outcomes], X = dataset[, var_tot], SL.library = libraries, 
                newX = newdata[, var_tot], verbose = verbose, family = family, method = method)
            risk[, i + 1] <- modelo$cvRisk
            coef[, i + 1] <- modelo$coef
            big_matrix[, i] <- modelo$SL.predict[, 1]
        }
        print(paste0("Finished - ", Sys.time()))
    }
    
    if (show_times == TRUE) {
        t2 <- Sys.time()
        print(t2 - t1)
        if (save_time == TRUE) 
            time3 <- c(time3, t2 - t1)
    }
    
    fecha <- format(Sys.time(), "%Y%m%d%X")
    fecha <- gsub(":", "", fecha)
    return(list(big_matrix, risk, coef, time_proc = time3))
    
}

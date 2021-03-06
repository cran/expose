% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interact.r
\name{interact}
\alias{interact}
\title{Extract the information from the simulation data frame to analyse the interaction effects}
\usage{
interact(allsim, dataset, exposures, confounders, squem)
}
\arguments{
\item{allsim}{dataset with all simulations values}

\item{dataset}{dataset with all variables}

\item{exposures}{a vector with exposures}

\item{confounders}{a vector with confounders}

\item{squem}{squeme of the values of the prediction values}
}
\value{
data frame with interaction values
}
\description{
Extract the information from the simulation data frame to analyse the interaction effects
}
\examples{
data(expose_data)
data(simu)
data(gen)
delta=c(1,0)
seku <- seq(0,1,0.05)
Exposures<- c('Var1','Var2','Var3','Var4','Var5')
summary_table_lines <- gen[[2]]
it <- interact (allsim = simu[[1]], dataset = expose_data,exposures = Exposures,
confounders = c('sex'), squem = summary_table_lines)
}

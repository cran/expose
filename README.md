# expose
 Introduction: 
 
Although the general population is exposed to multiple chemicals, most environmental epidemiology studies consider each chemical separately when estimating the adverse effects of environmental exposures, partly because of the lack of software accessible to environmental epidemiologists. 

Objective: 

We developed a statistical R package, accessible to environmental health scientists that can provide valid estimates for effects of environmental contaminants, dose-response relationships, and interactions in a multi-pollutant setting. 

Methods: 

This package implements the work of Oulhote et al (2017) that combines the G-formula, a maximum likelihood estimator, with the ensemble learning technique Super Learner. It incorporates a wide range of statistical techniques (e.g. generalized linear and additive models, random forest, extreme gradient boosting) to estimate and provide valid inference in multi-pollutant settings and reconstruct dose-response relationships non-parametrically. 

Results: 

This package will facilitate multi-pollutant analysis of environmental epidemiology. We also developed a user-friendly Shiny application that can be used by biomedical researchers modeling complex mixtures. We ran multiple simulations based on real scenarios and the proposed method yielded promising results across all the scenarios. Our approach was also able to estimate the true underlying structure of the data. We will present a tutorial describing the package functions and visualizations to illustrate our developed methods. 

Conclusion: 

This package will help to unravel the effects of chemical mixtures and their interactions in epidemiological studies. 


Youssef Oulhote, Marie-Abele Bind, Brent Coull, Chirag, Patel, Philippe Grandjean. 2017. Combining Ensemble Learning Techniques and G-Computation to Investigate Chemical Mixtures in Environmental Epidemiology Studies. Biorxiv. https://www.biorxiv.org/content/early/2017/06/30/147413.article-info


## Installation

You can install the released version of expose from www.github.com with:

``` r
library(devtools)
 devtools::install_github("itamuria/expose")
```

## Case study

This is a basic case study similar to the analysed on the paper. 

``` r
rm(list = ls())
time_doc1 <- Sys.time()


library(expose)

if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(ggplot2, repmis,RColorBrewer,SuperLearner, gam,splines,foreach,glmnet,Matrix,nnet,polspline,e1071,xgboost)


# Download data
source_data("https://github.com/itamuria/expose_dataset/blob/master/20180818_dtaset.RData?raw=true")
set.seed(22222)

# Define parameters

N <- dim(dtaset)[1]
Outcome="Y4"
seku <- seq(0,1,0.05)   #c(0,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,1)
our.num.sim <- 5
delta=c(1,0)


Exposures<- c("Var1","Var2","Var3","Var4","Var5")
Confounders<- c("sex")


###############################################################################################################
# Crossvalidation
###############################################################################################################

ourlibraries <- c("SL.glm", "SL.glm.interaction","SL.glmnet", "SL.gam","SL.nnet",
                  "SL.polymars","SL.svm","SL.xgboost")




Outcome 4
Outcome <- c("Y4")
m4 <- SuperLearner::CV.SuperLearner(Y = dtaset[,Outcome[1]],
                                   X = dtaset[, c(Confounders,Exposures)],V=10, SL.library = ourlibraries,
                                   family = "gaussian", method = "method.NNLS",
                                   verbose = FALSE)

mm <- m4$AllSL[[1]]

###############################################################################################################
# create basic data frame
###############################################################################################################

Outcome <- c("Y4")

gen <- general_function (dataset = dtaset, exposures = Exposures,
                         confounders = Confounders,
                         outcomes = Outcome[1], delta=delta, dr = seku)


summary_table_lines <- gen[[2]]



###############################################################################################################
# Simulations
###############################################################################################################

t1 <- Sys.time()
simu <- run_simulations (dataset = dtaset, exposures = Exposures,
                          confounders = Confounders, libraries = ourlibraries,
                          outcomes = Outcome[1], num.sim = our.num.sim, delta=delta, dr = seku,
                          newdata =gen[[1]], show_times = TRUE, verbose = FALSE, save_time = TRUE,
                         show_num_sim = TRUE, family = "gaussian", method ="method.NNLS")

t2 <- Sys.time()
(t2-t1)

h <- 4

###############################################################################################################
# Sensibility analysis
###############################################################################################################


################################# risk

ris <- simu[[2]]

len_sim <- dim(ris)[2]

lib <- gsub("SL.","",ris$libraries)
len_lib <- length(lib)

df1 <- data.frame(matrix(NA,1,2))
names(df1) <- c("Library","Value")

for(li in 1:len_lib)
{
  df2 <- data.frame(rep(lib[li],len_sim-1),t(ris[li,2:len_sim]))
  names(df2) <- c("Library","Value")
  df1 <- rbind(df1,df2)
}
df1 <- df1[-1,]

ggplot(df1, aes(x = Library, y = Value, fill = Library)) + geom_boxplot()  + theme_bw() +  
  coord_flip() + ggtitle("Risk per method (minimize)") + theme(legend.position="none")

# name to save
file_exp_name <- paste0("Y",h,"_Risk.jpg")

ggsave(file_exp_name)
#unlink(file_exp_name)


ris_table <- data.frame(lib,apply(ris[,-1],1,mean))
names(ris_table)<-c("Library","Performance")
knitr::kable(ris_table)


################################# Coefficient

coef <- simu[[3]]

coef_table <- data.frame(lib,apply(coef[,-1],1,mean))
names(coef_table)<-c("Library","Estimates")
knitr::kable(coef_table)

df1 <- data.frame(matrix(NA,1,2))
names(df1) <- c("Library","Value")

for(li in 1:len_lib)
{
  df2 <- data.frame(rep(lib[li],len_sim-1),t(coef[li,2:len_sim]))
  names(df2) <- c("Library","Value")
  df1 <- rbind(df1,df2)
}
df1 <- df1[-1,]

ggplot(df1, aes(x = Library, y = Value, fill = Library)) + geom_boxplot() + theme_bw() +
  coord_flip() + ggtitle("Coefficient per method (weight)") + theme(legend.position="none") 

file_exp_name <- paste0("Y",h,"_Coefficient.jpg")

ggsave(file_exp_name)


###############################################################################################################
# Naive ACE: Average causal effect
###############################################################################################################


ace.df.g <- naive_ace (allsim = simu[[1]], dataset = dtaset, ic_dis = "IC", st = summary_table_lines,
                       exposures = Exposures, delta = delta)

knitr::kable(ace.df.g)

# Use 95% confidence interval instead of SEM
pd <- position_dodge(0.1)
ggplot(ace.df.g, aes(x=Group, y=Mean, colour=Group)) + 
  geom_errorbar(aes(ymin=ICa, ymax=ICb), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) + ylab("Estimated (Boostrapped 95% CI") +
  xlab("Library") + geom_hline(yintercept = 0)+ theme_bw() +
  ggtitle("Average Treatment Effect") + theme(legend.position="none")


file_exp_name <- paste0("Y",h,"_NaiveAce.jpg")

ggsave(file_exp_name)

###############################################################################################################
# Dose Respond
###############################################################################################################


drr.grp <- dose_resp (allsim = simu[[1]], dataset = dtaset, st = summary_table_lines,
                     dr = seku, exposures = Exposures) 

knitr::kable(head(drr.grp))
drr.grp$dose <- as.numeric(gsub("DR_","",drr.grp$Quantile))

pd <- position_dodge(0.1)

ggplot(drr.grp, aes(x=dose, y=Mean)) + 
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.025, position=pd, size=0.5, color="blue") + geom_hline(yintercept = 0) + 
  geom_line(position=pd,col="blue") +
  geom_point(position=pd, size=2, shape=20, fill="black") + facet_grid(Exp ~ ., scales="free") +  
  ggtitle("Dose response") + theme(legend.position="none") + theme_bw()

file_exp_name <- paste0("Y",h,"_DoseResponse1.jpg")

ggsave(file_exp_name)

ggplot(drr.grp, aes(x=dose, y=Mean)) + 
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.05, position=pd, size=0.5, color="blue") +
  geom_line(position=pd,col="blue") +
  geom_point(position=pd, size=2, shape=20, fill="black") + facet_grid(Exp ~ .) + geom_hline(yintercept = 0) +  
  ggtitle("Dose response") + theme(legend.position="none") + theme_bw()

file_exp_name <- paste0("Y",h,"_DoseResponse2.jpg")

ggsave(file_exp_name)


###############################################################################################################
# ice: individual causal effect
###############################################################################################################

ice_res <- ice(allsim = simu[[1]], dataset = dtaset, dr = seku, squem = summary_table_lines, remove_extrem = FALSE)

# ploting
p=ggplot(ice_res,aes(X, Var1_pred,group=id))
p+geom_line()+theme_bw() + 
  scale_x_discrete(breaks=sort(unique(ice_res$X)),
                   labels=as.character(seku))

# name to save
file_exp_name <- paste0("Y",h,"_ICE_1.jpg")

ggsave(file_exp_name)

# coloring 
p=ggplot(ice_res,aes(X, Var2_pred,group=id, color=Var3)) 
p+geom_line()+ theme_bw() + theme(legend.position="none") + ggtitle("Individual effect") + 
  scale_x_discrete(breaks=sort(unique(ice_res$X)),
                   labels=as.character(seku))

file_exp_name <- paste0("Y",h,"_ICE_2.jpg")

ggsave(file_exp_name)

p=ggplot(ice_res,aes(X, Var1_pred,group=id,col=Var3)) + labs(x="Var1 percentiles",y="Predicted Y")
p+geom_line()+theme_bw()+scale_color_gradientn(name="Var3",colours=rev(brewer.pal(9,"YlOrRd"))) + 
  scale_x_discrete(breaks=sort(unique(ice_res$X)),
                   labels=as.character(seku)) + ggtitle("Individual effect")

file_exp_name <- paste0("Y",h,"_ICE_3.jpg")

ggsave(file_exp_name)


p=ggplot(ice_res,aes(X, Var3_pred,group=id,col=Var5)) + labs(x="Var3 percentiles",y="Predicted Y")
p+geom_line()+theme_bw()+scale_color_gradientn(name="Var5",colours=rev(brewer.pal(9,"YlOrRd"))) +
  scale_x_discrete(breaks=sort(unique(ice_res$X)),
                   labels=as.character(seku)) + ggtitle("Individual effect")

file_exp_name <- paste0("Y",h,"_ICE_4.jpg")

ggsave(file_exp_name)

# discrete variables

p=ggplot(ice_res,aes(X, Var3_pred,group=id,col=sex)) + labs(x="pp-Var3 percentiles",y="Predicted Y")
p+geom_line()+theme_bw()+scale_color_gradientn(name="Sex",colours=c("red","#FFFF00")) +
  theme(axis.text.x = element_text(face="bold",size=8, angle=0), legend.position = "none") + 
  scale_x_discrete(breaks=sort(unique(ice_res$X)),
                   labels=as.character(seku)) + ggtitle("Individual effect")

file_exp_name <- paste0("Y",h,"_ICE_5.jpg")

ggsave(file_exp_name)

###############################################################################################################
# Interactions
###############################################################################################################


it <- interact (allsim = simu[[1]], dataset = dtaset,exposures = Exposures,
                      confounders = c("sex"), squem = summary_table_lines)

knitr::kable(head(it))

pd <- position_dodge(0.1)
ggplot(it, aes(x=Interaction, y=Mean, colour=Interaction)) + 
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) + ylab("Estimated (Boostrapped 95% CI") +
  xlab("Interaction") + coord_flip() + geom_hline(yintercept = 0) + 
  ggtitle("Interaction") + theme_bw() + theme(legend.position="none") 

file_exp_name <- paste0("Y",h,"_Inter.jpg")

ggsave(file_exp_name)


time_doc3 <- Sys.time()

print(time_doc3 - time_doc1)

```


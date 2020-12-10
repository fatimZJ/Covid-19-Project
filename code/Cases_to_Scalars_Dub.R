## Load packages
library(deSolve)
library(DEoptim)

## Load the data
#source('code/1_loadData.r')

## Load the get beta function

source("code/function_IEMAGSEIR_Dublin.R")


nlm_fun <- function(scalars, Actual_Cc, population_df, contacts_list, interventions_df) {
  
 scalars_est <- exp(scalars)  
  
  names(scalars_est) <- unique(interventions_df$policy)
  Base <- simulation_SEIR_model(R0t = 3.4,
                                POP = population_df,
                                contacts_ireland = contacts_list,
                                interventions = interventions_df, 
                                dt = 1,
                                dateStart = interventions_df$start[1],
                                dateEnd = cumulative_cases$date[dim(cumulative_cases)[1]],# Time step (days)
                                scalars = scalars_est)
  
  Est_Cc <- rowSums(Base$sol_out[grepl('Cc_',names(Base$sol_out))])
  
  sum((Actual_Cc - (Est_Cc[-1]))^2) 
  
}


scalars_init <- c(1.362216848, 0.510116460, 0.106050230, 0.008395608, 0.122922232,
                  0.197667229, 0.248520907, 0.171967661, 0.091010276)
                #(0.971304374, 0.553483970, 0.101622333,
                #  0.009065504, 0.120585875, 0.258273993,
                #  0.225591140, 0.167105331, 0.220136291)


nlm_fun(scalars = c(1.000000, 0.97171108, 0.20812880, 0.07611468, 0.22169280, 
                    0.40394184, 0.44685362, 0.36467352, 0.42729580), 
        Actual_Cc = cumulative_cases$cases, interventions_df = interventions_info, population_df = population,
        contacts_list = contacts)


plot(cumulative_cases$cases, type = "l")

startt <- Sys.time()
ests_nlm <- nlm(nlm_fun,log(scalars_init),
                stepmax = 0.5,interventions_df = interventions_info, population_df = population,
                contacts_list = contacts, 
                iterlim = 1000, Actual_Cc = cumulative_cases$cases)
endt <- Sys.time()


scalars <- log(scalars_init)


################################################################################
exp(ests_nlm$estimate)
scalars_nlm <- c(1.362216848, 0.510116460, 0.106050230, 0.008395608, 0.122922232,
                  0.197667229, 0.248520907, 0.171967661, 0.091010276)
## with rounding t
c(6.858719e-02, 1.839589e+00, 9.594844e-02, 4.283682e-04, 2.721637e-06,
 3.610250e-01, 2.099779e-01, 1.736524e-01, 9.436696e-02)
##################################################################################
lower <- rep(0, length(scalars_init))
upper <- rep(3, length(scalars_init))

ests_DEoptim <- DEoptim(nlm_fun, Actual_Cc = cumulative_cases$cases,
                        lower, upper, DEoptim.control(itermax = 200, reltol = 1e-8,
                                                      steptol = 50))

startt <- Sys.time()
ests_optimx <-  optimx(c(1.000000, 0.97171108, 0.20812880, 0.07611468, 0.22169280, 
                         0.40394184, 0.44685362, 0.36467352, 0.42729580), 
                       Actual_Cc = cumulative_cases$cases, itnmax = 200,
                       nlm_fun, 
                       method = c("Nelder-Mead"),
                       control = list(maximize = FALSE))
endt <- Sys.time()


summaryRprof(filename="Profile.out")



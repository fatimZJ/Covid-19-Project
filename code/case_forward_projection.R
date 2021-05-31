##### Forward Projections (from May 17th)

library(readr)
library(zoo)
library(Matrix)
library(deSolve)
setwd('Joint_Shiny_App/')
source('code/02_get_beta.R')
source('code/04-function_SEIR_simulation.R')
load('data/contacts_IRL.Rdata')
linfo <- read_csv('data/Dublin_Interventions_forward_projection.csv')
linfo$start <- as.Date(linfo$start, format = '%d/%m/%Y')
linfo$end <- as.Date(linfo$end, format = '%d/%m/%Y')

pop <- read_csv('data/Dublin_pop_2019.csv')

startDate <- as.Date('2020-02-29')
endDate_tmax <- difftime(as.Date('2021-01-31'), startDate)[[1]]

def_pars <- c(3.7, 5.8, 11.6, 0.55, 0.2, 0.8, 0.1, 0.05, 7)
names(def_pars) <- c('L', 'Cv', 'Dv', 'h', 'f', 'tv', 'q', 'k', 'TT')
dub_def_beta <- getbeta(3.7, pars = def_pars, p_age = pop$propage,
            CONTACTMATRIX = contacts_IRL)

dub_xstart <- c(pop$popage - 15/16,
                  rep(1, 16),
                  rep(1/16, 16),
                  rep(0, 16 * 6))
names(dub_xstart) <- c(paste0('S_',1:16),
                       paste0('Ev_',1:16),
                       paste0('Ip_',1:16),
                       paste0('IA_',1:16),
                       paste0('Ii_',1:16),
                       paste0('It_',1:16),
                       paste0('Iti_',1:16),
                       paste0('Iq_',1:16),
                       paste0('R_',1:16))

# Run the SEIR model and extract susceptibles
tic()
res <- SEIR_model_simulation(pars = def_pars,
                             dateStart = startDate,
                             lockdown_information = linfo[c(1:2, 4)],
                             POP = pop,
                             contacts_ireland = contacts_IRL,
                             beta = dub_def_beta,
                             startval = dub_xstart,
                             dt = 1,  
                             tmax = endDate_tmax,
                             isolated_contacts = contacts_IRL)
toc()
res_lb <- SEIR_model_simulation(pars = def_pars,
                             dateStart = startDate,
                             lockdown_information = linfo[c(1:2, 6)],
                             POP = pop,
                             contacts_ireland = contacts_IRL,
                             beta = dub_def_beta,
                             startval = dub_xstart,
                             dt = 1,  
                             tmax = endDate_tmax,
                             isolated_contacts = contacts_IRL)
res_ub <- SEIR_model_simulation(pars = def_pars,
                             dateStart = startDate,
                             lockdown_information = linfo[c(1:2, 5)],
                             POP = pop,
                             contacts_ireland = contacts_IRL,
                             beta = dub_def_beta,
                             startval = dub_xstart,
                             dt = 1,  
                             tmax = endDate_tmax,
                             isolated_contacts = contacts_IRL)

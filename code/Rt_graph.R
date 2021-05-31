##### Effective Reproduction Rate Graph
library(readr)
library(zoo)
library(Matrix)
library(deSolve)
setwd('Joint_Shiny_App/')
source('code/02_get_beta.R')
source('code/04-function_SEIR_simulation.R')
load('data/contacts_IRL.Rdata')
linfo <- read_csv('data/Dublin_Interventions.csv')
linfo$start <- as.Date(linfo$start, format = '%d/%m/%Y')
linfo$end <- as.Date(linfo$end, format = '%d/%m/%Y')
linfo[[4]] <- c( 1.01952, 2.114842, 0.1554196, 0.05411131, 0.1555254, 
                 0.42444, 0.459333, 0.3194528, 0.1520108, 0.8319415,
                 0.7270618, 0.0173905 )
linfo[[5]] <- c( 0.748,
1.681,
0.131,
0.026,
0.098,
0.349,
0.396,
0.273,
0.105,
0.586,
0.542,
0.003 )
linfo[[6]] <- c( 1.471,
2.424,
0.179,
0.112,
0.259,
0.478,
0.556,
0.373,
0.192,
1.135,
0.962,
0.110 )

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
S_all <- res$solution[2:17]
S <- rowSums(S_all)
S_all_ub <- res_ub$solution[2:17]
S_ub <- rowSums(S_all_ub)
S_all_lb <- res_lb$solution[2:17]
S_lb <- rowSums(S_all_lb)

prop <- S / sum(pop$popage)
prop_lb <- S_lb / sum(pop$popage)
prop_ub <- S_ub / sum(pop$popage)

# Compute effective reproductive numbers
R0 <- 3.7
date_list <- vector("list", nrow(linfo))
for (i in seq_along(date_list)) {
  
  date_list[[i]] <- seq.Date(linfo[[1]][i], linfo[[2]][i], 1)
  
}

lens <- sapply(date_list, length)
theta <- rep(linfo[[4]], lens)
theta_lb <- rep(linfo[[5]], lens)
theta_ub <- rep(linfo[[6]], lens)

solid_line <- theta*R0*prop
upper <- theta_ub*R0*prop_ub
lower <- theta_lb*R0*prop_lb
xaax <- as.Date( unlist(date_list) )

pdf('Rt_plot_corrected.pdf', width = 10)
plot(xaax, rep(1,length(xaax)), type="n", xlab="Time (beginning 2020)", ylab="Effective Reproductive Number", 
     ylim=c(0, max(upper)), main = "Effective Reproductive Number over Time")
polygon(c(rev(xaax), xaax), c(rev(upper), lower), col = 'grey80', border = NA)
abline(v = xaax[c(2, 63, 124, 186, 247, 308)], lty = 3, col = 'gray')
abline(h = seq(0, 8, 2), lty = 3, col = 'gray')
lines(xaax, solid_line, type = 'l', col = 'dodgerblue', lwd = 2)
abline(h = 1, lty = 2, col = 'orange', lwd = 2)
dev.off()

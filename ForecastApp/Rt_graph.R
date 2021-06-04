##### Effective Reproduction Rate Graph
library(readr)
library(zoo)
library(Matrix)
library(deSolve)
setwd('ForecastApp/')
source('code/02_get_beta.R')
source('code/04-function_SEIR_simulation.R')
load('data/contacts_IRL.Rdata')
ol <- read_csv('data/ol_table_S2.csv')
linfo <- read_csv('data/Dublin_Interventions.csv')
linfo$start <- as.Date(linfo$start, format = '%d/%m/%Y')
linfo$end <- as.Date(linfo$end, format = '%d/%m/%Y')
linfo <- cbind(linfo, ol[c(2:3, 4:5)])

pop <- read_csv('data/Dublin_pop_2019.csv')

startDate <- as.Date('2020-02-29')
endDate_tmax <- difftime(as.Date('2021-01-31'), startDate)[[1]]

def_pars <- c(3.8, 5.8, 13.5, 0.55, 0.2, 0.8, 0.1, 0.05, 7)
names(def_pars) <- c('L', 'Cv', 'Dv', 'h', 'f', 'tv', 'q', 'k', 'TT')
dub_def_beta <- getbeta(3.4, pars = def_pars, p_age = pop$propage,
            CONTACTMATRIX = contacts_IRL)

num_exp <- 14.5344/16
num_inf <- 0.947286/16 

start_non_S <- c(rep(num_exp, 16), rep(num_inf, 16), rep(0, 16 * 7))

dub_xstart <- c(pop$popage - num_exp - num_inf, start_non_S)
names(dub_xstart) <- c(paste0('S_',1:16),
                       paste0('Ev_',1:16),
                       paste0('Ip_',1:16),
                       paste0('IA_',1:16),
                       paste0('Ii_',1:16),
                       paste0('It_',1:16),
                       paste0('Iti_',1:16),
                       paste0('Iq_',1:16),
                       paste0('R_',1:16),
                       paste0('Cc_',1:16))

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
                             lockdown_information = linfo[c(1:2, 7)],
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
R0 <- 3.4
date_list <- vector("list", nrow(linfo))
for (i in seq_along(date_list)) {
  
  date_list[[i]] <- seq.Date(linfo[[1]][i], linfo[[2]][i], 1)
  
}

lens <- sapply(date_list, length)
theta <- rep(linfo[[4]], lens)
theta_lb <- rep(linfo[[6]], lens)
theta_ub <- rep(linfo[[7]], lens)

solid_line <- theta*R0*prop
upper <- theta_ub*R0*prop_ub
lower <- theta_lb*R0*prop_lb
xaax <- as.Date( unlist(date_list) )

pdf('Rt_plot_corrected.pdf', width = 10)
plot(xaax, rep(1,length(xaax)), type="n", xlab="Date (beginning 2020)", ylab="Effective Reproductive Number", 
     ylim=c(0, 9), main = "Effective Reproductive Number over Time")
#polygon(c(rev(xaax), xaax), c(rev(upper), lower), col = 'grey80', border = NA)
abline(v = xaax[c(2, 63, 124, 186, 247, 308)], lty = 3, col = 'gray')
abline(h = seq(0, 8, 2), lty = 3, col = 'gray')
lines(xaax, solid_line, type = 'l', col = 'dodgerblue', lwd = 2)
lines(xaax, upper, type = 'l', col = 'green', lwd = 2)
lines(xaax, lower, type = 'l', col = 'yellow', lwd = 2)
legend('topright', col=c("dodgerblue", "green", "yellow", "orange"), lty=c(rep(1, 3), 2),
       legend=c("Optimised Scalers", "Upper Scalars", "Lower Scalars", "Epidemic Threshold"),
       bty = 'n')
abline(h = 1, lty = 2, col = 'orange', lwd = 2)
dev.off()

###

Cc <- rowSums(res$solution[146:161])
pdf('Cc_plot.pdf')
plot( Cc, col = 'dodgerblue')#, lwd = 2 )
dev.off()
max(Cc)

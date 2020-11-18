## Adapted from https://gist.github.com/jonocarroll/b17ce021b0637a31f584ed08a1fbe733
read.tcsv = function(file, header=TRUE, sep=",", ...) {
  
  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)[-1]
  
  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }
  
  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep)
  ## empty strings are converted to NA
  out = read.csv(text=x, sep=sep, header=header, na.strings = "", ...)
  return(out)
  
}

jg_dat <- read.tcsv("data/dat_seir_code.csv")
jg_dat$Date <- as.Date(jg_dat$Date, format = "%a %d %b %Y") # Reformat date column


tmax <- 225
dt <- 1
times <- seq(1, tmax, by = dt)

dateStartSchoolClosure <- as.Date('2020-03-12') # School closure
dateStartIntenseIntervention <- as.Date('2020-03-27')  #Intense intervention
dateEndIntenseIntervention <- as.Date('2020-05-18') #date we begin relaxing intense intervention 
dateStart <- as.Date('2020-02-28') #start date for epidemic in Ireland

numWeekStagger <- c(3,6,9,12,15)

tStartSchoolClosure <- (as.vector(dateStartSchoolClosure - dateStart)+1) #Time point to add the school closure effect
tStartIntenseIntervention <- as.vector(dateStartIntenseIntervention - dateStart)+1 # #Time point to add the intense lockdown effect
tEndIntenseIntervention <- as.vector(dateEndIntenseIntervention - dateStart)+1     # for pw = 0.1
tRelaxIntervention1 <- tEndIntenseIntervention + numWeekStagger[1]*7                               # for pw = 0.25
tRelaxIntervention2 <- tEndIntenseIntervention + numWeekStagger[2]*7                               # for pw = 0.5
tRelaxIntervention3 <- tEndIntenseIntervention + numWeekStagger[3]*7
tRelaxIntervention4 <- tEndIntenseIntervention + numWeekStagger[4]*7                               # for pw = 0.25
tRelaxIntervention5 <- tEndIntenseIntervention + numWeekStagger[5]*7                               # for pw = 0.5

## create step function from james gleeson's betas

## means

Initial_beta <- mean(jg_dat$Beta[1:tStartSchoolClosure - 1])

schoolClosure_beta <- mean(jg_dat$Beta[tStartSchoolClosure:tStartIntenseIntervention - 1])

lockdown_beta <- mean(jg_dat$Beta[tStartIntenseIntervention:tEndIntenseIntervention - 1])

RelaxInter1_beta <- mean(jg_dat$Beta[tEndIntenseIntervention:tRelaxIntervention1 - 1])

RelaxInter2_beta <- mean(jg_dat$Beta[tRelaxIntervention1:tRelaxIntervention2 - 1])

RelaxInter3_beta <- mean(jg_dat$Beta[tRelaxIntervention2:tRelaxIntervention3 - 1])

RelaxInter4_beta <- mean(jg_dat$Beta[tRelaxIntervention3:tRelaxIntervention4 - 1])

RelaxInter5_beta <- mean(jg_dat$Beta[tRelaxIntervention4:tRelaxIntervention5 - 1])

postlockdown_beta <- mean(jg_dat$Beta[tRelaxIntervention5:tmax])

## medians 
# 
# Initial_beta <- median(jg_dat$Beta[1:tStartSchoolClosure - 1])
# 
# schoolClosure_beta <- median(jg_dat$Beta[tStartSchoolClosure:tStartIntenseIntervention - 1])
# 
# lockdown_beta <- median(jg_dat$Beta[tStartIntenseIntervention:tEndIntenseIntervention - 1])
# 
# RelaxInter1_beta <- median(jg_dat$Beta[tEndIntenseIntervention:tRelaxIntervention1 - 1])
# 
# RelaxInter2_beta <- median(jg_dat$Beta[tRelaxIntervention1:tRelaxIntervention2 - 1])
# 
# RelaxInter3_beta <- median(jg_dat$Beta[tRelaxIntervention2:tRelaxIntervention3 - 1])
# 
# RelaxInter4_beta <- median(jg_dat$Beta[tRelaxIntervention3:tRelaxIntervention4 - 1])
# 
# RelaxInter5_beta <- median(jg_dat$Beta[tRelaxIntervention4:tRelaxIntervention5 - 1])
# 
# postlockdown_beta <- median(jg_dat$Beta[tRelaxIntervention5:tmax])


step_beta <- c(rep(Initial_beta, length(1:tStartSchoolClosure - 1)), 
                   rep(schoolClosure_beta, length(tStartSchoolClosure:tStartIntenseIntervention - 1)),
               rep(lockdown_beta, length(tStartIntenseIntervention:tEndIntenseIntervention - 1)), 
               rep(RelaxInter1_beta, length(tEndIntenseIntervention:tRelaxIntervention1 - 1)),
               rep(RelaxInter2_beta, length(tRelaxIntervention1:tRelaxIntervention2 - 1)), 
               rep(RelaxInter3_beta, length(tRelaxIntervention2:tRelaxIntervention3 - 1)),
               rep(RelaxInter4_beta, length(tRelaxIntervention3:tRelaxIntervention4 - 1)), 
               rep(RelaxInter5_beta, length(tRelaxIntervention4:tRelaxIntervention5 - 1)),
               rep(postlockdown_beta, length(tRelaxIntervention5:tmax)))

plot(jg_dat$Beta[1:tmax], type = "l", lwd = 2,
     panel.first = rect(c(tStartSchoolClosure, tStartIntenseIntervention, tEndIntenseIntervention,
                          tRelaxIntervention1, tRelaxIntervention2, tRelaxIntervention3,
                          tRelaxIntervention4), -1e6,
                        c(tStartIntenseIntervention,tEndIntenseIntervention, tRelaxIntervention1,
                          tRelaxIntervention2,tRelaxIntervention3,
                          tRelaxIntervention4,tRelaxIntervention5 ), 1e6,
                        col=c('gray63','gray48', 'gray70', 'gray75','gray80','gray85', 'gray90'), border=NA))
lines(step_beta, lwd  = 2, col = "red")

################################################################################
## Load the data
source('code/1_loadData.r')
contacts_ireland <- contacts

## Load the get beta function 
source("code/getbeta.R")

## Model Parameters

pars <- c(4.9, 5.9, 7.0, 0.25, 0.05, 0.05, 0.5, 0.75, 0.13, 3.6)
names(pars) <- c("L","Cv","Dv","h","i","j","f","tv","q","TT")

R <- vector()
Beta <- vector()

for (t in 1: tmax) {
R[t] <- getR0t(beta = step_beta[t], pars = pars, p_age = Irlpop$propage, CONTACTMATRIX = contacts_ireland)

Beta[t] <- getbeta(R0t = R[t], pars = pars, p_age = Irlpop$propage, CONTACTMATRIX = contacts_ireland)

}

R
Beta

################################################################################

## To get the our equivalent beta value- divide by the eigenvalue of the total average contact matrix
Beta0 <- getbeta(R0t = 3.5, pars = pars, p_age = Irlpop$propage, CONTACTMATRIX = contacts_ireland)
p_age <- Irlpop$propage

constraints <- list(home = diag(1,16,16),
           work = diag(1,16,16),
           school = diag(1,16,16),
           others = diag(1,16,16))

Csym <- lapply(contacts_ireland, function(x, p_age) (x + t(x)*((p_age)%*%t(1/p_age)))/2, p_age) # make sure contacts are reciprocal

C <- constraints[[1]]%*%Csym[[1]]+
  constraints[[2]]%*%Csym[[2]]+
  constraints[[3]]%*%Csym[[3]]+
  constraints[[4]]%*%Csym[[4]]

n <- dim(C)[1]
M <- C
for(i in 1:n)
{
  for(j in 1:n){
    M[i,j] = C[i,j]*p_age[i]/p_age[j]
  }
}
eig <- eigen(M)


beta_C <- Beta/max(Re(eig$values))

max(Re(eig$values))

max(Re(eigen(C)$values))

max(Re(eigen(contacts_ireland[[5]])$values))

all.equal(eigen(C)$vectors, eig$vectors)
all.equal(eigen(C)$values, eig$values)

all.equal(contacts_ireland[[5]], Csym[[5]])

contacts_ireland[[5]][1,2] %*% as.matrix(Irlpop$popage[2])

contacts_ireland[[5]][2,1] %*% as.matrix(Irlpop$popage[1])


Csym[[5]][1,2] %*% as.matrix(Irlpop$popage[2])

Csym[[5]][2,1] %*% as.matrix(Irlpop$popage[1])


################################################################################

plot(jg_dat$Beta[1:tmax], type = "l", lwd = 2, 
     xlab ="Time(days)",
     ylab = "Beta",
     panel.first = rect(c(tStartSchoolClosure, tStartIntenseIntervention, tEndIntenseIntervention,
                          tRelaxIntervention1, tRelaxIntervention2, tRelaxIntervention3,
                          tRelaxIntervention4), -1e6,
                        c(tStartIntenseIntervention,tEndIntenseIntervention, tRelaxIntervention1,
                          tRelaxIntervention2,tRelaxIntervention3,
                          tRelaxIntervention4,tRelaxIntervention5 ), 1e6,
                        col=c('gray63','gray48', 'gray70', 'gray75','gray80','gray85', 'gray90'), border=NA))

lines(step_beta, lwd  = 2, col = "red")
lines(beta_C, col = "blue", lwd = 2)
legend(120, 4,legend = c("JG's beta", "Averaged JG beta","Our Beta"),
       col = c("black", "red", "blue"), bty = 'n',lty = c(1,1),lwd = c(2,2), cex = 1)

################################################################################

Beta_scalar <- beta_C/beta_C[1]

scalars <- unique(Beta_scalar)
length(scalars)

################################################################################

jg_beta_C <- jg_dat$Beta/max(Re(eig$values))


plot(jg_dat$Beta[1:tmax], type = "l", lwd = 2)
lines(step_beta, lwd  = 2, col = "red")
lines(beta_C, col = "blue", lwd = 2)
lines(jg_beta_C, col = "green", lwd = 2)


length(jg_beta_C[1:225])
length(seq(1,225,by = 1))

################################################################################

# To get the our equivalent beta value- divide by the eigenvalue of the total average contact matrix
p_age <- Irlpop$propage

Csym <- lapply(contacts_ireland, function(x, p_age) (x + t(x)*((p_age)%*%t(1/p_age)))/2, p_age) # make sure contacts are reciprocal


eig = eigen(Csym[[5]])


beta_Csym <- Beta/max(Re(eig$values))

scalars_Csym <- unique(beta_Csym/beta_Csym[1])
length(scalars_Csym)

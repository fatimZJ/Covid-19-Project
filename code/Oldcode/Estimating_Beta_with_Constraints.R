
## Load the data
source('code/1_loadData.r')
contacts_ireland <- contacts

tmax <- 225
## Load the get beta function 
source("code/getbeta.R")

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

## Constraints

dateStartSchoolClosure = as.Date('2020-03-12') # School closure
dateStartIntenseIntervention = as.Date('2020-03-27')  #Intense intervention
dateEndIntenseIntervention = as.Date('2020-05-18') #date we begin relaxing intense intervention 
dateStart = as.Date('2020-02-28') #start date for epidemic in Ireland

numWeekStagger = c(3,6,9,12,15)

tStartSchoolClosure = as.vector(dateStartSchoolClosure - dateStart)+1 #Time point to add the school closure effect
tStartIntenseIntervention = as.vector(dateStartIntenseIntervention - dateStart)+1 # #Time point to add the intense lockdown effect
tEndIntenseIntervention = as.vector(dateEndIntenseIntervention - dateStart)+1     # for pw = 0.1
tRelaxIntervention1 = tEndIntenseIntervention + numWeekStagger[1]*7                               # for pw = 0.25
tRelaxIntervention2 = tEndIntenseIntervention + numWeekStagger[2]*7                               # for pw = 0.5
tRelaxIntervention3 = tEndIntenseIntervention + numWeekStagger[3]*7
tRelaxIntervention4 = tEndIntenseIntervention + numWeekStagger[4]*7                               # for pw = 0.25
tRelaxIntervention5 = tEndIntenseIntervention + numWeekStagger[5]*7                               # for pw = 0.5

## Specifying interventions and other model parameters
loadInterventions <- function(p_workopen)
{
  list(
    # constraints under a DO-NOTHING scenario 
    base =list(home = diag(1,16,16),
               work = diag(1,16,16),
               school = diag(1,16,16),
               others = diag(1,16,16)),
    # constraints under school closure + some social distancing for school-age going children but 100% workplace
    schcloseonly = list(home = diag(c(rep(1,4),rep(1,12))),
                        work = diag(1,16,16),
                        school = diag(0,16,16),
                        others = diag(c(rep(0.5,4),rep(1,12)))), 
    
    # constraints under work place distancing + schoolclosure 
    schcloseworkplacedist = list(home = diag(1,16,16),
                                 work = diag(p_workopen,16,16),
                                 school = diag(0,16,16),
                                 others = diag(c(rep(0.0,4),rep(0.0,12)))),
    
    # Post Outbeak, people still cautious 
    postoutbreak = list(home = diag(1,16,16),
                        work = diag(0.8,16,16),
                        school = diag(1.0,16,16),
                        others = diag(c(rep(0.8,4),rep(0.5,12)))))
  
}

## Model Parameters

pars <- c(4.9, 5.9, 7.0, 0.25, 0.05, 0.05, 0.5, 0.75, 0.13, 3.6)
names(pars) <- c("L","Cv","Dv","h","i","j","f","tv","q","TT")

Beta <- vector()
R <- vector()
## Get beta for each time point
for (t in 1:225){
constraintsIntervention = loadInterventions(p_workopen = 1)


if(t < tStartSchoolClosure)  
{
  
  INTERVENTION = "base"
  CONSTRAINT = constraintsIntervention$base
}
# I1:  When school winter break but before lockdown period, use 'schcloseonly'
if(t >= tStartSchoolClosure & t < tStartIntenseIntervention) 
{
  INTERVENTION = "schcloseonly"   
  CONSTRAINT = constraintsIntervention[[INTERVENTION]] 
}  

if(t >= tStartIntenseIntervention & t < tEndIntenseIntervention) 
{ 
  INTERVENTION = "schcloseworkplacedist"   
  CONSTRAINT = loadInterventions(p_workopen = 0.01)[[INTERVENTION]]
}  
if(t >= tEndIntenseIntervention & t < tRelaxIntervention1) 
{ 
  INTERVENTION = "schcloseworkplacedist"   
  CONSTRAINT = loadInterventions(p_workopen = 0.15)[[INTERVENTION]]
  #CONSTRAINT = constraintsIntervention[[INTERVENTION]] 
}  

if(t >= tRelaxIntervention1 & t < tRelaxIntervention2) 
{
  INTERVENTION = "schcloseworkplacedist"   
  CONSTRAINT = loadInterventions(p_workopen = 0.25)[[INTERVENTION]] 
}  

if(t >= tRelaxIntervention2 & t < tRelaxIntervention3) 
{
  INTERVENTION = "schcloseworkplacedist"   
  CONSTRAINT = loadInterventions(p_workopen = 0.40)[[INTERVENTION]]
}  

if(t >= tRelaxIntervention3 & t < tRelaxIntervention4) 
{
  INTERVENTION = "schcloseworkplacedist"   
  CONSTRAINT = loadInterventions(p_workopen = 0.55)[[INTERVENTION]]
}  

if(t >= tRelaxIntervention4 & t < tRelaxIntervention5) 
{
  INTERVENTION = "schcloseworkplacedist"   
  CONSTRAINT = loadInterventions(p_workopen = 0.70)[[INTERVENTION]]
}  
if(t >= tRelaxIntervention5)  
{
  INTERVENTION = "postoutbreak"
  CONSTRAINT = constraintsIntervention$postoutbreak
}

## Estimating Beta
R[t] <- getR0t(beta = jg_dat$Beta[1], constraints = CONSTRAINT, pars = pars, p_age = Irlpop$propage, CONTACTMATRIX = contacts_ireland)

Beta[t] <- getbeta(R0t = 3.4, pars = pars, p_age = Irlpop$propage, CONTACTMATRIX = contacts_ireland)

}


plot(Beta, type = "l", lwd = 2)

#lines(jg_dat$Beta[1:tmax])

plot(jg_dat$Beta[1:tmax], type = "l", lwd = 2)
lines(Beta, col = "red", lwd = 2)


## Fitting a model with these beta values and no contacts matrices should give the 
## same output as fitting a model with a fixed beta value and these contact matrices

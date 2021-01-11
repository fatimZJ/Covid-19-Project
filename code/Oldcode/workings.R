


scalars_init <- c(1.362216848, 0.510116460, 0.106050230, 0.008395608, 0.122922232,
  0.197667229, 0.248520907, 0.171967661, 0.091010276)
names(scalars_init) <- names(scalars)


Base$tstart_intervention -> tstart_intervention
Base$tend_intervention -> tend_intervention
scalars <- Base$scalars

CONSTRAINT <- vector()

for ( j in 1:length(Base$sol_out$time)) {
  t <- Base$sol_out$time[j]
  for( i in 1:length(tstart_intervention)){
    if( (t > ((tstart_intervention[[i]]) - 1)) & (t <= (tend_intervention[[i]])))
    {
      INTERVENTION <- noquote(names(tstart_intervention)[i])
      CONSTRAINT[j]  <- as.numeric(scalars[INTERVENTION])
    } else {
    INTERVENTION <- noquote(names(tstart_intervention)[1])
    CONSTRAINT[j]  <- as.numeric(scalars[INTERVENTION])
  }
}
}

length(Base$sol_out$time)
length(CONSTRAINT)
plot(3.4*CONSTRAINT, type = "l")
abline( h = 1, col= "red")


if ((t > (as.numeric(tstart_intervention) - 1)) & (t <= (tend_intervention))) {
  INTERVENTION <- noquote(names(tend_intervention))
  CONSTRAINT[j]  <- as.numeric(scalars[INTERVENTION])
}


for ( j in 1:length(Base$sol_out$time)) {
  t <- Base$sol_out$time[j]
  if (t == 0) {
    INTERVENTION <- noquote(names(tstart_intervention)[1])
    CONSTRAINT[j]  <- as.numeric(scalars[INTERVENTION])
  } else { for( i in 1:length(tstart_intervention)){
    if( (t > ((tstart_intervention[[i]]) - 1)) & (t <= (tend_intervention[[i]])))
    {
      INTERVENTION <- noquote(names(tstart_intervention)[i])
      CONSTRAINT[j]  <- as.numeric(scalars[INTERVENTION])
    } 
  }
  }
}

length(Base$sol_out$time)
length(CONSTRAINT)
plot(3.4*CONSTRAINT, type = "l")
abline( h = 1, col= "red")


################################################################################
dateStart <- interventions_info$start[1]
dateEnd <- interventions_info$end[dim(interventions_info)[1]]

#tstart_intervention -> interventions_info$start
#tend_intervention -> interventions_info$end

tstart_intervention <- list()
tend_intervention <- list()

for ( i in 1:length(interventions$policy)) {
  
  tstart_intervention[[i]] <- as.vector(interventions$start[i] - dateStart) 
  names(tstart_intervention)[[i]] <- interventions$policy[i]
  tend_intervention[[i]] <- as.vector(interventions$end[i] - dateStart) 
  names(tend_intervention)[[i]] <- interventions$policy[i]
}

CONSTRAINT <- vector()
INTERVENTION <- vector()

for ( j in 1:length(Base$sol_out$time)) {
  t <- Base$sol_out$time[j]
for( i in 1:length(tstart_intervention)){
  if( (t >= ((tstart_intervention[[i]]))) & (t < ((tend_intervention[[i]]) + 1)))
  {
    INTERVENTION[j] <- noquote(names(tstart_intervention)[i])
    CONSTRAINT[j]  <- as.numeric(scalars[INTERVENTION[j]])
  } else {
    INTERVENTION[j] <- noquote(names(tstart_intervention)[1])
    CONSTRAINT[j]  <- as.numeric(scalars[INTERVENTION[j]])
  }
}
}



test1 <- cbind(INTERVENTION, CONSTRAINT)


CONSTRAINT <- vector()
INTERVENTION <- vector()

tstart_intervention <- as.numeric(difftime(interventions$start, dateStart, units = "days"))
names(tstart_intervention) <- interventions$policy

tend_intervention <- as.numeric(difftime(interventions$end, dateStart, units = "days"))
names(tend_intervention) <-  interventions$policy


for ( j in 1:length(Base$sol_out$time)) {
  t <- Base$sol_out$time[j]
ind <- ((t >= ((tstart_intervention))) & (t < (as.numeric(tend_intervention) + 1)))

INTERVENTION[j] <- noquote(names(ind[ind == 1]))

CONSTRAINT[j] <- ifelse( !any(ind), scalars["No Intervention"],  scalars[INTERVENTION[j]])


}

test2 <- cbind(INTERVENTION, CONSTRAINT)

all.equal(test1, test2)


### Use R0 to get beta
getbeta <- function(R0t, pars, constraints, p_age, calculate_transmission_probability = TRUE, CONTACTMATRIX = contacts) {
  
  ### Return arbitrary beta if no calculation is desired
  if (!calculate_transmission_probability) { return(0.025) }
  
  ### Extract Parameters
  h <- pars["h"]
  i <- pars["i"]
  j <- pars["j"]
  L <- pars["L"]
  Cv <- pars["Cv"]
  Dv <- pars["Dv"]
  f <- pars["f"]
  q <- pars["q"]
  tv <- pars["tv"]
  TT <- pars["TT"]
  
  ### Prepare population contact matrix and apply constraints
  n <- dim(CONTACTMATRIX[[1]])[1]
  if (missing(constraints)) {
    constraints <- list(home = diag(1,n),
                        work = diag(1,n), 
                        school = diag(1,n), 
                        others = diag(1,n))
  }
  
  Csym <- lapply(CONTACTMATRIX, function(x, p_age) (x + t(x)*((p_age)%*%t(1/p_age)))/2, p_age) # make sure contacts are reciprocal
  
  C <- constraints[[1]]%*%Csym[[1]]+
    constraints[[2]]%*%Csym[[2]]+
    constraints[[3]]%*%Csym[[3]]+
    constraints[[4]]%*%Csym[[4]]
  
  ### Create the N matrix
  # (this part is horribly inefficient, lots of room for improvement)
  transmission_rates <- c(0, 1, h, i, 1, j, 1)
  N_vals <- matrix(0, nrow = n, ncol = n*7)
  col_ind <- seq(1, n*7, 7)
  for(i in 1:n) {
    count <- 1
    for(j in 1:n) {
      k <- col_ind[count]
      N_vals[i,k:(k+6)] = transmission_rates * C[i,j] * p_age[i]/p_age[j]
      count <- count + 1
    }
  }
  N <- sparseMatrix( dims = rep(n*7, 2), i = rep(seq(1, 7*n, 7), n*7), j = rep(seq(1, n*7), each = n), x = as.vector(N_vals) )
  
  ### Create the inverse V matrix
  V_inv_vals <- c(L, (1-f)*(Cv-L), Cv-L, f*Dv, Dv, (1-f)*q*(Dv-Cv+L),
                  q*(Dv-Cv+L), Dv-Cv+L, TT*tv*(1-f), TT*tv, TT,
                  tv*(1-f)*(Dv-Cv+L-TT), tv*(Dv-Cv+L-TT), Dv-Cv+L-TT,
                  Dv-Cv+L-TT, (1-q-tv)*(1-f)*(Dv-Cv+L), (1-q-tv)*(Dv-Cv+L),
                  Dv-Cv+L)
  V_inv_i <- c(1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7)
  V_inv_j <- c(1, 1, 2, 1, 3, 1, 2, 4, 1, 2, 5, 1, 2, 5, 6, 1, 2, 7)
  V_inv_list <- rep( list(spMatrix(nrow = 7, ncol = 7, i = V_inv_i, j = V_inv_j, x = V_inv_vals)), n)
  
  V_inv <- bdiag(V_inv_list)
  
  ### Calculate and return beta
  eig <- eigen(N%*%V_inv)
  R0t / max(Re(eig$values))
  
  #0.3
}


Irl_pop <- read.csv("Data/Ireland_pop.csv")    #population of Ireland! versus population of County Dublin
JG_SEIR_out <- read.csv("Data/dat_seir_code.csv", row.names = NULL)
dim(JG_SEIR_out)
JG_SEIR_out <- t(JG_SEIR_out)
colnames(JG_SEIR_out)
rownames(JG_SEIR_out)

JG_SEIR_out <- cbind(rownames(JG_SEIR_out), data.frame(JG_SEIR_out, row.names=NULL))
colnames(JG_SEIR_out) <- c(JG_SEIR_out[1,])
JG_SEIR_out <- JG_SEIR_out[-1,]

head(JG_SEIR_out)


Rt <- as.numeric(JG_SEIR_out[,9])
pars <- c(4.9, 5.9, 7.0, 0.25, 0.05, 0.05, 0.5, 0.75, 0.13, 3.6)
names(pars) <- c("L","Cv","Dv","h","i","j","f","tv","q","TT") 

constraints <- list(home = diag(1,16,16),
                    work = diag(1,16,16),
                    school = diag(1,16,16),
                    others = diag(1,16,16))

p_age <- Irl_pop$propage
CONTACTMATRIX <- contacts

Bt_d <- vector()
for (i in 1:length(Rt)){
 
beta <- getbeta(R0t = Rt[i], pars = pars, constraints = constraints, p_age = p_age, 
        calculate_transmission_probability = TRUE, CONTACTMATRIX = CONTACTMATRIX) 
Bt_d[i] <- beta
}
  
plot(as.numeric(JG_SEIR_out$Beta), col = "red", type = "l", lwd = 1.5)
lines(Bt_d, type = "l", lwd = 1.5)

par(mfrow=c(1,1))

################################################################################
Csym <- lapply(CONTACTMATRIX, function(x, p_age) (x + t(x)*((p_age)%*%t(1/p_age)))/2, p_age) # make sure contacts are reciprocal

C <- constraints[[1]]%*%Csym[[1]]+
  constraints[[2]]%*%Csym[[2]]+
  constraints[[3]]%*%Csym[[3]]+
  constraints[[4]]%*%Csym[[4]]
n <- dim(C)[1]
M = C
for(i in 1:n)
{
  for(j in 1:n){
    M[i,j] = C[i,j]*p_age[i]/p_age[j]
  }
}
eig = eigen(M)
  
BTest = as.numeric(JG_SEIR_out$Beta)/max(Re(eig$values))


plot(as.numeric(JG_SEIR_out$Beta), col = "red", type = "l", lwd = 1.5)
lines(Bt_d, type = "l", lwd = 1.5)
lines(BTest, type = "l", col = "grey", lwd = 1.5)

################################################################################
### Extract Parameters
h <- pars["h"]
i <- pars["i"]
j <- pars["j"]
L <- pars["L"]
Cv <- pars["Cv"]
Dv <- pars["Dv"]
f <- pars["f"]
q <- pars["q"]
tv <- pars["tv"]
TT <- pars["TT"]
g <- 1 - f
# Compartment Contributions
IaC <- prod(f, Dv, h)
IpC <- prod(g, Cv - L)
IiC <- prod(g, q, Dv - Cv + L, i)
It1C <- prod(g, tv, TT)
It2C <- prod(g, tv, Dv - Cv + L - TT, j)
InC <- prod(g, 1 - q - tv, Dv - Cv + L)
allC <- sum(IaC, IpC, IiC, It1C, It2C, InC)

#R0t = 2
bt = Rt/(max(Re(eig$values))*allC)

plot(as.numeric(JG_SEIR_out$Beta), col = "red", type = "l", lwd = 2)
lines(Bt_d, type = "l", lwd = 2)
lines(BTest, type = "l", col = "grey", lwd = 2)
lines(bt, type = "l", col = "green", lwd = 2)


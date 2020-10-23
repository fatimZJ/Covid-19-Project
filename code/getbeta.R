### get beta

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
  
}

### Use beta to get R0
getR0t <- function(beta, pars, constraints, p_age, calculate_transmission_probability = TRUE, CONTACTMATRIX = contacts) {
  
  ### Return arbitrary beta if no calculation is desired
  if (!calculate_transmission_probability) { return(0.999) }
  
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
  beta * max(Re(eig$values))
  
}


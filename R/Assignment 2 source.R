# Rejection_Density <- Vectorize(function(y, data=Poisson){
#   x <- data$x
#   z <- data$z
#   prod(exp(y*z*x-exp(y*x)))
# })
#
# Rejection_Density2 <- function(y, data=Poisson){
#   x <- data$x
#   z <- data$z
#   exp(sum(y*z*x-exp(y*x)))
# }
#
# rejec_dens1 <- Vectorize(function(y){
#   xy <- y * poisson$x
#   prod(exp(xy * poisson$z - exp(xy)))
# })
#
# rejec_dens5 <- function(y){
#   term1 <- y*sum(poisson$x * poisson$z)
#   term2 <- rowSums(exp(outer(y, poisson$x)))
#   exp(term1 - term2)
# }

rejec_dens1 <- Vectorize(function(y, data=Poisson){
  x <- data$x
  z <- data$z
  prod(exp(y*z*x-exp(y*x)))
})

rejec_dens3 <- Vectorize(function(y){
  xy <- y * poisson$x
  prod(exp(xy * poisson$z - exp(xy)))
})

rejec_dens5 <- function(y){
  term1 <- y*sum(poisson$x * poisson$z)
  term2 <- rowSums(exp(outer(y, poisson$x)))
  exp(term1 - term2)
}


### Rejection sampler ###

sekvens <- seq(0,0.5,length.out=50)
out <- sapply(sekvens, function(x) rejec_dens5(x))
mu <- function(a, b, n){
  sekvens <- seq(a,b,length.out=n)
  out <- sapply(sekvens, function(x) rejec_dens5(x))
  sekvens[which.max(out)]
}
p <- function(y){
  mu <- mu(0,0.5,50)
  exp(-0.5*(y-mu)^2) #mu=0.24 is the maksimum value of the "density"
}
q <- function(x) rejec_dens5(x)
alpha_prime <- min(sapply(sekvens, function(y) p(y)/q(y)))


reject_sampler = function(n, rProposal = rnorm, dProposal = p, dTarget = q) {
  result = numeric(n)
  mu <- mu(0,0.5,50)
  for(i in 1:n) {
    U = 1
    ratio = 0
    while(U > ratio) {
      U = runif(1)
      Y = rProposal(1,mean=mu)
      ratio = alpha_prime * dTarget(Y) / dProposal(Y)
    }
    result[i] = Y
  }
  result
}

### Optimiseringsproces ###
Rcpp::sourceCpp("C:/Users/Jakob/OneDrive - University of Copenhagen/Universitetet/Forsikringsmatematik/År 5/CompStat/C++ scripts/Sapply_cpp.cpp")

mu_opt <- function(a, b, n){
  sekvens <- seq(a,b,length.out=n)
  out <- sapply_cpp(sekvens, function(x) rejec_dens5(x))
  sekvens[which.max(out)]
}
p_opt <- function(y){
  mu <- mu_opt(0,0.5,50)
  exp(-0.5*(y-mu)^2) #mu=0.24 is the maksimum value of the "density"
}

alpha_prime_opt <- min(sapply_cpp(sekvens, function(y) p(y)/q(y)))

ratio_vec <- c(0)

reject_sampler_opt = function(n, rProposal = rnorm, dProposal = p_opt, dTarget = q, trace = FALSE) {
  count <- 0
  result = numeric(n)
  mu <- mu_opt(0,0.5,50)
  for(i in 1:n) {
    U <- 1
    ratio <- 0
    while(U > ratio) {
      count <- count + 1
      U <- runif(1)
      Y <- rProposal(1,mean=mu)
      ratio <- alpha_prime_opt * dTarget(Y) / dProposal(Y)
      ratio_vec <<- append(ratio_vec, ratio)
    }
    result[i] <- Y
  }
  if(trace){
    cat("Rejection frequency =", (count-n)/count, "\n")
  }
  result
}


fkt_stream <- function(m, fkt, ...) {
  args <- list(...)
  cache <- do.call(fkt, c(m, args)) # Calls the function with n=m and arguments at                                     # told. E.g. ... could be mean and sd for                                         # fkt=rnorm
  j <- 0
  fact <- 1
  next_rn <- function(r = m) {
    j <<- j + 1
    if(j > m) {                     # When all the generated RVs are used we                                          # generates new ones.
      if(fact == 1 && r < m) fact <<- m / (m - r)
      m <<- floor(fact * (r + 1))
      cache <<- do.call(fkt, c(m, args))
      j <<- 1
    }
    cache[j]                        # The fkt_stream function reterns the jth                                         # element of the RV vector
  }
  next_rn
}

reject_sampler_opt2 = function(n, rProposal = rnorm, dProposal = p_opt, dTarget = q, trace = FALSE) {
  count <- 0
  result = numeric(n)
  mu <- mu_opt(0,0.5,50)

  #Samling the random variables
  U <- fkt_stream(m=n, fkt=runif)
  Y <- fkt_stream(m=n, rProposal, mean=mu)

  for(i in 1:n) {
    reject <- TRUE
    while(reject) {
      count <- count + 1
      y <- Y(n-i)
      ratio <- alpha_prime_opt * dTarget(y) / dProposal(y)
      reject <- (U(n-i) > ratio)
      ratio_vec <<- append(ratio_vec, ratio)
    }
    result[i] <- y
  }
  if(trace){
    cat("Rejection frequency =", (count-n)/count, "\n")
  }
  result
}




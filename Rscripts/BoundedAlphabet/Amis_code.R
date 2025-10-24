Rcpp::sourceCpp("../../src/RcppFunctions.cpp")

## Ami's code from matlab to R
set.seed(123)  # for reproducibility

# Define vectors
n_vec <- c(seq(500, 2000, by = 100), seq(2000, 10000, by = 500))
m <- 10000
p <- runif(m)
alph <- 0.05

benchmark <- log(m / alph) / n_vec

mine_analytical <- numeric(length(n_vec))
mine_analytical_cpp <- numeric(length(n_vec))
lb <- numeric(length(n_vec))
lb_cpp <- numeric(length(n_vec))
ok_vec <- numeric(length(n_vec))
num_of_exp <- 100

n_ind = 1
for (n_ind in seq_along(n_vec)) {
  n <- n_vec[n_ind]
  
  exp_ind = 1
  for (exp_ind in 1:num_of_exp) {
    # Generate r: binomial random numbers normalized by n
    r <- rbinom(m, size = n, prob = p) / n
    
    # Lower bound computation
    current_lb <- 0
    current_lb_cpp <- 0
    if(any(r == 0)) {
      current_lb <- log(sum(r == 0) / alph) / n
      current_lb_cpp <- compute_LB_analytical(n,r*n,alph)
      
    }
    lb[n_ind] <- lb[n_ind] + current_lb
    lb_cpp[n_ind] <- lb_cpp[n_ind] + current_lb_cpp
    
    # Analytical bound computation
    a1 <- 0.99 * alph
    a2 <- 0.01 * alph
    # a <- (log(n))^2
    a <- log(n)
    m_p <- sum((1 - r)^a) + a * sqrt(m * log(1 / a2) / n)
    current_analytical <- log(m_p / a1) / (n - a)
    mine_analytical[n_ind] <- mine_analytical[n_ind] + current_analytical
    
    current_analytical_cpp <- compute_UB_analytical( n, r*n, m, a, alph, FALSE)
    mine_analytical_cpp[n_ind] <- mine_analytical_cpp[n_ind] + current_analytical_cpp
    
    # Check condition
    if (any(r == 0)) {
      ok <- 0
      if (max(p[r == 0]) <= current_analytical) {
        ok <- 1
      }
      ok_vec[n_ind] <- ok_vec[n_ind] + ok
    } else {
      ok_vec[n_ind] <- ok_vec[n_ind] + 1
    }
  }
}



plot(x = n_vec, y = benchmark, type = "l", col = "black", lwd = 2, ylim = c(0,0.025))
points(x = n_vec, y = mine_analytical/num_of_exp, type = "l", lwd = 2, col = "red")
points(x = n_vec, y = mine_analytical_cpp/num_of_exp, type = "l", lwd = 2, col = "red")
points(x = n_vec, y = lb/num_of_exp, type = "l", lwd = 2, col = "green")
points(x = n_vec, y = lb_cpp/num_of_exp, type = "l", lwd = 2, col = "green")


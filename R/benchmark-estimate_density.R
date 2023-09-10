
benchmark_estimate_density <- function(){

  # Instantiate
  set.seed(0)
  n_seq <- 3^seq(4, 10)
  n <- max(n_seq)
  sim_bern <- rbinom(n, 1, 0.5)
  sim_normal <- rnorm(n)
  sim_mix <- 6*sim_bern + sim_normal
  set.seed(NULL)

  # Create containers
  benchmarks <- vector("list", length(n_seq))

  # Do benchmarks
  for (i in seq_along(n_seq)){
    print("Finished a Benchmark")
    benchmarks[[i]] <- microbenchmark::microbenchmark(
      iter_dens_est(sim_normal[1:n_seq[i]], maxiter = 3),
      iter_dens_est(sim_mix[1:n_seq[i]], maxiter = 3),
      iter_dens_est(sim_normal[1:n_seq[i]], maxiter = 10),
      iter_dens_est(sim_mix[1:n_seq[i]], maxiter = 10),
      iter_dens_est(sim_normal[1:n_seq[i]], maxiter = 30),
      iter_dens_est(sim_mix[1:n_seq[i]], maxiter = 30),
      times = 3
    )
  }

  return(benchmarks)

  # Draw plots
  # Later

}

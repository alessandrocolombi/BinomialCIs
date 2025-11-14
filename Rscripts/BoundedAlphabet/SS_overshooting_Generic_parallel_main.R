wd = "C:/Users/colom/BinomialCIs/Rscripts/BoundedAlphabet/"
# wd = "~/Lucia/Ale/BinomialCIs/Rscripts/BoundedAlphabet/"
setwd(wd)

# The path to the script to execute in parallel
scriptpath = "./SS_overshooting_Generic.R"

# Parameters
params = c(0.25,0.5,1.02) # zipfs case
# params = c(1,0.1,0.001)   # Uniform case


SimNo <- seq_along(params)
for(i in seq_along(SimNo)){
  # Parallel execution 
  param <- params[SimNo[i]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
}



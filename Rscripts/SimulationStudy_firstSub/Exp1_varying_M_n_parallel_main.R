wd = "C:/Users/colom/BinomialCIs/Rscripts/SimulationStudy_firstSub/"
setwd(wd)

# The path to the script to execute in parallel
scriptpath = "./Exp1_varying_M_n.R"

# Constant
params = c(0.25,0.5,1.02) # zipfs case
# params = c(0.005,0.1,0.25) # geom case
# params = c(0.001,0.05,0.5) # constant case


SimNo <- seq_along(params)
for(i in seq_along(SimNo)){
  # Parallel execution 
  param <- params[SimNo[i]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
}


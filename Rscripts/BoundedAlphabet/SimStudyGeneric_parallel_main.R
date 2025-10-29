wd = "C:/Users/colom/BinomialCIs/Rscripts/BoundedAlphabet/"
# wd = "~/Lucia/Ale/BinomialCIs/Rscripts/BoundedAlphabet/"
setwd(wd)

# The path to the script to execute in parallel
scriptpath = "./SimStudyGeneric.R"

# Constant
params = c(0.1,0.25,0.5,0.75,1.02) # zipfs case
params = c(0.001,0.005,0.015,0.05,0.1,0.25) # geom case
params = c(0.001,0.01,0.05,0.25,0.5) # constant case


SimNo <- seq_along(params)
for(i in seq_along(SimNo)){
  # Parallel execution 
  param <- params[SimNo[i]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
}



param <- params[SimNo[2]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
param <- params[SimNo[3]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
param <- params[SimNo[4]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)

param <- params[SimNo[5]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
param <- params[SimNo[6]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)

param <- params[SimNo[7]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
param <- params[SimNo[8]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
param <- params[SimNo[9]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
param <- params[SimNo[11]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
param <- params[SimNo[12]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
param <- params[SimNo[13]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
param <- params[SimNo[14]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
param <- params[SimNo[15]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
param <- params[SimNo[16]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)

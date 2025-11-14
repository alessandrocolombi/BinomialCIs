wd = "C:/Users/colom/BinomialCIs/Rscripts/UnboundedAlphabet/"
# wd = "~/Lucia/Ale/BinomialCIs/Rscripts/UnboundedAlphabet/"
setwd(wd)

# The path to the script to execute in parallel
scriptpath = "./SSUbb_3PBetaPr.R"

idxs = c(1)
for(i in seq_along(idxs)){
  # Parallel execution 
  idx <- idxs[i]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
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

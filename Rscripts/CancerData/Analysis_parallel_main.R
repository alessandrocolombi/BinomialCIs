wd = "C:/Users/colom/BinomialCIs/Rscripts/CancerData/"
wd = "~/Lucia/Ale/BinomialCIs/Rscripts/CancerData/"
setwd(wd)

# The path to the script to execute in parallel
scriptpath = "./Analysis_exec.R"
Datasets = 1:32
SimNo <- 17:32

# Parallel execution 
idx <- Datasets[SimNo[1]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- Datasets[SimNo[2]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- Datasets[SimNo[3]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- Datasets[SimNo[4]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- Datasets[SimNo[5]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- Datasets[SimNo[6]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- Datasets[SimNo[7]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- Datasets[SimNo[8]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- Datasets[SimNo[9]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- Datasets[SimNo[11]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- Datasets[SimNo[12]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- Datasets[SimNo[13]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- Datasets[SimNo[14]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- Datasets[SimNo[15]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- Datasets[SimNo[16]]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)

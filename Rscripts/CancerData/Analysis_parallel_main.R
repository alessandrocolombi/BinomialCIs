wd = "C:/Users/colom/BinomialCIs/Rscripts/CancerData/"
wd = "~/Lucia/Ale/BinomialCIs/Rscripts/CancerData"
setwd(wd)

# The path to the script to execute in parallel
scriptpath = "./Analysis_exec.R"
SimNo <- 1:16

# Parallel execution 
idx <- SimNo[1]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[2]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[3]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[4]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[5]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[6]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[7]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[8]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[9]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[11]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[12]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[13]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[14]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[15]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[16]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)

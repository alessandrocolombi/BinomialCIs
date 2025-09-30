wd = "C:/Users/colom/BinomialCIs/Rscripts/CancerData/"
setwd(wd)

# The path to the script to execute in parallel
scriptpath = "./Analysis_exec.R"
SimNo <- 1:10

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
idx <- SimNo[10]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)

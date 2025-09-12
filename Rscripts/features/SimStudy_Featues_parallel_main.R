## Source script to execute samplers on different datasets in parallel ##
## Once set, from R console type source("<path/to/this_file>") ##

# Set proper path
# path <- <path/to/this/file>
# setwd(path)

# The path to the script to execute in parallel
setwd("C:/Users/colom/bnp_upperbounds/Rscripts/features")
# setwd("/home/lucia.paci/Lucia/Ale/bnp_upperbounds/Rscripts/features")

scriptpath = "SimStudy_Featues_exec.R"
SimNo <- 1:15

# Parallel execution (5 cores)
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
idx <- SimNo[11]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[12]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[13]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[14]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[15]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)

## Source script to execute samplers on different datasets in parallel ##
## Once set, from R console type source("<path/to/this_file>") ##

# Set proper path
# path <- <path/to/this/file>
# setwd(path)

# The path to the script to execute in parallel
setwd("G:/.shortcut-targets-by-id/1Ck2MctcmCBWOueSMeW4Zr8IOvOaOzdqz/BicoccaDrive/g-masses/R/Scripts")

scriptpath = "SimStudy_Featues_exec.R"
SimNo <- 6:10

# Parallel execution (5 cores)
idx <- SimNo[1]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[2]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[3]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[4]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)
idx <- SimNo[5]; rstudioapi::jobRunScript(path = scriptpath, importEnv = TRUE)

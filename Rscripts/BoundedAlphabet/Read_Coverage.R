
setwd("C:/Users/colom/BinomialCIs/Rscripts/BoundedAlphabet/save/Coverage")

# Zipfs ---------------------------------------------------------
files = list.files("Zipfs/")
for(i in seq_along(files)){
  cat("\n --- \n")
  nome_file = paste0("Zipfs/",files[i])
  load(nome_file)
  print(round( cov_mat, 3))
}


# Geom ---------------------------------------------------------
files = list.files("Geom/")
for(i in seq_along(files)){
  cat("\n --- \n")
  nome_file = paste0("Geom/",files[i])
  load(nome_file)
  print(round( cov_mat, 3))
}

# Constant ---------------------------------------------------------
files = list.files("Constant/")
for(i in seq_along(files)){
  cat("\n --- \n")
  nome_file = paste0("Constant/",files[i])
  load(nome_file)
  print(round( cov_mat, 3))
}

# Constant ---------------------------------------------------------
files = list.files("Uniform/")
for(i in seq_along(files)){
  cat("\n --- \n")
  nome_file = paste0("Uniform/",files[i])
  load(nome_file)
  print(round( cov_mat, 3))
}
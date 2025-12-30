# Install neonUtilities and neonPlantEcology
# remotes::install_github("NEONStats/NEON-utilities")
# remotes::install_github("admahood/neonPlantEcology")
wd = "C:/Users/colom/BinomialCIs/Rscripts/"
setwd(wd)

library(neonUtilities)
library(neonPlantEcology)

# Download raw plant diversity data for a NEON domain or site
# (e.g., domain "D14" or another domain code)
dpID <- "DP1.10058.001"  # Example DP ID for plant diversity product
raw <- neonUtilities::loadByProduct(dpID = dpID)

# Convert to community matrix (occurrence = binary)
comm_matrix <- neonPlantEcology::npe_community_matrix(raw, binary = TRUE)

# comm_matrix is now your n×T incidence matrix
head(comm_matrix)
dim(comm_matrix)

save(comm_matrix,file = "Neon.Rdat")

# Load and read all --------------------------------------------------------------------
wd = "C:/Users/colom/BinomialCIs/Rscripts/BBSData"
setwd(wd)

Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../R/Rfunctions.R")

# https://github.com/trashbirdecology/bbsAssistant

# Download data -----------------------------------------------------------


# 2. Install the BBS assistant package
remotes::install_github("https://github.com/trashbirdecology/bbsAssistant")

# 3. Load the library
library(bbsAssistant)

library(Matrix)
library(tidyverse)
bbs <- grab_bbs_data()
names(bbs)

# Focus on the observation data
obs_raw <- bbs$observations

save(obs_raw, file = "DataRaw_all.Rdat")
# Select all roads in a single year ----------------------------------------------------
year = 2019
# Select a single year (e.g., 2019) to get a clean 'n x T' snapshot
matrix_prep <- obs_raw %>%
  filter(Year == year) %>%
  select(RouteDataID, AOU) %>% # RouteDataID = Plot (n), AOU = Species ID (T)
  distinct() # Remove duplicates to ensure binary 0/1

# 2. Pivot the data into a Wide Matrix
# This will turn RouteDataID into rows (n) and AOU into columns (T)
incidence_matrix <- matrix_prep %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = AOU, values_from = presence, values_fill = 0) %>%
  column_to_rownames("RouteDataID")

# 3. Check your dimensions
# n = Number of survey routes, T = Number of species
print(dim(incidence_matrix))

save(incidence_matrix, file = "Data2019_allRoutes.Rdat")

data = incidence_matrix
dim(data)
n = nrow(data)
Nj = c(apply(data, 2, sum))
Kobs = length( which(Nj > 0) )
TabNj = tabulate(Nj, nbins = n) 

# Run
ExtrCurve = rarefaction.array(object = data, n_reorderings = 10, seed = 1234)

# Plot
par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
plot(x = 0, y = 0, type = "n",
     main = "", xlab = "#obs.", ylab = "#species",
     ylim = c(0,Kobs+1),
     xlim = c(0,n+1),
     pch = 1) # init plot
polygon( c(1:n, rev(1:n)),
         c(ExtrCurve[1,], rev(ExtrCurve[3,])),
         col = "grey75",
         border = NA) # plot in-sample bands
points(x = 1:n, y = ExtrCurve[2,], type = "l", lwd = 3) # plot mean obs




# Bounded
b_n <- log(n)
Mguess = 10 * Kobs
Nj_guess = c(Nj, rep(0,Mguess - length(Nj) ))
alpha = 0.05
U_bounded <- compute_UB_analytical(n, Nj_guess, Mguess, b_n, alpha, FALSE)
U_bounded

# Unbounded
beta = 1e-5
Shat  <- sum(Nj) / n
Sstar <- ( sqrt( -log(beta) / (2 * n) ) +
             sqrt( Shat + (-log(beta) / (2 * n)) ) )^2
r_n   <- log( Sstar / (-log(1 - alpha + beta)) ) +
         log(n) - log(log(n))
U_unbounded <- compute_UB_rnorm(n, alpha, beta, r_n, Shat)
U_unbounded

# Coverage
C_hat <- coverage_ChaoJost(n, Nj)
C_hat

# Read and save specific data in a single state for a single year ----------------------------------------------------

library(tidyverse)
library(Matrix)

# Configuration

# list of available states
bbsAssistant::region_codes
target_year <- 2019
target_state <- 14  # Select California


# 1. Filter and Reshape
stop_level_data <- bbs$observations %>%
  filter(Year == target_year, StateNum == target_state) %>%
  # Select the ID, Species, and all 50 stop columns
  select(RouteDataID, AOU, Stop1:Stop50) %>%
  # Pivot the stops from columns into rows
  pivot_longer(cols = Stop1:Stop50, 
               names_to = "StopNumber", 
               values_to = "Count") %>%
  # Convert to binary: keep only rows where bird was present
  filter(Count > 0) %>%
  # Create a unique ID for each Plot (Route + Stop)
  mutate(PlotID = paste0(RouteDataID, "_", StopNumber)) %>%
  select(PlotID, AOU) %>%
  distinct()

# 2. Create the n x T Sparse Matrix
# This ensures R doesn't crash from memory usage
incidence_matrix_large <- stop_level_data %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = AOU, 
              values_from = presence, 
              values_fill = 0) %>%
  column_to_rownames("PlotID") %>%
  as.matrix()

# 3. Verify the size
print(dim(incidence_matrix_large))

save(incidence_matrix_large, file = paste0("data/BBS_2019_",target_state,".Rdat"))



# Extr Curves -------------------------------------------------------------
target_state = 14
load(paste0("data/BBS_2019_",target_state,".Rdat"))
data = incidence_matrix_large
dim(data)
n = nrow(data)
Nj = c(apply(data, 2, sum))
Kobs = length( which(Nj > 0) )
TabNj = tabulate(Nj, nbins = n) 

# Run
ExtrCurve = rarefaction.array(object = data, n_reorderings = 10, seed = 1234)

# Plot
par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
plot(x = 0, y = 0, type = "n",
     main = "", xlab = "#obs.", ylab = "#species",
     ylim = c(0,Kobs+1),
     xlim = c(0,n+1),
     pch = 1) # init plot
polygon( c(1:n, rev(1:n)),
         c(ExtrCurve[1,], rev(ExtrCurve[3,])),
         col = "grey75",
         border = NA) # plot in-sample bands
points(x = 1:n, y = ExtrCurve[2,], type = "l", lwd = 3) # plot mean obs




# Bounded
b_n <- log(n)
Mguess = 10 * Kobs
Nj_guess = c(Nj, rep(0,Mguess - length(Nj) ))
alpha = 0.05
U_bounded <- compute_UB_analytical(n, Nj_guess, Mguess, b_n, alpha, FALSE)
U_bounded

# Unbounded
beta = 1e-5
Shat  <- sum(Nj) / n
Sstar <- ( sqrt( -log(beta) / (2 * n) ) +
             sqrt( Shat + (-log(beta) / (2 * n)) ) )^2
r_n   <- log( Sstar / (-log(1 - alpha + beta)) ) +
  log(n) - log(log(n))
U_unbounded <- compute_UB_rnorm(n, alpha, beta, r_n, Shat)
U_unbounded

# Coverage
C_hat <- coverage_ChaoJost(n, Nj)
C_hat

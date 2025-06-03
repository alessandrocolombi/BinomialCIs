setwd("C:/Users/colom/bnp_upperbounds/Rscripts")
# setwd("/home/lucia.paci/Lucia/Ale/bnp_upperbounds/Rscripts")

# Librerie ----------------------------------------------------------------

suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
Rcpp::sourceCpp("../src/RcppFunctions.cpp")
source("../R/Rfunctions.R")


# Species - Zipfs ---------------------------------------------------------
load("save/SimStudySpecies_coverages_zipfs_all.dat")
round( t(coverages), 3)
# Species - NegBin ---------------------------------------------------------
load("save/SimStudySpecies_coverages_NegBin_all.dat")
round( t(coverages), 3)
# Species - Geom ---------------------------------------------------------
load("save/SimStudySpecies_coverages_geom_all.dat")
round( t(coverages), 3)
# Species - Unif ---------------------------------------------------------
load("save/SimStudySpecies_coverages_unif_all.dat")
round( t(coverages), 3)
# Features - Zipfs ---------------------------------------------------------
load("save/SimStudyFeatures_coverages_zipfs_all.dat")
round( t(coverages), 3)
# Features - Zipfs ---------------------------------------------------------
load("save/SimStudyFeatures_coverages_geom_all.dat")
round( t(coverages), 3)


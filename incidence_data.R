
# SpadeR - Incidence data -------------------------------------------------

library(SpadeR)
data(DiversityData)
DiversityData$Inci
# The single-column data include the observed incidence-based frequencies of tropical rainforest ants 
# collected by Berlese extraction of soil samples (217 sampling units) in Costa Rica (Longino et al. 2002).
# In the input data, the entry in the first row denotes the number of sampling units (217); the
# subsequent 117 rows denote species incidence frequencies of the observed species.

DiversityData$Inci_count
# The seed-bank data consist of the observed species incidence-based frequency counts of seedlings
# that germinated from soil samples (Butler and Chazdon, 1998); here each soil sample is regarded
# as a sampling unit. A total of 34 species of seedlings were found in the 121 soil samples. 
# The incidence frequency counts are read as = (121 1 3 2 2 3 3 ... 61 1); each number needs to be separated
# by at least one blank space or by separated by rows. The first entry, indicating that there are 121
# soil samples, is followed by the 18 pairs (1, 3), (2, 2), (3, 3), (4, 3), (5, 1), (6, 5), and so on, up to
# (61, 1). Here (1, 3) indicates that there are 3 unique species, (2, 2) indicates there are 2 duplicate
# species, and so on, with (61, 1) indicating that there is one species found in 61 soil samples.


DiversityData$Inci_raw
# The data consist of raw incidence data of the seed-bank records, described above for the incidence
# frequency counts data. The input data include a 34 x 121 (species-by-sampling-unit) matrix. For
# each element of the matrix, "1" means a detection and "0" means a non-detection.


library(SpadeR)
data(SimilarityMultData)
SimilarityMultData$Inci_raw
dim(SimilarityMultData$Inci_raw)


library(SpadeR)
data(ChaoSpeciesData)
?ChaoSpeciesData
ChaoSpeciesData$Inci_raw


data(ChaoSharedData)
?ChaoSharedData
ChaoSharedData$Inci_raw

data(SimilarityPairData)
SimilarityPairData$Inci_raw



# vegan - Abundance data --------------------------------------------------
library(vegan)
data(dune)
data(BCI)


# iNEXT - Abundance data --------------------------------------------------
library(iNEXT)
data(bird)
data(spider)
?spider


# unmarked data -----------------------------------------------------------
library(unmarked)
data(package = "unmarked")
data(crossbill)
?crossbill
data(cruz)
?cruz
data(frogs)
data(birds)
?birds
data(gf)
?gf


library(unmarked)
data(salamanders)
?salamanders


# mcemGLM -----------------------------------------------------------------
library(mcemGLM)
data("salamander")
?salamander
salamander



# DAAG --------------------------------------------------------------------
library(DAAG)
data("frogs")
?frogs
frogs

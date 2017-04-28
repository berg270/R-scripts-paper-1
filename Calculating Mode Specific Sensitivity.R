########## Calculating Mode Specific Sensitivity (MSS)


##### Set-up
rm(list = ls()) # clear environment
## Set general working directory
setwd("D:/UserData/berg270/My Documents/PhD/ESMU project/Sensitivity and ecoregions paper/Sensitivity/Implementation in R//MSS") 


##### Load packages
library(xlsx)
library(dplyr)


##### Load required databases, check if they are fully loaded (with release notes ECOTOX)!!
MOA <- read.csv('MOAtox.csv')
## For large ECOTOX files, use read.table for initial run of script, or when a new version of ECOTOX is downloaded
#tests <- read.table("D:/UserData/berg270/My Documents/PhD/ESMU project/Sensitivity and ecoregions paper/Sensitivity/ecotox_ascii_12_15_2016/tests.txt", sep = '|', fill = TRUE, quote = "", header = TRUE)
#save(tests, file = 'tests.Rda') # For repetition of the script, save ecotox files as Rdata
load('tests.Rda') # And quickly load ECOTOX file
#results <- read.table("D:/UserData/berg270/My Documents/PhD/ESMU project/Sensitivity and ecoregions paper/Sensitivity/ecotox_ascii_12_15_2016/results.txt", sep = '|', fill = TRUE, quote = "", header = TRUE)
#save(results, file = 'results.Rda')
load('results.Rda')
#species <- read.delim("D:/UserData/berg270/My Documents/PhD/ESMU project/Sensitivity and ecoregions paper/Sensitivity/ecotox_ascii_12_15_2016/validation/species.txt", sep = '|')
#save(species, file = 'species.Rda')
load('species.Rda')


##### MOA
## Remove dashes and spaces from CAS numbers MOAtox to make them comparable with ECOTOX
MOA$CAS <- gsub('-','', MOA$CAS)
## Filter on broad and/or specific MOA
MOA <- filter(MOA, Broad.Mode.of.Action == "Narcosis" )
MOA <- filter(MOA, Specific.Mode.of.Action != "NA") # remove CASnrs for which no specific mode of action is known


##### Extract all test results from ECOTOX with the correct MOA; add species information
names(MOA)[1] <- "test_cas" # change column name so that it can be merged
ECOTOX <- merge(MOA, tests, by = 'test_cas')
ECOTOX <- merge(ECOTOX, results, by = "test_id")
ECOTOX <- merge(ECOTOX, species, by = "species_number") 
#count(distinct(ECOTOX, CAS)) # count number of CAS numbers used in analysis


##### Selection of tox data
ECOTOX <- filter(ECOTOX, kingdom == 'Animalia')
ECOTOX <- filter(ECOTOX, subphylum_div != 'Vertebrata')
ECOTOX <- filter(ECOTOX, organism_habitat == 'Water')
ECOTOX <- filter(ECOTOX, media_type == 'FW')
ECOTOX <- filter(ECOTOX, test_location == 'LAB')
ECOTOX <- filter(ECOTOX, conc1_mean != "NA")
ECOTOX <- filter(ECOTOX, conc1_mean != "NR")
### These parts could be replaced with loops:
## Filter observed test duration
ECOTOX1 <- filter(ECOTOX, obs_duration_mean == 1 & obs_duration_unit == 'd')
ECOTOX2 <- filter(ECOTOX, obs_duration_mean == 2 & obs_duration_unit == 'd')
ECOTOX3 <- filter(ECOTOX, obs_duration_mean == 3 & obs_duration_unit == 'd')
ECOTOX4 <- filter(ECOTOX, obs_duration_mean == 4 & obs_duration_unit == 'd')
ECOTOX5 <- filter(ECOTOX, obs_duration_mean == 24 & obs_duration_unit == 'h')
ECOTOX6 <- filter(ECOTOX, obs_duration_mean == 48 & obs_duration_unit == 'h')
ECOTOX7 <- filter(ECOTOX, obs_duration_mean == 72 & obs_duration_unit == 'h')
ECOTOX8 <- filter(ECOTOX, obs_duration_mean == 96 & obs_duration_unit == 'h')
## Combine selected data
ECOTOX <- rbind(ECOTOX1,ECOTOX2)
ECOTOX <- rbind(ECOTOX, ECOTOX3)
ECOTOX <- rbind(ECOTOX, ECOTOX4)
ECOTOX <- rbind(ECOTOX, ECOTOX5)
ECOTOX <- rbind(ECOTOX, ECOTOX6)
ECOTOX <- rbind(ECOTOX, ECOTOX7)
ECOTOX <- rbind(ECOTOX, ECOTOX8)
## Filter observed endpoints and effect
ECOTOX1 <- filter(ECOTOX, endpoint == 'EC50' & effect == 'MOR')
ECOTOX2 <- filter(ECOTOX, endpoint == 'EC50*' & effect == 'MOR')
ECOTOX3 <- filter(ECOTOX, endpoint == 'EC50' & effect == '~MOR')
ECOTOX4 <- filter(ECOTOX, endpoint == 'EC50*' & effect == '~MOR')
ECOTOX5 <- filter(ECOTOX, endpoint == 'LC50')
ECOTOX6 <- filter(ECOTOX, endpoint == 'LC50*')
ECOTOX7 <- filter(ECOTOX, endpoint == 'LC50/')
## Combine selected data
ECOTOX <- rbind(ECOTOX1,ECOTOX2)
ECOTOX <- rbind(ECOTOX, ECOTOX3)
ECOTOX <- rbind(ECOTOX, ECOTOX4)
ECOTOX <- rbind(ECOTOX, ECOTOX5)
ECOTOX <- rbind(ECOTOX, ECOTOX6)
ECOTOX <- rbind(ECOTOX, ECOTOX7)


##### Correct units of concentration
## Import chemical properties from table Andreas (database of Dick de Zwart?)
Chemical_properties <- read.csv("Chemical_properties_Table_Andreas.csv")
names(Chemical_properties)[1] <- "test_cas"
#names(Chemical_properties)[2] <- "CAS."
ECOTOX <- merge(ECOTOX, Chemical_properties, by = "test_cas", all.x = TRUE)
## Check which units are in the database, think about if and how transformation is possible
units <- unique(ECOTOX$conc1_unit) # to see which units are present
ECOTOX$conc1_mean <- as.numeric(ECOTOX$conc1_mean) # Change conc1_mean to numeric to be able to preform calculations
## Unit transformation
ECOTOX1 <- filter(ECOTOX, conc1_unit == "ug/L" | conc1_unit == "AI ug/L") # (AI:active ingredient) 
ECOTOX1$Standardized <- ECOTOX1$conc1_mean
ECOTOX2 <- filter(ECOTOX, conc1_unit == "mg/L" | conc1_unit == "AI mg/L" | conc1_unit == "mg/dm3" | conc1_unit == "ae mg/L" | conc1_unit == "ug/ml" | conc1_unit == "AI ug/ml") # * 1.000 (ae:acid equivalents)
ECOTOX2$Standardized <- log(ECOTOX2$conc1_mean*1000)
ECOTOX3 <- filter(ECOTOX, conc1_unit == "mg/ml" | conc1_unit == "g/L") # * 1.000.000
ECOTOX3$Standardized <- log(ECOTOX3$conc1_mean * 1000000)
ECOTOX4 <- filter(ECOTOX, conc1_unit == "umol/L" | conc1_unit == "uM" | conc1_unit == "nmol/ml") # * MW (molair weight:g/mol)
ECOTOX4$Standardized <- log(ECOTOX4$conc1_mean * ECOTOX4$MW.g.Mol.)
ECOTOX5 <- filter(ECOTOX, conc1_unit == "mmol/L" | conc1_unit == "mM") # * MW * 1.000
ECOTOX5$Standardized <- log(ECOTOX5$conc1_mean * ECOTOX5$MW.g.Mol. * 1000)
ECOTOX6 <- filter(ECOTOX, conc1_unit == "nM" | conc1_unit == "nmol/L") # (nmol/L) * MW / 1.000
ECOTOX6$Standardized <- log(ECOTOX6$conc1_mean * ECOTOX6$MW.g.Mol. / 1000)
ECOTOX7 <- filter(ECOTOX, conc1_unit == "M" | conc1_unit == "mol/L") # (M:mol/L) * MW * 1.000.000
ECOTOX7$Standardized <- log(ECOTOX7$conc1_mean * ECOTOX7$MW.g.Mol. * 1000000)
## These units can be used if data on density is available (41 obs)
#ECOTOX8 <- filter(ECOTOX, conc1_unit == "ul/L") # * p (density g/L)  
#ECOTOX8$Standardized <- ECOTOX8$MW.g.Mol..x * density
#ECOTOX9 <- filter(ECOTOX, conc1_unit == "AI ml/L") # * p * 1.000.000 
#ECOTOX9$Standardized <- ECOTOX9$conc1_mean * density
#ECOTOX10 <- filter(ECOTOX, conc1_unit == "ml/L") # * p * 1.000
#ECOTOX10$Standardized <- ECOTOX10$conc1_mean * 1000 * densities
## These units are not usable, because cannot be transformed to ug/L (1154 obs)
#ECOTOX11 <- filter(ECOTOX, conc1_unit == "ppm") # mass/volume/mole fraction? --> not usable
#ECOTOX12 <- filter(ECOTOX, conc1_unit == "ppb") # mass/volume/mole fraction? --> not usable
#ECOTOX13 <- filter(ECOTOX, conc1_unit == "L/ha") # not usable?
#ECOTOX14 <- filter(ECOTOX, conc1_unit == "%") # percentage of what? --> not usable
#ECOTOX15 <- filter(ECOTOX, conc1_unit == "AI ppm") # mass/volume/mole fraction? --> not usable
#ECOTOX16 <- filter(ECOTOX, conc1_unit == "ppt") # mass/volume/mole fraction? --> not usable
#ECOTOX17 <- filter(ECOTOX, conc1_unit == "% v/v") # (percent volume per volume) --> what are the units? --> not usable
#ECOTOX18 <- filter(ECOTOX, conc1_unit == "ug/kg/d") # impossible to convert to ug/L
## Combine selected data
ECOTOX <- rbind(ECOTOX1,ECOTOX2)
ECOTOX <- rbind(ECOTOX, ECOTOX3)
ECOTOX <- rbind(ECOTOX, ECOTOX4)
ECOTOX <- rbind(ECOTOX, ECOTOX5)
ECOTOX <- rbind(ECOTOX, ECOTOX6)
ECOTOX <- rbind(ECOTOX, ECOTOX7)


##### Calculate MSS values
## Create list of unique CAS numbers 
uniqueCAS <- unique(ECOTOX$CAS)
## Select only CASnumbers for which > 3 species have been tested
cnt <- list() # Can I turn this into a dataframe? ***
subs <- select(ECOTOX, CAS, species_number) # select only necessary columns
for(i in 1:length(uniqueCAS)) {
  aux <- subset(subs, subs$CAS == uniqueCAS[i]) # aux (temporarily) contains all the rows of one CASnr
  cnt[[i]] <- length(unique(aux$species_number)) > 3 # cnt registers if for each CASnr more (TRUE) or fewer (FALSE) than 3 species have been tested
}
selectedCASnrs <- as.data.frame(uniqueCAS) # *** Because then I can lose these two sentences
selectedCASnrs$count <- cnt # ***
final_selectedCASnrs <- filter(selectedCASnrs, count == 'TRUE') # select only those CASnrs for which > 3 species have been tested
## Within each CASnr, calculate the relative sensitivity (RS) of each species
relative_sensitivity <- data.frame()
for(i in 1:NROW(final_selectedCASnrs)){
  aux = subset(ECOTOX, ECOTOX$CAS == final_selectedCASnrs[i,1])
  aux2 <- aggregate(aux$conc1_mean, by = list(aux$species_number), FUN = mean) # Calculate average if more than 1 entry per species-chemical combination is present
  aux2$RS <- (aux2$x - mean(aux2$x)) / sd(aux2$x) # Calculate relative sensitivity (RS)       
  relative_sensitivity <- rbind(relative_sensitivity, aux2)  
}
## Calculate MSS
MSS <- aggregate(relative_sensitivity$RS, by = list(relative_sensitivity$Group.1), FUN = mean) # Over all chemicals, calculate the average RS (group1:species_number)
names(MSS)[1] <- "species_number"
MSS_with_taxonomy <- merge(MSS, species, by = "species_number")
write.xlsx(MSS_with_taxonomy, 'MSS_with_taxonomy.xlsx')



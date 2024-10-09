# Installing packages
library(sf)
library(dplyr)
library(mFD)
library(picante)
library(tidyverse)
library(corrplot)
library(caret)
library(missForest)
library(tibble)
library(biscale)
library(ggspatial)
library(cowplot)
library(FD)


# clean environment
rm(list = ls())

# load Spanish UTM grid
malla <- st_read("Spatial_Data/Malla_municipios/Malla10x10_Ter_p.shp")

# load atlas data (2014-2018)
atlas16 <- read.csv("Data/SteppedBirdLIst Atlas 2014_­2018.csv", stringsAsFactors = F)

#Converting atlas16 into a spatial layer
Atlas16_shp <- left_join(malla[, c("UTMCODE", "CUADRICULA", "XCENTROIDE", "YCENTROIDE")], atlas16, 
                     by = c("UTMCODE"="UTM"))

# Loading Ebird data 2019-2023 (26 csv files) - Peninsular Spain
file_list <- list.files(path="Data_Ebird/", pattern = "\\.csv$", full.names = TRUE)
# Loading each file in a separate dataframe
data_list <- lapply(file_list, read.csv)
# There are some Ebird records reported as X, which are presence observations, therefore, I am assigning a 1 to those rows 
data_list <- lapply(data_list, function(df) {
  df$max_indivs <- ifelse(df$max_indivs == "X", 1, df$max_indivs)
  return(df)
})
# Some values in "max_indivs" are characters, therefore transforming into numeric. 
data_list <- lapply(data_list, function(df) {
  df$max_indivs <- as.numeric(df$max_indivs)
  return(df)
})
#Merging all 26 datasets 
Steppebirds_Ebird <- bind_rows(data_list)

# We are only using breeding season data, therefore selecting only records from february to August 
Steppebirds_EbirdFebAug <- filter(Steppebirds_Ebird, month %in% c("2", "3", "4", "5", "6", "7", "8")) 

# I only need to know if species were present or not in each UTM, therefore, assigning a 0 for absence and a 1 for presence. 
Steppebirds_EbirdFebAug$presence<- ifelse(Steppebirds_EbirdFebAug$max_indivs > 0, 1, 0)
Steppebirds_present <- Steppebirds_EbirdFebAug[Steppebirds_EbirdFebAug$presence == 1, ]

Atlas16  <- as.data.frame(Atlas16_shp)
names (Atlas16)[2] = "UTM"
names (Atlas16)[9] = "species"

Atlas16_1=select(Atlas16,UTM,species)
Steppebirds_pres=select(Steppebirds_present,UTM, species)

#Changing Curruca conspicillata to Sylvia conspicillata
Steppebirds_pres$species <- gsub("Curruca conspicillata", "Sylvia conspicillata", Steppebirds_pres$species)

#Merging data from Atlas 16 and Ebird
Atlas_Ebird <- rbind(Atlas16_1, Steppebirds_pres)

#Removing duplicates based on UTM and species
Communities1 <- distinct(Atlas_Ebird, UTM, species)

#Removing NAs
Communities1 <- na.omit(Communities1)

# Number of species per UTM
species_count <- aggregate(Communities1$species, by=list(Communities1$UTM), FUN=length)

# Selecting only UTM cells with 3 or more species. This is done because the number of species per UTM 
# must be higher than the number of axes to compute the convex hull in FD analyses. I am using 2 axes,therefore 
# FD indices are not calculated in UTM cells with less than 3 species.
Communities1_1 <- Communities1[Communities1$UTM %in% species_count[species_count$x >= 5,]$Group.1,]

# Transposing species column to become names of columns, and assigning a 1 or 0 for presence or absence in each cell.
Communities2 <- Communities1_1 %>%
  mutate(presence = 1) %>%
  spread(key = species, value = presence, fill = 0) 

# UTM column as row names 
rownames(Communities2) <- Communities2$UTM
Communities2$UTM <- NULL 

Communities2 <- as.matrix(Communities2)

# Load steppe birds trait data:
Species_traits <- read.csv("Data/TraitsSteppebirdsOct_2023.csv", stringsAsFactors = F)
str(Species_traits)
names(Species_traits)

# Removing columns related to common names, family, order, references and some traits that wont be worth to analyse 
# such as: Trophic level (Trophic niche is more specific), Maximum longevity (most values are the same as longevity), 
# relative_brain_size (removing now to calculate after imputation of brain size NAs values)
Species_traits <- Species_traits[ -c(2:4,19,23,25,30,35,36,38,41,42,46) ]

Species_traits$Habitat <- as.factor(Species_traits$Habitat)
Species_traits$Habitat.Density <- as.factor(Species_traits$Habitat.Density)
Species_traits$Migration <- as.factor(Species_traits$Migration)
Species_traits$Trophic.Niche <- as.factor(Species_traits$Trophic.Niche)
Species_traits$Primary.Lifestyle <- as.factor(Species_traits$Primary.Lifestyle)
Species_traits$Degree_of_development <- as.factor(Species_traits$Degree_of_development)
Species_traits$Gregariousness <- as.factor(Species_traits$Gregariousness)
Species_traits$Mating_system <- as.factor(Species_traits$Mating_system)

# Counting the number of NAs per each variable
nasTraits = colSums(is.na(Species_traits))

nasTraits
# Just a few NAs: PopDensity_ind_km2(3), birth_or_hatching_weight_g(10), brain size(8), hab_breadth(2), 
# degree of development(1)               

# Imputation for NAs
set.seed(42)
imp <- missForest(Species_traits[,c(2:33)], verbose=FALSE)
Species_traits_imp <- cbind(Species_traits$ESP_LAT, imp$ximp)

names (Species_traits_imp)[1] = "ESP_LAT"

# Now calculating relative brain size to body mass using residuals. Positive values relate to brains that are
# bigger than expected
model <- lm(brain_size_g ~ Mass, data = Species_traits_imp)
residuals  <-  model$residuals
Species_traits_imp$relative_brain <- residuals
Species_traits_imp$brain_size_g  <- NULL

#Loading metadata: name of each variable and type of variable: Q (quantitative), N (nominal)
traits_cat <- read.csv("Data/Traits_categories.csv", stringsAsFactors = F)

#Setting species names as row names
rownames(Species_traits_imp) <- Species_traits_imp$ESP_LAT
Species_traits_imp$ESP_LAT <- NULL 

# Species traits summary:
species_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = traits_cat,   
  sp_tr      = Species_traits_imp, 
  stop_if_NA = TRUE)

# Summary of the assemblages * species dataframe:
asb_communities_summ <- mFD::asb.sp.summary(asb_sp_w = Communities2)
asb_sp_occ <- asb_communities_summ$"asb_sp_occ"
asb_communities_summ$"sp_tot_w"              # Species total presences in all assemblages
asb_communities_summ$"asb_sp_richn"           # Species richness per community
asb_communities_summ$"asb_sp_nm"[[1]]             # Names of species present in the first assemblage

#Computing distances between species based on functional traits

sp_dist         <- mFD::funct.dist(
  sp_tr         = Species_traits_imp,
  tr_cat        = traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

summary(as.matrix(sp_dist))
# This function returns a dist object with traits-based distances between all pairs of species
round(sp_dist, 3)   
traits_based_distances <- as.matrix(round(sp_dist, 3)  )

#Computing functional spaces & their quality
fspaces_quality      <- mFD::quality.fspaces(
  sp_dist             = sp_dist,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = TRUE,
  fdendro             = "average")

# Quality metrics of spaces
round(fspaces_quality$"quality_fspaces", 3) 

sp_faxes_coord_birds <- fspaces_quality$details_fspaces$sp_pc_coord

fspaces_quality$details_fspaces$pc_eigenvalues

#Calculating FD indexes
alpha_fd_indices <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_birds[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = Communities2,
  ind_vect         = c("fdis", "fmpd", "fnnd", "feve", "fric", "fdiv", "fori", 
                       "fspe", "fide"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values <- alpha_fd_indices$"functional_diversity_indices"

# Converting row names (UTM) into a column:
FDindexes <- rownames_to_column(fd_ind_values, "UTMCODE")

#Converting Index values into a spatial layer
FDindex <- left_join(malla[, c("UTMCODE", "CUADRICULA", "XCENTROIDE", "YCENTROIDE")],  FDindexes, 
                      by = c("CUADRICULA"="UTMCODE"))
sum(is.na(FDindex$fric))

# Saving FD indexes as spatial data
write_sf(FDindex, "Spatial_Data/Richness/FuncDiversity.shp")

# Including Spain provinces
provinces <- mapSpain::esp_get_prov()
provinces <- provinces[!provinces$iso2.prov.name.es %in% c("Santa Cruz de Tenerife", "Las Palmas"),]

mapTD_FD1 <- mapTD_FD + 
  geom_sf(data = provinces, fill = NA, color = "black", linewidth = 0.5) +
  annotation_north_arrow(location = "tr", which_north = "grid", 
                         pad_x = unit(0, "in"), pad_y = unit(0.05, "in"),
                         style = north_arrow_fancy_orienteering())

legendTDFD <- bi_legend(pal = "DkBlue2",
                      dim = 3,
                      xlab = "TD ",
                      ylab = "FD ",
                      size = 13)

TDFDmap <- ggdraw() +
  draw_plot(mapTD_FD1, 0, 0, 0.8, 1) +
  draw_plot(legendTDFD, 0.8, 0, 0.2, 1) 

#Exporting figure
ggsave("Figures/TDvsFDmap.png", TDFDmap, wi = 20, he = 20, un = "cm", dpi = 300)


###################################
###  Calculating FRic WITH FD PACKAGE
###################################

FDindexes_FD <- dbFD(
  Species_traits_imp, Communities2, w.abun=FALSE, 
  corr="sqrt", calc.FRic=TRUE, m=4, 
  stand.FRic=TRUE, scale.RaoQ=FALSE, calc.FGR=FALSE,
  calc.CWM=FALSE, print.pco=FALSE
)

FDindexes_FD$x.values
FD_package  <- as.data.frame(FDindexes_FD)

unique(length(FDindexes_FD$FRic))
###################################
###  Null model permutations to calculate SESFRic######
###################################

set.seed(1)
null.model.type <- "independentswap" # available methods: "independentswap", "frequency", "richness", "trialswap"
nperm <- 1000
FRic <- matrix(NaN, nrow=nperm, ncol=length(FDindexes_FD$FRic))
abund.null <- vector(mode="list", nperm)
t1 <- Sys.time()
for(i in seq(nperm)){ # loop for each permutation
  # Create permutated abundance matrix
  abund.i <- randomizeMatrix(Communities2, null.model = null.model.type)
  traits.i <- Species_traits_imp
  
  abund.null[[i]] <- abund.i
  
  # Record results
  res.i <- dbFD(
    traits.i, abund.i, w.abun=FALSE, 
    corr="sqrt", calc.FRic=TRUE, m=4,
    stand.FRic=TRUE, calc.FGR=FALSE,
    calc.CWM=FALSE, print.pco=FALSE
  )
  # Save result
  FRic[i,] <- res.i$FRic
  
  rm(abund.i); rm(res.i)
  print(paste0(round(i/nperm*100), "% of permutations completed"))
}
t2 <- Sys.time()
t2 - t1 # time elapsed
(t2 - t1)/nperm * 1000 # time required for 1000 permutations

### Summary stats
# FRic
SESFRic = (FDindexes_FD$FRic - apply(FRic, 2, mean)) / apply(FRic, 2, sd) # standardized effect size
qFRic <- NaN*FDindexes_FD$FRic
for(i in seq(qFRic)){
  qFRic[i] <- sum(FDindexes_FD$FRic[i] > FRic[,i]) / length(FRic[,i])
}
sigFRic <- qFRic < 0.05 | qFRic > 0.95 # test if outside distribution

FRic_index <- FD_package
SESFRic_index <- as.data.frame(SESFRic)

# Converting row names (UTM) into a column:
SESFRic_index <- rownames_to_column(SESFRic_index, "UTMCODE")

# Combine SESFRic and FRic into one dataframe
final_index <- cbind(SESFRic_index, FRic_index)

cor(final_index$SESFRic, final_index$FRic)

#Converting Index values into a spatial layer
FDindex_SES <- left_join(malla[, c("UTMCODE", "CUADRICULA", "XCENTROIDE", "YCENTROIDE")],  final_index, 
                     by = c("CUADRICULA"="UTMCODE"))

# Saving FD indexes as spatial data
write_sf(FDindex_SES, "Spatial_Data/Richness/FDindex_SES.shp")

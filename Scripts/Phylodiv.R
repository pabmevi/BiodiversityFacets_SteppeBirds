#Script to analyse phylogenetic diversity

library(ape)
library(phytools)
library(sf)
library(dplyr)
library(picante)
library(tidyr)
library(tibble)
library(rgl)
library(ggplot2)
library(ggspatial)
library(tmap)
library(mapSpain)
library(units)
library(viridis)

# clean environment
rm(list = ls())
setwd("~/GitHub/BiodiversityFacets")

tree <- read.nexus("Data/consensustree150_05credibility.nex")

# Species to keep
species <- c("Circus_pygargus", "Circus_cyaneus", "Burhinus_oedicnemus", "Alauda_arvensis", "Chersophilus_duponti", 
             "Otis_tarda", "Anthus_campestris", "Asio_flammeus", "Cisticola_juncidis", "Melanocorypha_calandra", 
             "Glareola_pratincola", "Falco_naumanni", "Coturnix_coturnix", "Galerida_cristata", "Galerida_theklae", 
             "Oenanthe_oenanthe", "Oenanthe_hispanica", "Cursorius_cursor", "Sylvia_conspicillata", "Pterocles_alchata", 
             "Pterocles_orientalis", "Alectoris_rufa", "Tetrax_tetrax", "Calandrella_brachydactyla", "Calandrella_rufescens", 
             "Miliaria_calandra")

tree <- keep.tip(tree, species)
plot(tree, show.tip.label = TRUE, cex = 0.8)
axisPhylo()

# load Spanish UTM grid
malla <- st_read("Spatial_Data/Malla_municipios/Malla10x10_Ter_p.shp") #peninsular

# load atlas data (2014-2018)
atlas16 <- read.csv("Data/Atlas_CP.csv", stringsAsFactors = F)

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

Atlas16_1=select(Atlas16, UTM, species)
Steppebirds_pres=select(Steppebirds_present,UTM, species)

#Changing Curruca conspicillata to Sylvia conspicillata
Steppebirds_pres$species <- gsub("Curruca conspicillata", "Sylvia conspicillata", Steppebirds_pres$species)
Steppebirds_pres$species <- gsub("Emberiza calandra", "Miliaria calandra", Steppebirds_pres$species)
Steppebirds_pres$species <- gsub("Alaudala rufescens", "Calandrella rufescens", Steppebirds_pres$species)

Atlas16_1$species <- gsub("Curruca conspicillata", "Sylvia conspicillata", Atlas16_1$species)
Atlas16_1$species <- gsub("Emberiza calandra", "Miliaria calandra", Atlas16_1$species)
Atlas16_1$species <- gsub("Alaudala rufescens", "Calandrella rufescens", Atlas16_1$species)

#Merging data from Atlas 16 and Ebird
Atlas_Ebird <- rbind(Atlas16_1, Steppebirds_pres)

#Removing duplicates based on UTM and species
Communities1 <- distinct(Atlas_Ebird, UTM, species)

#Removing NAs
Communities1 <- na.omit(Communities1)

Communities1$species <- gsub(" ", "_", Communities1$species)

# Number of species per UTM
species_count <- aggregate(Communities1$species, by=list(Communities1$UTM), FUN=length)

## Selecting only UTM cells with 5 or more species. which we define as minimal habitat area for steppe birds.
Communities1_0 <- Communities1[Communities1$UTM %in% species_count[species_count$x >= 5,]$Group.1,]

# Transposing species column to become names of columns, and assigning a 1 or 0 for presence or absence in each cell.
Communities2 <- Communities1_0 %>%
  mutate(presence = 1) %>%
  spread(key = species, value = presence, fill = 0) 

# UTM column as row names 
rownames(Communities2) <- Communities2$UTM
Communities2$UTM <- NULL 

Communities2 <- as.data.frame(Communities2)
#Phylogenetic diversity per cell
phylodiv <- pd(Communities2, tree, include.root = TRUE)

# Converting row names (UTM) into a column:
phylodiv <- rownames_to_column(phylodiv, "UTMCODE")

# PD index can be highly influenced by TD, then I am calculating the effect size of PD to avoid biases towards TD
set.seed(123)
phylodiv_null <- ses.pd(Communities2, tree, null.model = "taxa.labels", runs = 1000, include.root = TRUE)

# Converting row names (UTM) into a column:
phylodiv_null <- rownames_to_column(phylodiv_null, "UTMCODE")


#Converting Index values into a spatial layer
phylodiv1 <- left_join(malla[, c("UTMCODE", "CUADRICULA")],  phylodiv_null, 
                     by = c("CUADRICULA"="UTMCODE"))

#Loading the 3 diversity facets
TD <- st_read("Spatial_Data/Richness/TDindex.shp")
FD <- st_read("Spatial_Data/Richness/FDindex_SES.shp")
PD <- phylodiv1

FD$geometry  <- NULL
PD$geometry  <- NULL

TD <- TD %>%
  filter(species >= 5)

FD <- FD %>%
  filter(nbsp >= 5)

PD <- PD %>%
  filter(ntaxa >= 5)

datos <- merge(TD, FD, by = "UTMCODE")
datosMAP <- merge(datos, PD, by = "UTMCODE")

datosMAP <- datosMAP[, c("UTMCODE","species", "FRic", "SESFRic", "pd.obs.z", "pd.obs")]
datosMAP <- st_transform(datosMAP, 25830) 

datosMAP$area_cell <- st_area(datosMAP)

#Removing units from rows
datosMAP <- drop_units(datosMAP)
names(datosMAP)[5] <- "SESPD"
names(datosMAP)[6] <- "PD"
#Now I want to analyse the correlation of TD, FD, and PD. 

forcorrmatrix <- datosMAP[, c("species", "FRic", "SESFRic", "SESPD", "PD")]
forcorrmatrix$geometry <- NULL
corrmatrix <- cor(forcorrmatrix, use = "complete.obs") #Very high correlation between TD - FD = 81%, and TD - PD = 94%.

write_sf(datosMAP, "Spatial_Data/3facets.shp")

##################################################################
# Calculating the contribution of each species to FRic
##################################################################

# Calculate total PD
total_PD <- sum(tree$edge.length)

# Species list
species_list <- tree$tip.label

contribuciones <- data.frame(
  species = species_list,
  PD_total = NA,
  PD_sin_especie = NA,
  contrib = NA
)

# Calculate each species' contribution
for (i in seq_along(species_list)) {
  sp <- species_list[i]
  tree_sin_sp <- drop.tip(tree, sp)
  PD_sin_sp <- sum(tree_sin_sp$edge.length)
  
  contribuciones$PD_total[i] <- total_PD
  contribuciones$PD_sin_especie[i] <- PD_sin_sp
  contribuciones$contrib[i] <- total_PD - PD_sin_sp
}

# Remove underscores from species names
contribuciones$species_clean <- gsub("_", " ", contribuciones$species)

# Order data by contribution for better visualization
contribuciones_sorted <- contribuciones[order(contribuciones$contrib), ]

phy <- ggplot(contribuciones_sorted, aes(x = reorder(species_clean, contrib), y = contrib)) +
  geom_bar(stat = "identity", fill = "grey") +
  coord_flip() +  
  labs(
    x = "Species",
    y = "PD Contribution",
    title = ""
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    axis.text.y = element_text(family = "Arial", face = "italic", size = 12),
  )

ggsave("Figures/Species contribution to PD.png", phy, wi = 20, he = 20, un = "cm", dpi = 300)

#Plotting the phylogenetic tree
png("Figures/phylotree.png", width = 1000, height = 1000, res = 150, family = "Arial")
plot(tree, show.tip.label = TRUE, cex = 1.2, font = 3)
dev.off()

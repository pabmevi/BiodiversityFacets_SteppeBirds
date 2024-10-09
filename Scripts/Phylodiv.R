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

tree <- read.nexus("Data/consensustree150_05credibility.nex")

# Species to keep
species <- c("Circus_pygargus", "Circus_cyaneus", "Burhinus_oedicnemus", "Alauda_arvensis", "Chersophilus_duponti", 
             "Otis_tarda", "Anthus_campestris", "Asio_flammeus", "Cisticola_juncidis", "Melanocorypha_calandra", 
             "Glareola_pratincola", "Falco_naumanni", "Coturnix_coturnix", "Galerida_cristata", "Galerida_theklae", 
             "Oenanthe_oenanthe", "Oenanthe_hispanica", "Cursorius_cursor", "Sylvia_conspicillata", "Pterocles_alchata", 
             "Pterocles_orientalis", "Alectoris_rufa", "Tetrax_tetrax", "Calandrella_brachydactyla", "Calandrella_rufescens", 
             "Miliaria_calandra")

tree <- keep.tip(tree, species)

# load Spanish UTM grid
malla <- st_read("Spatial_Data/Malla_municipios/Malla10x10_Ter_p.shp") #peninsular

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
phylodiv <- left_join(malla[, c("UTMCODE", "CUADRICULA")],  phylodiv_null, 
                     by = c("CUADRICULA"="UTMCODE"))

# Saving FD indexes as spatial data
write_sf(phylodiv, "Spatial_Data/Richness/PhylDiversity_sizeeffect.shp")

####################################################################
############### Overlapping the 3 indexes############################
TD <- st_read("Spatial_Data/Richness/At_Eb2023_p.shp")
FD <- st_read("Spatial_Data/Richness/FDindex_SES.shp")
PD <- st_read("Spatial_Data/Richness/PhylDiversity_sizeeffect.shp")

# Adding Autonomous communities from Spain
comm <-esp_get_ccaa()
comm <- comm[!comm$iso2.ccaa.name.es %in% c("Canarias"),]

TD$geometry  <- NULL
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

#Converting atlas16 into a spatial layer
datosMAP <- left_join(malla[, c("UTMCODE", "CUADRICULA")], datosMAP, 
                      by = c("UTMCODE"))

datos_seleccionados <- datosMAP[, c("species", "FRic", "SESFRic", "pd_bs_z", "pd_obs")]
datos_seleccionados$geometry <- NULL
matriz_correlacion <- cor(datos_seleccionados, use = "complete.obs")

#Very high correlation between TD - FD = 80%, and TD - PD = 93%.

Fric <- ggplot() +
  geom_sf(data = datosMAP, aes(fill = FRic)) +
  geom_sf(data = comunidades_autonomas, fill = "transparent", color = "black") +
  scale_fill_viridis(option = "cividis", na.value = "grey95", direction = -1 )+
  labs(x = "Longitude", y = "Latitude", fill = "Fric" ) +
  annotation_scale(location = "bl", width_hint = .3) +
  theme_bw() +
  theme(title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.position = c(0.85, 0.15),
        axis.title = element_text(size = 15, face = "plain")) + 
  guides(colour = "none", size = "none", alpha = "none")

ggsave("Figures/FDindex.png", Fric, wi = 45, he = 20, un = "cm", dpi = 300)

SESFric <- ggplot() +
  geom_sf(data = datosMAP, aes(fill = SESFRic)) +
  scale_fill_viridis(option = "cividis", na.value = "grey95", direction = -1 )+
  labs(x = "Longitude", y = "Latitude", fill = "SESFric" ) +
  annotation_scale(location = "bl", width_hint = .3) +
  theme_bw() +
  theme(title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.position = c(0.85, 0.15),
        axis.title = element_text(size = 15, face = "plain")) + 
  guides(colour = "none", size = "none", alpha = "none")

ggsave("Figures/SES_FD.png", SESFric, wi = 45, he = 20, un = "cm", dpi = 300)

Richn <- ggplot() +
  geom_sf(data = datosMAP, aes(fill = species)) +
  geom_sf(data = comunidades_autonomas, fill = "transparent", color = "black") +
  scale_fill_viridis(option = "cividis", na.value = "grey95", direction = -1 )+
  labs(x = "Longitude", y = "Latitude", fill = "Richness" ) +
  annotation_scale(location = "bl", width_hint = .3) +
  theme_bw() +
  theme(title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.position = c(0.85, 0.15),
        axis.title = element_text(size = 15, face = "plain")) + 
  guides(colour = "none", size = "none", alpha = "none")

ggsave("Figures/TDindex.png", Richn, wi = 45, he = 20, un = "cm", dpi = 300)

Phylo <- ggplot() +
  geom_sf(data = datosMAP, aes(fill = pd_obs)) +
  geom_sf(data = comunidades_autonomas, fill = "transparent", color = "black") +
  scale_fill_viridis(option = "cividis", na.value = "grey95", direction = -1 )+
  labs(x = "Longitude", y = "Latitude", fill = "PD" ) +
  annotation_scale(location = "bl", width_hint = .3) +
  theme_bw() +
  theme(title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.position = c(0.85, 0.15),
        axis.title = element_text(size = 15, face = "plain")) + 
  guides(colour = "none", size = "none", alpha = "none")

ggsave("Figures/PDindex.png", Phylo, wi = 45, he = 20, un = "cm", dpi = 300)

SES_PD <- ggplot() +
  geom_sf(data = datosMAP, aes(fill = pd_bs_z)) + #pd_bs_z corresponds to the effect size
  geom_sf(data = comunidades_autonomas, fill = "transparent", color = "black") +
  scale_fill_viridis(option = "cividis", na.value = "grey95", direction = -1 )+
  labs(x = "Longitude", y = "Latitude", fill = "SESPD" ) +
  annotation_scale(location = "bl", width_hint = .3) +
  theme_bw() +
  theme(title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.position = c(0.85, 0.15),
        axis.title = element_text(size = 15, face = "plain")) + 
  guides(colour = "none", size = "none", alpha = "none")

ggsave("Figures/SES_PD.png", SES_PD, wi = 45, he = 20, un = "cm", dpi = 300)

# Calculating percentil 70 for each index
cutoff_sp_richn <- quantile(datosMAP$species, probs = 0.70, na.rm = TRUE)
cutoff_fric <- quantile(datosMAP$SESFRic, probs = 0.70, na.rm = TRUE)
cutoff_PD <- quantile(datosMAP$pd_bs_z, probs = 0.70, na.rm = TRUE)

# New variable for each index. It is TRUE if the value is inside the top 30%
datosMAP$top_sp_richn <- datosMAP$species > cutoff_sp_richn
datosMAP$top_SESfric <- datosMAP$SESFRic > cutoff_fric
datosMAP$top_SESPD <- datosMAP$pd_bs_z > cutoff_PD

# Creating a new variable to represent overlapping
datosMAP$overlap <- with(datosMAP, ifelse(top_sp_richn & top_SESfric & top_SESPD, "Hotspot",
                                          ifelse(top_sp_richn & top_SESfric, "Top 30% sp_richn and SES_fric",
                                                 ifelse(top_sp_richn & top_SESPD, "Top 30% sp_richn and SES_PD",
                                                        ifelse(top_SESfric & top_SESPD, "Top 30% SES_fric and SES_PD",
                                                               ifelse(top_sp_richn, "Top 30% sp_richn",
                                                                      ifelse(top_SESfric, "Top 30% SES_fric",
                                                                             ifelse(top_SESPD, "Top 30% SES_PD", "Outside hotspots"))))))))

# Colours for each overlapping category 
colores <- c("Hotspot" = "red", "Top 30% sp_richn and SES_fric" = "blue", "Top 30% sp_richn and SES_PD" = "green", 
             "Top 30% SES_fric and SES_PD" = "yellow", "Top 30% sp_richn" = "purple", "Top 30% SES_fric" = "lightblue",
             "Top 30% SES_PD" = "cyan", "Outside hotspots" = "white")

# Calculating percentil 70 for each index
cutoff_sp_richn <- quantile(datosMAP$species, probs = 0.83, na.rm = TRUE)
cutoff_fric <- quantile(datosMAP$SESFRic, probs = 0.83, na.rm = TRUE)
cutoff_PD <- quantile(datosMAP$pd_bs_z, probs = 0.83, na.rm = TRUE)

# New variable for each index. It is TRUE if the value is inside the top 17%
datosMAP$top_sp_richn <- datosMAP$species > cutoff_sp_richn
datosMAP$top_SESfric <- datosMAP$SESFRic > cutoff_fric
datosMAP$top_SESPD <- datosMAP$pd_bs_z > cutoff_PD

# Creating a new variable to represent overlapping
datosMAP$overlap <- with(datosMAP, ifelse(top_sp_richn & top_SESfric & top_SESPD, "Hotspot",
                                          ifelse(top_sp_richn & top_SESfric, "Top 17% sp_richn and SES_fric",
                                                 ifelse(top_sp_richn & top_SESPD, "Top 17% sp_richn and SES_PD",
                                                        ifelse(top_SESfric & top_SESPD, "Top 17% SES_fric and SES_PD",
                                                               ifelse(top_sp_richn, "Top 17% sp_richn",
                                                                      ifelse(top_SESfric, "Top 17% SES_fric",
                                                                             ifelse(top_SESPD, "Top 17% SES_PD", "Outside hotspots"))))))))

# Colours for each overlapping category 
colores <- c("Hotspot" = "red", "Top 17% sp_richn and SES_fric" = "blue", "Top 17% sp_richn and SES_PD" = "green", 
             "Top 17% SES_fric and SES_PD" = "yellow", "Top 17% sp_richn" = "purple", "Top 17% SES_fric" = "lightblue",
             "Top 17% SES_PD" = "cyan", "Outside hotspots" = "white")

 # Final map
mapa <- tm_shape(datosMAP) +
  tm_fill(col = "overlap", palette = colores, title = "Index overlapping", style = "cat")

# Adding Spain borders to final map
mapa <- tm_shape(datosMAP) +
  tm_fill(col = "overlap", palette = colores, title = "Index overlapping", style = "cat") +
  tm_shape(comunidades_autonomas) +
  tm_borders()  

mapa

# Creating a new variable to represent overlapping
datosMAP$hotspots <- with(datosMAP, ifelse(top_sp_richn & top_SESfric & top_SESPD, "Hotspot",
                                           "Outside hotspots"))

# Asignar "Outside hotspots" a los NA en la columna hotspots
datosMAP$hotspots[is.na(datosMAP$hotspots)] <- "Outside hotspots"

# Crear un nuevo objeto solo con las filas donde hotspots es "Hotspot"
hotspots_data <- subset(datosMAP, hotspots == "Hotspot")

#Exporting figure
tmap_save(mapa, filename = "Figures/hotspot3indxSES_30perc_25june.png", width = 10, height = 10, dpi = 300)

######################CONFLICT PV###################

# load energy infrastructures 

#Loading Atlas_Ebird 2018-2023 shp
Atlas_Ebird23 <- st_read("Spatial_Data/Richness/At_Eb2023_p.shp")
#Loading PV power per municipality
PV <- st_read("Spatial_Data/Energy/PV_p.shp")
#Transforming to have the same coordinates system
PV <- st_transform(PV, 25830)

# Removing provinces with no MW 
PV1 <- PV[complete.cases(PV$MW),]

#Calculating the area of each intersected section
intersect0 <- st_intersection(PV1, Atlas_Ebird23)
intersect0$areaintsct_km2 <- st_area(intersect0)/1000000
#Removing units from rows
intersect0 <- drop_units(intersect0)

#Calculating the proportional MW for each intersection based on the size of the containing cell
intersect0$Proport_MW <- intersect0$MW * intersect0$areaintsct_km2 / intersect0$Area_km2

#Transforming the shp into a data frame 
intersect  <- select(intersect0, UTMCODE, MW, Area_km2, species, areaintsct_km2, Proport_MW)
intersect$geometry <- NULL

#Summing the total PV power per cell 
PV_MW_pcell <- intersect %>% 
  group_by(UTMCODE) %>% 
  summarize(MW_pcell = sum(Proport_MW)) 

#Merging PropMW_pcell and Atlas2014-2018+Ebird 2019-2023
PVpcell <- left_join(Atlas_Ebird23[, c("UTMCODE", "XCENTROIDE", "YCENTROIDE", "Area_km2", "species")], 
                     PV_MW_pcell, by = c("UTMCODE"))

PVpcell1  <-  PVpcell[complete.cases(PVpcell$MW_pcell), ]

# Calcula los cuartiles (25% y 75%)
PVpcellQ <- quantile(PVpcell1$MW_pcell, probs = c(0.25, 0.75))

# Asigna categorías según los cuartiles
PVpcell1$category <- cut(PVpcell1$MW_pcell,
                    breaks = c(-Inf, PVpcellQ[1], PVpcellQ[2], Inf),
                    labels = c("low", "medium", "high"))

# Creating a new variable to represent overlapping
datosMAP$hotspots <- with(datosMAP, ifelse(top_sp_richn & top_SESfric & top_SESPD, "Hotspot",
                                          "Outside hotspots"))

datosMAP1 <- datosMAP %>%
  filter(hotspots == "Hotspot")

# Asegúrate de que las variables estén definidas previamente
# comunidades_autonomas, datosMAP y PVpcell

conflictPV <- ggplot() +
  geom_sf(data = PVpcell1, aes(fill = category), color = NA) + 
  geom_sf(data = datosMAP1, fill = "transparent", color = "red") +
  geom_sf(data = comunidades_autonomas, fill = "transparent", color = "black") +# Utiliza la columna real de tus datos
  scale_fill_manual(values = c("low" = "thistle", "medium" = "lightblue", "high" = "blue")) +  # Asigna colores a cada categoría
  labs(x = "Longitude", y = "Latitude", fill = "") +
  annotation_scale(location = "bl", width_hint = 0.3) +
  theme_bw() +
  theme(
    title = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.position = c(0.85, 0.25),
    axis.title = element_text(size = 12, face = "plain"),
    legend.key.size = unit(0.5, "cm")
  ) +
  guides(colour = "none", size = "none", alpha = "none") +
  labs(title = "PV conflicts")

write_sf(PVpcell, "Spatial_Data/Richness/PVconflictcell.shp")

ggsave("Figures/Gainloss.png", Gainloss, wi = 45, he = 20, un = "cm", dpi = 300)


# FOR WF POWER--------------------------------------------

#Loading Atlas_Ebird 2018-2023 shp
Atlas_Ebird23 <- st_read("Spatial_Data/Richness/At_Eb2023_p.shp")
#Loading WF power per municipality
WF <- st_read("Spatial_Data/Energy/WF_p.shp")
#Transforming to have the same coordinates system
WF <- st_transform(WF, 25830)

# Removing provinces with no MW 
WF1 <- WF[complete.cases(WF$MW),]

#Calculating the area of each intersected section
intersect0WF <- st_intersection(WF1, Atlas_Ebird23)
intersect0WF$areaintsct_km2 <- st_area(intersect0WF)/1000000
#Removing units from rows
intersect0WF <- drop_units(intersect0WF)

#Calculating the proportional MW for each intersection based on the size of the containing cell
intersect0WF$Proport_MW <- intersect0WF$MW * intersect0WF$areaintsct_km2 / intersect0WF$Area_km2

#Transforming the shp into a data frame 
intersectWF  <- select(intersect0WF, UTMCODE, MW, Area_km2, species, areaintsct_km2, Proport_MW)
intersectWF$geometry <- NULL

#Summing the total PV power per cell 
WF_MW_pcell <- intersectWF %>% 
  group_by(UTMCODE) %>% 
  summarize(MW_pcell = sum(Proport_MW)) 

#Merging PropMW_pcell and Atlas2014-2018+Ebird 2019-2023
WFpcell <- left_join(Atlas_Ebird23[, c("UTMCODE", "XCENTROIDE", "YCENTROIDE", "Area_km2", "species")], 
                     WF_MW_pcell, by = c("UTMCODE"))

WFpcell1  <-  WFpcell[complete.cases(WFpcell$MW_pcell), ]

# Calcula los cuartiles (25% y 75%)
WFpcellQ <- quantile(WFpcell1$MW_pcell, probs = c(0.25, 0.75))

# Asigna categorías según los cuartiles
WFpcell1$category <- cut(WFpcell1$MW_pcell,
                         breaks = c(-Inf, WFpcellQ[1], WFpcellQ[2], Inf),
                         labels = c("low", "medium", "high"))

conflictWF <- ggplot() +
  geom_sf(data = WFpcell1, aes(fill = category), color = NA) + 
  geom_sf(data = datosMAP1, fill = "transparent", color = "red") +
  geom_sf(data = comunidades_autonomas, fill = "transparent", color = "black") +# Utiliza la columna real de tus datos
  scale_fill_manual(values = c("low" = "thistle", "medium" = "lightblue", "high" = "blue")) +  # Asigna colores a cada categoría
  labs(x = "Longitude", y = "Latitude", fill = "") +
  annotation_scale(location = "bl", width_hint = 0.3) +
  theme_bw() +
  theme(
    title = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.position = c(0.85, 0.25),
    axis.title = element_text(size = 12, face = "plain"),
    legend.key.size = unit(0.5, "cm")
  ) +
  guides(colour = "none", size = "none", alpha = "none") +
  labs(title = "WF conflicts")

# FOR PV POWER IN PENINSULAR SPAIN--------------------------------------------

#Loading Atlas_Ebird 2018-2023 shp
Atlas_Ebird23 <- st_read("Spatial_Data/Richness/At_Eb2023_p.shp")
#Loading PV power per municipality
PV <- st_read("Spatial_Data/Energy/PV_planned_p.shp")
#Transforming to have the same coordinates system
PV <- st_transform(PV, 25830)

# Removing provinces with no MW 
PV1 <- PV[complete.cases(PV$MW),]

#Calculating the area of each intersected section
intersect0 <- st_intersection(PV1, Atlas_Ebird23)
intersect0$areaintsct_km2 <- st_area(intersect0)/1000000
#Removing units from rows
intersect0 <- drop_units(intersect0)

#Calculating the proportional MW for each intersection based on the size of the containing cell
intersect0$Proport_MW <- intersect0$MW * intersect0$areaintsct_km2 / intersect0$Area_km2

#Transforming the shp into a data frame 
intersect  <- select(intersect0, UTMCODE, MW, Area_km2, species, areaintsct_km2, Proport_MW)
intersect$geometry <- NULL

#Summing the total PV power per cell 
PV_MW_pcell <- intersect %>% 
  group_by(UTMCODE) %>% 
  summarize(MW_pcell = sum(Proport_MW)) 

#Merging PropMW_pcell and Atlas2014-2018+Ebird 2019-2023
PVpcell_p <- left_join(Atlas_Ebird23[, c("UTMCODE", "XCENTROIDE", "YCENTROIDE", "Area_km2", "species")], 
                     PV_MW_pcell, by = c("UTMCODE"))

PVpcell_p  <-  PVpcell_p[complete.cases(PVpcell_p$MW_pcell), ]

# Calcula los cuartiles (25% y 75%)
PVpcellQ <- quantile(PVpcell_p$MW_pcell, probs = c(0.25, 0.75))

# Asigna categorías según los cuartiles
PVpcell_p$category <- cut(PVpcell_p$MW_pcell,
                         breaks = c(-Inf, PVpcellQ[1], PVpcellQ[2], Inf),
                         labels = c("low", "medium", "high"))

conflictPVp <- ggplot() +
  geom_sf(data = PVpcell_p, aes(fill = category), color = NA) + 
  geom_sf(data = datosMAP1, fill = "transparent", color = "red") +
  geom_sf(data = comunidades_autonomas, fill = "transparent", color = "black") +# Utiliza la columna real de tus datos
  scale_fill_manual(values = c("low" = "thistle", "medium" = "lightblue", "high" = "blue")) +  # Asigna colores a cada categoría
  labs(x = "Longitude", y = "Latitude", fill = "") +
  annotation_scale(location = "bl", width_hint = 0.3) +
  theme_bw() +
  theme(
    title = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.position = c(0.85, 0.25),
    axis.title = element_text(size = 12, face = "plain"),
    legend.key.size = unit(0.5, "cm")
  ) +
  guides(colour = "none", size = "none", alpha = "none") +
  labs(title = "Planned PV conflicts")


#Loading Atlas_Ebird 2018-2023 shp
Atlas_Ebird23 <- st_read("Spatial_Data/Richness/At_Eb2023_p.shp")
#Loading WF power per municipality
WF <- st_read("Spatial_Data/Energy/WF_planned_p.shp")
#Transforming to have the same coordinates system
WF <- st_transform(WF, 25830)

# Removing provinces with no MW 
WF1 <- WF[complete.cases(WF$MW),]

#Calculating the area of each intersected section
intersect0WF <- st_intersection(WF1, Atlas_Ebird23)
intersect0WF$areaintsct_km2 <- st_area(intersect0WF)/1000000
#Removing units from rows
intersect0WF <- drop_units(intersect0WF)

#Calculating the proportional MW for each intersection based on the size of the containing cell
intersect0WF$Proport_MW <- intersect0WF$MW * intersect0WF$areaintsct_km2 / intersect0WF$Area_km2

#Transforming the shp into a data frame 
intersectWF  <- select(intersect0WF, UTMCODE, MW, Area_km2, species, areaintsct_km2, Proport_MW)
intersectWF$geometry <- NULL

#Summing the total PV power per cell 
WF_MW_pcell <- intersectWF %>% 
  group_by(UTMCODE) %>% 
  summarize(MW_pcell = sum(Proport_MW)) 

#Merging PropMW_pcell and Atlas2014-2018+Ebird 2019-2023
WFpcell_p <- left_join(Atlas_Ebird23[, c("UTMCODE", "XCENTROIDE", "YCENTROIDE", "Area_km2", "species")], 
                     WF_MW_pcell, by = c("UTMCODE"))

WFpcell_p1  <-  WFpcell_p[complete.cases(WFpcell_p$MW_pcell), ]

# Calcula los cuartiles (25% y 75%)
WFpcellpQ <- quantile(WFpcell_p1$MW_pcell, probs = c(0.25, 0.75))

# Asigna categorías según los cuartiles
WFpcell_p1$category <- cut(WFpcell_p1$MW_pcell,
                          breaks = c(-Inf, WFpcellpQ[1], WFpcellpQ[2], Inf),
                          labels = c("low", "medium", "high"))

conflictWFp <- ggplot() +
  geom_sf(data = WFpcell_p1, aes(fill = category), color = NA) + 
  geom_sf(data = datosMAP1, fill = "transparent", color = "red") +
  geom_sf(data = comunidades_autonomas, fill = "transparent", color = "black") +# Utiliza la columna real de tus datos
  scale_fill_manual(values = c("low" = "thistle", "medium" = "lightblue", "high" = "blue")) +  # Asigna colores a cada categoría
  labs(x = "Longitude", y = "Latitude", fill = "") +
  annotation_scale(location = "bl", width_hint = 0.3) +
  theme_bw() +
  theme(
    title = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.position = c(0.85, 0.25),
    axis.title = element_text(size = 12, face = "plain"),
    legend.key.size = unit(0.5, "cm")
  ) +
  guides(colour = "none", size = "none", alpha = "none") +
  labs(title = "Planned WF conflicts")


library(sf)
library(ggplot2)

# Intersección entre datosMAP1 y WFpcell_p1
interseccion <- st_intersection(datosMAP1, WFpcell_p1)

# Graficar la intersección
ggplot() +
  geom_sf(data = WFpcell_p1, aes(fill = category), color = NA) + 
  geom_sf(data = interseccion, aes(fill = category), color = "red") +
  geom_sf(data = comunidades_autonomas, fill = "transparent", color = "black") +
  scale_fill_manual(values = c("low" = "thistle", "medium" = "lightblue", "high" = "blue")) +
  labs(x = "Longitude", y = "Latitude", fill = "") +
  annotation_scale(location = "bl", width_hint = 0.3) +
  theme_bw() +
  theme(
    title = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.position = c(0.85, 0.25),
    axis.title = element_text(size = 12, face = "plain"),
    legend.key.size = unit(0.5, "cm")
  ) +
  guides(colour = "none", size = "none", alpha = "none") +
  labs(title = "Planned WF conflicts")

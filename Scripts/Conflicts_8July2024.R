library(tmap)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(sf)
library(mapSpain)
library(units)
library(writexl)
setwd("~/GitHub/BiodiversityFacets")

####################################################################
############### Overlapping the 3 indexes############################

# clean environment
rm(list = ls())

TD <- st_read("Spatial_Data/Richness/TDindex.shp")
FD <- st_read("Spatial_Data/Richness/FDindex_SES.shp")
PD <- st_read("Spatial_Data/Richness/PhylDiversity_sizeeffect.shp")

# Adding Autonomous communities from Spain
comm <-esp_get_ccaa()
comm <- comm[!comm$iso2.ccaa.name.es %in% c("Canarias"),]
comm <- st_transform(comm, 25830) 

prov <-esp_get_prov()
prov <- prov[!prov$iso2.prov.name.es %in% c("Las Palmas", "Santa Cruz de Tenerife"),]
prov <- st_transform(prov, 25830) 

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

datosMAP <- datosMAP[, c("UTMCODE","species", "FRic", "SESFRic", "pd_bs_z", "pd_obs")]
datosMAP <- st_transform(datosMAP, 25830) 

datosMAP$area_cell <- st_area(datosMAP)

#Removing units from rows
datosMAP <- drop_units(datosMAP)

# load Spanish UTM grid
malla <- st_read("Spatial_Data/Malla_municipios/Malla10x10_Ter_p.shp")
malla <- st_transform(malla, 25830) 

write_sf(datosMAP, "Spatial_Data/3facets.shp")

###########Using photovoltaic polygons areas to identify conflicts #################################
###################################################################################################

#Loading solar plants
solar <- st_read("Spatial_Data/Energy/allsolar_cut_diss.shp")
solar <- st_transform(solar, 25830) 
names(solar)[1] <- "area_PVplant"

#Loading TD, FD, PD data per cell
datosMAP <- st_read("Spatial_Data/3facets.shp")
datosMAP$area_cell <- st_area(datosMAP)

#Removing units from rows
datosMAP <- drop_units(datosMAP)

# Intersecting solar plants polygons and cells with 3 biodiv facets 
PVint <- st_intersection(solar, datosMAP) 
PVint$area_int <- st_area(PVint)
PVint <- drop_units(PVint)
PVint$geometry <- NULL

#Summing the total PV area per cell 
PVint1 <- PVint %>% 
  group_by(UTMCODE) %>% 
  summarize(PVarea_pcell = sum(area_int)) 

PVint1 <- merge(PVint, PVint1, by= "UTMCODE")
PVint1 <- PVint1 %>%
  distinct(UTMCODE, .keep_all = TRUE)

#Calculating the percentage of UTM area occupied by PV plants
PVint1$PercPV_cell <- PVint1$PVarea_pcell * 100 / PVint1$area_cell

PVint1_subset <- PVint1[, c("UTMCODE", "area_PVplant", "area_int", "PVarea_pcell", "PercPV_cell")]

# Merging the dataframe with intersection data to the one with TD, FD, and PD
PVint11 <- merge(datosMAP, PVint1_subset, by = "UTMCODE", all.x = TRUE)
names(PVint11)[5] <- "SESPD"
names(PVint11)[6] <- "PD"

###########Using photovoltaic power per municipality to identify conflicts #################################
###################################################################################################

PVint2 <- PVint11[, c("UTMCODE","species", "FRic", "SESFRic", "SESPD", "PD", "PercPV_cell")]

# Calculating percentil 70 for each index
cutoff_sp_richn <- quantile(PVint2$species, probs = 0.70, na.rm = TRUE)
cutoff_SESfric <- quantile(PVint2$SESFRic, probs = 0.70, na.rm = TRUE)
cutoff_SESPD <- quantile(PVint2$SESPD, probs = 0.70, na.rm = TRUE)
cutoff_fric <- quantile(PVint2$FRic, probs = 0.70, na.rm = TRUE)
cutoff_PD <- quantile(PVint2$PD, probs = 0.70, na.rm = TRUE)

# New variable for each index. It is TRUE if the value is inside the top 30%
PVint2$top_sp_richn <- PVint2$species > cutoff_sp_richn
PVint2$top_SESfric <- PVint2$SESFRic > cutoff_SESfric
PVint2$top_SESPD <- PVint2$SESPD > cutoff_SESPD
PVint2$top_fric <- PVint2$FRic > cutoff_fric
PVint2$top_PD <- PVint2$PD > cutoff_PD

# Creating a new variable to represent overlapping using SES values
PVint2$overlapSES <- with(PVint2, ifelse(top_sp_richn & top_SESfric & top_SESPD, "Hotspot",
                                          ifelse(top_sp_richn & top_SESfric, "Top 30% sp_richn and SES_fric",
                                                 ifelse(top_sp_richn & top_SESPD, "Top 30% sp_richn and SES_PD",
                                                        ifelse(top_SESfric & top_SESPD, "Top 30% SES_fric and SES_PD",
                                                               ifelse(top_sp_richn, "Top 30% sp_richn",
                                                                      ifelse(top_SESfric, "Top 30% SES_fric",
                                                                             ifelse(top_SESPD, "Top 30% SES_PD", "Outside hotspots"))))))))
# Sorting overlapSES levels
PVint2$overlapSES <- factor(PVint2$overlap, levels = c("Hotspot", "Top 30% sp_richn and SES_fric", 
                                                    "Top 30% sp_richn and SES_PD", "Top 30% SES_fric and SES_PD", 
                                                    "Top 30% sp_richn", "Top 30% SES_fric", "Top 30% SES_PD", 
                                                    "Outside hotspots"))

# Creating a new variable to represent overlapping using raw values
PVint2$overlap <- with(PVint2, ifelse(top_sp_richn & top_fric & top_PD, "Hotspot",
                                      ifelse(top_sp_richn & top_fric, "Top 30% sp_richn and fric",
                                             ifelse(top_sp_richn & top_PD, "Top 30% sp_richn and PD",
                                                    ifelse(top_fric & top_PD, "Top 30% fric and PD",
                                                           ifelse(top_sp_richn, "Top 30% sp_richn",
                                                                  ifelse(top_fric, "Top 30% fric",
                                                                         ifelse(top_PD, "Top 30% PD", "Outside hotspots"))))))))

# Sorting overlap levels
PVint2$overlap <- factor(PVint2$overlap, levels = c("Hotspot", "Top 30% sp_richn and fric", 
                                                       "Top 30% sp_richn and PD", "Top 30% fric and PD", 
                                                       "Top 30% sp_richn", "Top 30% fric", "Top 30% PD", 
                                                       "Outside hotspots"))

# Colours for each overlapping category 
colores <- c("Hotspot" = "blue4", "Top 30% sp_richn and SES_fric" = "blue", "Top 30% sp_richn and fric" = "blue","Top 30% sp_richn and SES_PD" = "dodgerblue", 
             "Top 30% sp_richn and PD" = "dodgerblue","Top 30% SES_fric and SES_PD" = "lightskyblue", "Top 30% fric and PD" = "lightskyblue", "Top 30% sp_richn" = "cyan3", "Top 30% SES_fric" = "paleturquoise3",
             "Top 30% fric" = "paleturquoise3", "Top 30% SES_PD" = "azure3", "Top 30% PD" = "azure3", "Outside hotspots" = "white")

PVint2$overlapSES[is.na(PVint2$overlapSES)] <- "Outside hotspots"
PVint2$overlap[is.na(PVint2$overlap)] <- "Outside hotspots"

# Final map
mapaSES <- tm_shape(PVint2) +
  tm_fill(col = "overlapSES", palette = colores, title = "Index overlapping", style = "cat")

# Adding Spain borders to final map
mapaSES <- tm_shape(PVint2) +
  tm_fill(col = "overlapSES", palette = colores, title = "Index overlapping", style = "cat") +
  tm_shape(prov) +
  tm_borders() +
  tm_scale_bar(breaks = c(0, 50, 100, 150, 200), position = c("left", "bottom")) +  # Agregar barra de escala con intervalos personalizados
  tm_compass(type = "arrow", position = c("right", "top"))

tmap_save(mapaSES, filename = "Figures/3facetsSEShotsp.png")

# Final map
mapa <- tm_shape(PVint2) +
  tm_fill(col = "overlap", palette = colores, title = "Index overlapping", style = "cat")

# Adding Spain borders to final map
mapa <- tm_shape(PVint2) +
  tm_fill(col = "overlap", palette = colores, title = "Index overlapping", style = "cat") +
  tm_shape(prov) +
  tm_borders() +
  tm_scale_bar(breaks = c(0, 50, 100, 150, 200), position = c("left", "bottom")) +  # Agregar barra de escala con intervalos personalizados
  tm_compass(type = "arrow", position = c("right", "top"))

tmap_save(mapa, filename = "Figures/3facetshotsp.png")

# Variables in terciles
terciles <- quantile(PVint2$PercPV_cell, probs = c(0, 0.33, 0.67, 1), na.rm = TRUE)

# Creating the new column to show high, medium and low occupacy of PV plants inside each cell based on terciles.
PVint2 <- PVint2 %>%
  mutate(PV_occupancy = cut(PercPV_cell, breaks = terciles, labels = c("low", "medium", "high"), include.lowest = TRUE))

# New variable to represent overlapping (high is for the overlap of the 3 indexes, medium for the overlap of any 2 indexes, and low for 1 index)
PVint2$div_cat <- with(PVint2, ifelse(top_sp_richn & top_SESfric & top_SESPD, "high",
                                          ifelse(top_sp_richn & top_SESfric, "medium",
                                                 ifelse(top_sp_richn & top_SESPD, "medium",
                                                        ifelse(top_SESfric & top_SESPD, "medium",
                                                               ifelse(top_sp_richn, "low",
                                                                      ifelse(top_SESfric, "low",
                                                                             ifelse(top_SESPD, "low", "NotA"))))))))

# I am assigning 3 levels of conflict based on the presence or absence of PV plants in each cell. I am not taking into account
# the area of each PV, just the presence or absence. If there is a PV in a cell with high TD, FD, and PD (hotspot), 
# then strong conflict is assigned, high TD-PD, TD-FD, or FD-PD cells with PV plants are assigned as moderate conflict,
# while cells with just one high diversity facet and the presence of PV plant are assigned the term "conflict". 
PVint2 <- PVint2 %>%
  mutate(conflict = case_when(
    div_cat == "high" & !is.na(PV_occupancy) ~ "strong conflict",
    div_cat == "medium" & !is.na(PV_occupancy) ~ "moderate conflict",
    div_cat == "low" & !is.na(PV_occupancy) ~ "conflict",
    div_cat == "high" & is.na(PV_occupancy) ~ "no-go",
    TRUE ~ "Outside hotspots"  # Optional: fill with NA if no condition is met
  ))
# Ordenar los niveles de la variable overlap
PVint2$conflict <- factor(PVint2$conflict, levels = c("strong conflict", "moderate conflict", 
                                                    "conflict", "no-go", "Outside hotspots"))
# Colours for each overlapping category 
colors <- c("strong conflict" = "black", "moderate conflict" = "orange", "conflict" = "wheat", "no-go" = "honeydew4", "Outside hotspots" = "white")

# Final map
mapconflict <- tm_shape(PVint2) +
  tm_fill(col = "conflict", palette = colors, title = "", style = "cat")

# Adding Spain borders to final map
mapconflict <- tm_shape(PVint2) +
  tm_fill(col = "conflict", palette = colors, title = "", style = "cat") +
  tm_shape(prov) +
  tm_borders() +
  tm_scale_bar(breaks = c(0, 50, 100, 150, 200), position = c("left", "bottom")) +  # Agregar barra de escala con intervalos personalizados
  tm_compass(type = "arrow", position = c("right", "top"))

tmap_save(mapconflict, filename = "Figures/3facetsSESconflict.png")


# I am assigning 9 levels of conflict based on 3, 2 and 1 facet diverse cells, and the occupancy levels of PV plants per cell (high, medium and low). 

PVint2 <- PVint2 %>%
  mutate(extendedconflicts = case_when(
    div_cat == "high" & PV_occupancy == "high" ~ "Strong conflict 3", # 3 facet biodiverse areas and high PV plants occupancy
    div_cat == "high" & PV_occupancy == "medium" ~ "Strong conflict 2", # 3 facet biodiverse areas and medium PV plants occupancy
    div_cat == "high" & PV_occupancy == "low" ~ "Strong conflict 1", # 3 facet biodiverse areas and low PV plants occupancy
    div_cat == "medium" & PV_occupancy == "high" ~ "Moderate conflict 3", # 2 facet biodiverse areas and high PV plants occupancy
    div_cat == "medium" & PV_occupancy == "medium" ~ "Moderate conflict 2", # 2 facet biodiverse areas and medium PV plants occupancy
    div_cat == "medium" & PV_occupancy == "low" ~ "Moderate conflict 1", # 3 facet biodiverse areas and high PV plants occupancy
    div_cat == "low" & PV_occupancy == "high" ~ "Conflict 3", # 3 facet biodiverse areas and high PV plants occupancy
    div_cat == "low" & PV_occupancy == "medium" ~ "Conflict 2", # 3 facet biodiverse areas and high PV plants occupancy
    div_cat == "low" & PV_occupancy == "low" ~ "Conflict 1", # 3 facet biodiverse areas and high PV plants occupancy
    TRUE ~ "Outside hotspots"  # Optional: fill with NA if no condition is met
  ))

# Sorting conflict levels
PVint2$extendedconflicts <- factor(PVint2$extendedconflicts, levels = c("Strong conflict 3", "Strong conflict 2", 
                                                        "Strong conflict 1", "Moderate conflict 3", "Moderate conflict 2", "Moderate conflict 1",
                                                        "Conflict 3", "Conflict 2", "Conflict 1", "Outside hotspots"))
# Colours for each overlapping category 
colors <- c("Strong conflict 3" = "black", "Strong conflict 2" = "gray21", "Strong conflict 1"  = "gray41", "Moderate conflict 3" = "sienna1", 
            "Moderate conflict 2" = "lightsalmon", "Moderate conflict 1" = "tan", "Conflict 3" = "khaki3", "Conflict 2" = "peachpuff3", 
            "Conflict 1" = "wheat", "Outside hotspots" = "white")

# Crear una variable numérica para el conflicto
PVint2 <- PVint2 %>%
  mutate(conflict_level = case_when(
    extendedconflicts == "Strong conflict 3" ~ 9,
    extendedconflicts == "Strong conflict 2" ~ 8,
    extendedconflicts == "Strong conflict 1" ~ 7,
    extendedconflicts == "Moderate conflict 3" ~ 6,
    extendedconflicts == "Moderate conflict 2" ~ 5,
    extendedconflicts == "Moderate conflict 1" ~ 4,
    extendedconflicts == "Conflict 3" ~ 3,
    extendedconflicts == "Conflict 2" ~ 2,
    extendedconflicts == "Conflict 1" ~ 1,
    TRUE ~ 0  # Para "Outside hotspots"
  ))

# Crear el mapa usando una escala continua de rojos
mapextendedconflict <- tm_shape(PVint2) +
  tm_fill(col = "conflict_level", palette = "Reds", title = "Conflict Level", style = "cont",
  breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)) +
  tm_shape(prov) +
  tm_borders() +
  tm_scale_bar(breaks = c(0, 50, 100, 150, 200), position = c("left", "bottom")) +
  tm_compass(type = "arrow", position = c("right", "top"))


# Guardar el mapa
tmap_save(mapextendedconflict, filename = "Figures/extendedconflicts.png")













#Now creating risk maps for no-go areas. This is where no PV plants exist, but risk levels are assigned according to the diversity of the 3 facets
#High diversity is where the 3 facets overlap, medium where 2 indexes overlap, and one represents for the top values of just one index.

PVint2 <- PVint2 %>%
  mutate(nogo_risk = case_when(
    div_cat == "high" & is.na(PV_occupancy) ~ "no go: high risk", # 3 facet biodiverse areas where no PV plants exist
    div_cat == "medium" & is.na(PV_occupancy) ~ "no go: moderate risk", # 2 facet biodiverse areas where no PV plants exist
    div_cat == "low" & is.na(PV_occupancy) ~ "no go: risk", # 1 facet biodiverse areas where no PV plants exist
    TRUE ~ "Outside hotspots"  # Optional: fill with NA if no condition is met
  ))

# Ordenar los niveles de la variable overlap
PVint2$nogo_risk <- factor(PVint2$nogo_risk, levels = c("no go: high risk", "no go: moderate risk", 
                                                      "no go: risk", "Outside hotspots"))
# Colours for each overlapping category 
colors <- c("no go: high risk" = "red", "no go: moderate risk" = "orange", "no go: risk" = "wheat", "Outside hotspots" = "white")

# Final map
mapnogo_risk <- tm_shape(PVint2) +
  tm_fill(col = "nogo_risk", palette = colors, title = "", style = "cat")

# Adding Spain borders to final map
mapnogo_risk <- tm_shape(PVint2) +
  tm_fill(col = "nogo_risk", palette = colors, title = "", style = "cat") +
  tm_shape(prov) +
  tm_borders() +
  tm_scale_bar(breaks = c(0, 50, 100, 150, 200), position = c("left", "bottom")) +
  tm_compass(type = "arrow", position = c("right", "top"))

tmap_save(mapnogo_risk, filename = "Figures/nogo_risks.png")






forcorrmatrix <- PVint2[, c("species", "FRic", "SESFRic", "SESPD", "PD")]
forcorrmatrix$geometry <- NULL
corrmatrix <- cor(forcorrmatrix, use = "complete.obs")

#Very high correlation between TD - FD = 80%, and TD - PD = 93%.

#I want separate maps for the top 30% of each diversity facet

#TOP TD

PVint2$top_sp_richn[is.na(PVint2$top_sp_richn)] <- FALSE

# Colours for each overlapping category 
colors1 <- c("TRUE" = "cyan3", "FALSE" = "white", "NA" = "white")

# Final map
mapTopTD <- tm_shape(PVint2) +
  tm_fill(col = "top_sp_richn", palette = colors1, title = "Top 30% TD", style = "cat")

# Adding Spain borders to final map
mapTopTD <- tm_shape(PVint2) +
  tm_fill(col = "top_sp_richn", palette = colors1, title = "Top 30% TD", style = "cat") +
  tm_shape(prov) +
  tm_borders() +
  tm_scale_bar(breaks = c(0, 50, 100, 150, 200), position = c("left", "bottom")) +  # Agregar barra de escala con intervalos personalizados
  tm_layout(legend.show = FALSE)+
  tm_credits(text = "Top 30% TD", position = c("center", "bottom"), size = 1.2) 

tmap_save(mapTopTD, filename = "Figures/mapTopTD.png")

#TOP SESFD

PVint2$top_SESfric[is.na(PVint2$top_SESfric)] <- FALSE

# Colours for each overlapping category 
colors2 <- c("TRUE" = "paleturquoise3", "FALSE" = "white")

# Final map
mapTopFD <- tm_shape(PVint2) +
  tm_fill(col = "top_SESfric", palette = colors2, title = "", style = "cat")

# Adding Spain borders to final map
mapTopFD <- tm_shape(PVint2) +
  tm_fill(col = "top_SESfric", palette = colors2, title = "", style = "cat") +
  tm_shape(prov) +
  tm_borders() +
  tm_scale_bar(breaks = c(0, 50, 100, 150, 200), position = c("left", "bottom")) +  # Agregar barra de escala con intervalos personalizados
  tm_layout(legend.show = FALSE)+
  tm_credits(text = "Top 30% SESFRic", position = c("center", "bottom"), size = 1.2) 

tmap_save(mapTopFD, filename = "Figures/mapTopSESFD.png")

#TOP FD

PVint2$top_fric[is.na(PVint2$top_fric)] <- FALSE

# Colours for each overlapping category 
colors2 <- c("TRUE" = "paleturquoise3", "FALSE" = "white")

# Final map
mapTopFD1 <- tm_shape(PVint2) +
  tm_fill(col = "top_fric", palette = colors2, title = "", style = "cat")

# Adding Spain borders to final map
mapTopFD <- tm_shape(PVint2) +
  tm_fill(col = "top_fric", palette = colors2, title = "", style = "cat") +
  tm_shape(prov) +
  tm_borders() +
  tm_scale_bar(breaks = c(0, 50, 100, 150, 200), position = c("left", "bottom")) +  # Agregar barra de escala con intervalos personalizados
  tm_layout(legend.show = FALSE)+
  tm_credits(text = "Top 30% SESFRic", position = c("center", "bottom"), size = 1.2) 

tmap_save(mapTopFD, filename = "Figures/mapTopFD.png")

#TOP SESPD

PVint2$top_SESPD[is.na(PVint2$top_SESPD)] <- FALSE

# Colours for each overlapping category 
colors3 <- c("TRUE" = "azure3", "FALSE" = "white")

# Final map
mapTopPD <- tm_shape(PVint2) +
  tm_fill(col = "top_SESPD", palette = colors3, title = "", style = "cat")

# Adding Spain borders to final map
mapTopPD <- tm_shape(PVint2) +
  tm_fill(col = "top_SESPD", palette = colors3, title = "", style = "cat") +
  tm_shape(prov) +
  tm_borders() +
  tm_scale_bar(breaks = c(0, 50, 100, 150, 200), position = c("left", "bottom")) +  
  tm_layout(legend.show = FALSE)+
  tm_credits(text = "Top 30% SESPD", position = c("center", "bottom"), size = 1.2) 

tmap_save(mapTopPD, filename = "Figures/mapTopSESPD.png")

#TOP PD

PVint2$top_PD[is.na(PVint2$top_PD)] <- FALSE

# Colours for each overlapping category 
colors3 <- c("TRUE" = "azure3", "FALSE" = "white")

# Final map
mapTopPD <- tm_shape(PVint2) +
  tm_fill(col = "top_PD", palette = colors3, title = "", style = "cat")

# Adding Spain borders to final map
mapTopPD <- tm_shape(PVint2) +
  tm_fill(col = "top_PD", palette = colors3, title = "", style = "cat") +
  tm_shape(prov) +
  tm_borders() +
  tm_scale_bar(breaks = c(0, 50, 100, 150, 200), position = c("left", "bottom")) +  
  tm_layout(legend.show = FALSE)+
  tm_credits(text = "Top 30% PD", position = c("center", "bottom"), size = 1.2) 

tmap_save(mapTopPD, filename = "Figures/mapTopPD.png")

# I want to know the number of hotspot cells inside each province, and their corresponding percentage in relation to 
# the number of cells per province and in relation to the total hotspot cells in Spain. Also the corresponding areas. 

#Selecting only hotspot cells

Hotspotcells <- subset(PVint2, overlapSES == "Hotspot") #There are 242 hotspot cells in Spain

# Spatial intersect between provinces and hotspot cells
hotsp_prov_int <- st_intersection(Hotspotcells, prov)

#Calculating the area of each intersection section
hotsp_prov_int$areainthots_prov <- st_area(hotsp_prov_int)

# Spatial intersect between provinces and malla to know the number of cells per province
prov_int <- st_intersection(malla, prov)

cells_perprovince <- prov_int %>% 
  group_by(prov.shortname.en) %>% 
  summarize(UTMCODE = n())

names(cells_perprovince)[2]="N cells"
cells_perprovince$geometry <- NULL

hotspots_perprovince <- hotsp_prov_int %>% 
  group_by(prov.shortname.en) %>% 
  summarize(UTMCODE = n())

names(hotspots_perprovince)[2]="N hotspot cells"

#Percentage of hotspot cells in relation to the number of cells per province
Province_stats <- merge(cells_perprovince, hotspots_perprovince, by="prov.shortname.en")

Province_stats$Perc_hotp_prov <- Province_stats$`N hotspot cells`*100/Province_stats$`N cells`

#Percentage of hotspot cells in relation to the total number of hotspot cells in Spain
Province_stats$Perc_hotp_Spain <- Province_stats$`N hotspot cells`*100/242

#As one hotspot cell may fall in different provinces, calculating areas could be more precise. 

# Calculating the total area from hotspots cells
areaHots <- st_area(Hotspotcells)
sum(areaHots) #The 242 hotspot cells equal to 24018154101 [m^2] or 24018.15 km2

# Calculating the area from each province 
prov$areaprovince <- st_area(prov)

#Summing the total hotspot area per province 
areas_prov_hotsp <- hotsp_prov_int %>% 
  group_by(prov.shortname.en) %>% 
  summarize(Hotsarea_pprov = sum(areainthots_prov)) 

areas_prov_hotsp$geometry <- NULL

areas_prov_hotsp1 <- merge(areas_prov_hotsp, prov[,c("prov.shortname.en", "areaprovince")], by="prov.shortname.en" ) 
areas_prov_hotsp1 <- drop_units(areas_prov_hotsp1)
Province_stats <- merge(Province_stats, areas_prov_hotsp1, by="prov.shortname.en" ) 
Province_stats <- drop_units(Province_stats)
Province_stats$geometry.x <- NULL

Province_stats$per_hotareaprov <- Province_stats$Hotsarea_pprov*100/Province_stats$areaprovince
Province_stats$per_hotareaSpain <- Province_stats$Hotsarea_pprov*100/24018154101 

# I want to know the number of conflict cells inside each province, and their corresponding percentage in relation to 
# the number of cells per province and in relation to the total conflict cells in Spain. Also the corresponding areas. 

#Selecting only strong conflict cells

stconflictcells <- subset(PVint2, conflict == "strong conflict") #There are 143 strong conflict cells in Spain

# Spatial intersect between provinces and hotspot cells
stconflict_prov_int <- st_intersection(stconflictcells, prov)

#Calculating the area of each intersection section
stconflict_prov_int$areaintconf_prov <- st_area(stconflict_prov_int)

stconflict_perprovince <- stconflict_prov_int %>% 
  group_by(prov.shortname.en) %>% 
  summarize(UTMCODE = n())

names(stconflict_perprovince)[2]="N strong conflict cells"

#Percentage of conflict cells in relation to the number of cells per province
Province_stats1 <- merge(Province_stats, stconflict_perprovince, by="prov.shortname.en", all.x=TRUE)
Province_stats1$Perc_stconf_prov <- Province_stats1$`N strong conflict cells`*100/Province_stats$`N cells`

#Percentage of conflict cells in relation to the total number of conflict cells in Spain
Province_stats1$Perc_stconf_Spain <- Province_stats1$`N strong conflict cells`*100/143

Province_stats1$geometry.x <- NULL
Province_stats1$geometry.y <- NULL

#As one conflict cell may fall in different provinces, calculating areas could be more precise. 

# Calculating the total area from hotspots cells
areaconflict <- st_area(stconflictcells)
sum(areaconflict) #The 143 strong conflict cells equal to 14285812972 [m^2] or 14285.812972 km2

#Summing the total hotspot area per province 
areas_prov_confl <- stconflict_prov_int %>% 
  group_by(prov.shortname.en) %>% 
  summarize(Confarea_pprov = sum(areaintconf_prov)) 

areas_prov_confl$geometry <- NULL

areas_prov_confl1 <- merge(areas_prov_confl, prov[,c("prov.shortname.en", "areaprovince")], by="prov.shortname.en" ) 
areas_prov_confl1 <- drop_units(areas_prov_confl1)
Province_stats1 <- merge(Province_stats1, areas_prov_confl1[,c("prov.shortname.en", "Confarea_pprov")], by="prov.shortname.en", all.x=TRUE ) 
Province_stats1 <- drop_units(Province_stats1)

Province_stats1$per_confareaprov <- Province_stats1$Confarea_pprov*100/Province_stats1$areaprovince
Province_stats1$per_confareaSpain <- Province_stats1$Confarea_pprov*100/14285812972 

# No-go areas

#Selecting only no-go cells
nogocells <- subset(PVint2, conflict == "no-go") #There are 99 no-go cells in Spain

# Spatial intersect between provinces and no-go cells
nogo_prov_int <- st_intersection(nogocells, prov)
#Calculating the area of each intersection section
nogo_prov_int$areaintnogo_prov <- st_area(nogo_prov_int)

nogo_perprovince <- nogo_prov_int %>% 
  group_by(prov.shortname.en) %>% 
  summarize(UTMCODE = n())

names(nogo_perprovince)[2]="no-go cells"

#Percentage of no-go cells in relation to the number of cells per province
Province_stats2 <- merge(Province_stats1, nogo_perprovince, by="prov.shortname.en", all.x=TRUE)
Province_stats2$Perc_nogo_prov <- Province_stats2$`no-go cells`*100/Province_stats2$`N cells`

#Percentage of no-go cells in relation to the number of cells per province
Province_stats2 <- merge(Province_stats1, nogo_perprovince, by="prov.shortname.en", all.x=TRUE)
Province_stats2$Perc_nogo_prov <- Province_stats2$`no-go cells`*100/Province_stats2$`N cells`
Province_stats2$Perc_nogo_Spain <- Province_stats2$`no-go cells`*100/99

Province_stats2$geometry <- NULL

#As one no-go cell may fall in different provinces, calculating areas could be more precise. 

# Calculating the total area from hotspots cells
areanogo <- st_area(nogocells)
sum(areanogo) #The 99 no-go cells equal to 9732341129 [m^2] or 9732.341129 km2

#Summing the total nogo area per province 
areas_prov_nogo <- nogo_prov_int %>% 
  group_by(prov.shortname.en) %>% 
  summarize(nogoarea_pprov = sum(areaintnogo_prov)) 

areas_prov_nogo$geometry <- NULL

areas_prov_nogo1 <- merge(areas_prov_nogo, prov[,c("prov.shortname.en", "areaprovince")], by="prov.shortname.en") 
areas_prov_nogo1 <- drop_units(areas_prov_nogo1)
Province_stats3 <- merge(Province_stats2, areas_prov_nogo1[,c("prov.shortname.en", "nogoarea_pprov")], by="prov.shortname.en", all.x=TRUE ) 
Province_stats3 <- drop_units(Province_stats3)

Province_stats3$per_nogoareaprov <- Province_stats3$nogoarea_pprov*100/Province_stats3$areaprovince
Province_stats3$per_nogoareaSpain <- Province_stats3$nogoarea_pprov*100/9732341129 

write_xlsx(Province_stats3, 'Data/Province_stats_22july.xlsx') 

library(tmap)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(sf)
library(mapSpain)
library(units)
library(writexl)
library(purrr)
library(raster)
library(reshape2)
library(ggspatial)
library(tidyr)

setwd("~/GitHub/BiodiversityFacets")

####################################################################
############### Overlapping the 3 indexes############################

# clean environment
rm(list = ls())

# Adding Autonomous communities from Spain
AC <-esp_get_ccaa()
AC <- AC[!AC$iso2.ccaa.name.es %in% c("Canarias"),]
AC <- st_transform(AC, 25830) 

# load Spanish UTM grid
malla <- st_read("Spatial_Data/Malla_municipios/Malla10x10_Ter_p.shp")
malla <- st_transform(malla, 25830) 

###########Using photovoltaic polygons areas to identify conflicts #################################
###################################################################################################

#Loading solar plants
solar <- st_read("Spatial_Data/Energy/allsolar_cut_diss.shp")
solar <- st_transform(solar, 25830) 
names(solar)[1] <- "area_PVplant"

#Loading TD, FD, PD data per cell
datosMAP <- st_read("Spatial_Data/3facets.shp")

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
                                          ifelse(top_sp_richn & top_SESfric, "Top 30% TD and FD",
                                                 ifelse(top_sp_richn & top_SESPD, "Top 30% TD and PD",
                                                        ifelse(top_SESfric & top_SESPD, "Top 30% FD and PD",
                                                               ifelse(top_sp_richn, "Top 30% TD",
                                                                      ifelse(top_SESfric, "Top 30% FD",
                                                                             ifelse(top_SESPD, "Top 30% PD", "No overlap"))))))))
# Sorting overlapSES levels
PVint2$overlapSES <- factor(PVint2$overlap, levels = c("Hotspot", "Top 30% TD and FD", 
                                                    "Top 30% TD and PD", "Top 30% FD and PD", 
                                                    "Top 30% TD", "Top 30% FD", "Top 30% PD", 
                                                    "No overlap"))

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
colores <- c("Hotspot" = "blue4", "Top 30% TD and FD" = "blue", "Top 30% sp_richn and fric" = "blue","Top 30% TD and PD" = "dodgerblue", 
             "Top 30% sp_richn and PD" = "dodgerblue","Top 30% FD and PD" = "lightskyblue", "Top 30% fric and PD" = "lightskyblue", "Top 30% TD" = "cyan3", "Top 30% FD" = "paleturquoise3",
             "Top 30% fric" = "paleturquoise3", "Top 30% PD" = "azure3", "Top 30% PD" = "azure3", "No overlap" = "white")

PVint2$overlapSES[is.na(PVint2$overlapSES)] <- "No overlap"
PVint2$overlap[is.na(PVint2$overlap)] <- "Outside hotspots"

# Final map with TD, SESFRIC and SESPD
mapaSES <- tm_shape(PVint2) +
  tm_fill(col = "overlapSES", palette = colores, title = "Index overlapping", style = "cat")

# Adding Spain borders to final map
mapaSES <- tm_shape(PVint2) +
  tm_fill(col = "overlapSES", palette = colores, title = "Index overlapping", style = "cat") +
  tm_shape(AC) +
  tm_borders() +
  tm_scale_bar(breaks = c(0, 50, 100, 150, 200), position = c("left", "bottom")) +  # Agregar barra de escala con intervalos personalizados
  tm_compass(type = "arrow", position = c("right", "top"))

tmap_save(mapaSES, filename = "Figures/3facetsSES_23decemb.png")
tmap_save(mapaSES, filename = "Figures/3facetsSES_23decemb.pdf")

# Final map with TD, FRIC and PD
mapa <- tm_shape(PVint2) +
  tm_fill(col = "overlap", palette = colores, title = "Index overlapping", style = "cat")

# Adding Spain borders to final map
mapa <- tm_shape(PVint2) +
  tm_fill(col = "overlap", palette = colores, title = "Index overlapping", style = "cat") +
  tm_shape(AC) +
  tm_borders() +
  tm_scale_bar(breaks = c(0, 50, 100, 150, 200), position = c("left", "bottom")) +  # Agregar barra de escala con intervalos personalizados
  tm_compass(type = "arrow", position = c("right", "top"))

tmap_save(mapa, filename = "Figures/3facetshotsp.png")

# PV occupancy percentages in terciles
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

# Alternatively (possibly the best option), I am assigning 9 levels of conflict based on 3, 2 and 1 facet diverse cells, and the occupancy levels of PV plants per cell (high, medium and low). 

PVint2 <- PVint2 %>%
  mutate(extendedconflicts = case_when(
    div_cat == "high" & PV_occupancy == "high" ~ "Conflict 9", # 3 facet biodiverse areas and high PV plants occupancy
    div_cat == "high" & PV_occupancy == "medium" ~ "Conflict 8", # 3 facet biodiverse areas and medium PV plants occupancy
    div_cat == "high" & PV_occupancy == "low" ~ "Conflict 7", # 3 facet biodiverse areas and low PV plants occupancy
    div_cat == "medium" & PV_occupancy == "high" ~ "Conflict 6", # 2 facet biodiverse areas and high PV plants occupancy
    div_cat == "medium" & PV_occupancy == "medium" ~ "Conflict 5", # 2 facet biodiverse areas and medium PV plants occupancy
    div_cat == "medium" & PV_occupancy == "low" ~ "Conflict 4", # 2 facet biodiverse areas and low PV plants occupancy
    div_cat == "low" & PV_occupancy == "high" ~ "Conflict 3", # 1 facet biodiverse areas and high PV plants occupancy
    div_cat == "low" & PV_occupancy == "medium" ~ "Conflict 2", # 1 facet biodiverse areas and medium PV plants occupancy
    div_cat == "low" & PV_occupancy == "low" ~ "Conflict 1", # 1 facet biodiverse areas and low PV plants occupancy
    TRUE ~ "Outside hotspots"  # Optional: fill with NA if no condition is met
  ))

# Sorting conflict levels
PVint2$extendedconflicts <- factor(PVint2$extendedconflicts, levels = c("Conflict 9", "Conflict 8", 
                                                        "Conflict 7", "Conflict 6", "Conflict 5", "Conflict 4",
                                                        "Conflict 3", "Conflict 2", "Conflict 1", "Outside hotspots"))

# Creating a numeric varible for the conflict
PVint2 <- PVint2 %>%
  mutate(conflict_level = case_when(
    extendedconflicts == "Conflict 9" ~ 9,
    extendedconflicts == "Conflict 8" ~ 8,
    extendedconflicts == "Conflict 7" ~ 7,
    extendedconflicts == "Conflict 6" ~ 6,
    extendedconflicts == "Conflict 5" ~ 5,
    extendedconflicts == "Conflict 4" ~ 4,
    extendedconflicts == "Conflict 3" ~ 3,
    extendedconflicts == "Conflict 2" ~ 2,
    extendedconflicts == "Conflict 1" ~ 1,
    TRUE ~ 0  # Para "Outside hotspots"
  ))

# Crear el mapa usando una escala continua de rojos
mapextendedconflict <- tm_shape(PVint2) +
  tm_fill(col = "conflict_level", palette = "Reds", title = "Conflict Level", style = "cont",
  breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)) +
  tm_shape(AC) +
  tm_borders() +
  tm_scale_bar(breaks = c(0, 50, 100, 150, 200), position = c("left", "bottom")) +
  tm_compass(type = "arrow", position = c("right", "top"))

# Guardar el mapa
tmap_save(mapextendedconflict, filename = "Figures/extendedconflicts.png")

#Now creating risk maps for no-go areas. This is where no PV plants exist, but risk levels are assigned according to the diversity of the 3 facets
#High diversity is where the 3 facets overlap, medium where 2 indexes overlap, and one represents for the top values of just one index.

PVint2 <- PVint2 %>%
  mutate(nogo_risk = case_when(
    div_cat == "high" & is.na(PV_occupancy) ~ "no go: very high risk", # 3 facet biodiverse areas where no PV plants exist
    div_cat == "medium" & is.na(PV_occupancy) ~ "no go: high risk", # 2 facet biodiverse areas where no PV plants exist
    div_cat == "low" & is.na(PV_occupancy) ~ "no go: moderate risk", # 1 facet biodiverse areas where no PV plants exist
    TRUE ~ "Outside hotspots"  # Optional: fill with NA if no condition is met
  ))

# Ordenar los niveles de la variable overlap
PVint2$nogo_risk <- factor(PVint2$nogo_risk, levels = c("no go: very high risk", "no go: high risk", 
                                                      "no go: moderate risk", "Outside hotspots"))
# Colours for each overlapping category 
colors <- c("no go: very high risk" = "red", "no go: high risk" = "orange", "no go: moderate risk" = "wheat", "Outside hotspots" = "white")

# Final map
mapnogo_risk <- tm_shape(PVint2) +
  tm_fill(col = "nogo_risk", palette = colors, title = "", style = "cat")

# Adding Spain borders to final map
mapnogo_risk <- tm_shape(PVint2) +
  tm_fill(col = "nogo_risk", palette = colors, title = "", style = "cat") +
  tm_shape(AC) +
  tm_borders() +
  tm_scale_bar(breaks = c(0, 50, 100, 150, 200), position = c("left", "bottom")) +
  tm_compass(type = "arrow", position = c("right", "top"))

tmap_save(mapnogo_risk, filename = "Figures/nogo_risks.png")

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
  tm_shape(AC) +
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
  tm_shape(AC) +
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
  tm_shape(AC) +
  tm_borders() +
  tm_scale_bar(breaks = c(0, 50, 100, 150, 200), position = c("left", "bottom")) +  
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
  tm_shape(AC) +
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
  tm_shape(AC) +
  tm_borders() +
  tm_scale_bar(breaks = c(0, 50, 100, 150, 200), position = c("left", "bottom")) +  
  tm_layout(legend.show = FALSE)+
  tm_credits(text = "Top 30% PD", position = c("center", "bottom"), size = 1.2) 

tmap_save(mapTopPD, filename = "Figures/mapTopPD.png")

# I want to know the number of hotspot cells inside each AC, and their corresponding percentage in relation to 
# the number of cells per AC and in relation to the total hotspot cells in Spain. Also the corresponding areas. 

#Selecting only hotspot cells

Hotspotcells <- subset(PVint2, div_cat == "high") #There are 240 hotspot cells in Spain
Twofacetdiv <- subset(PVint2, div_cat == "medium") #There are 742 Twofacetdiv cells in Spain
Onefacetdiv <- subset(PVint2, div_cat == "low") #There are 1037 Onefacetdiv cells in Spain

# Spatial intersect between ACs and hotspot cells
hotsp_AC_int <- st_intersection(Hotspotcells, AC)
twofacet_AC_int <- st_intersection(Twofacetdiv, AC)
onefacet_AC_int <- st_intersection(Onefacetdiv, AC)

#Calculating the area of each intersection section
hotsp_AC_int$areainthots_AC <- st_area(hotsp_AC_int)
twofacet_AC_int$areainttwof_AC <- st_area(twofacet_AC_int)
onefacet_AC_int$areaintonef_AC <- st_area(onefacet_AC_int)

# Spatial intersect between ACs and malla to know the number of cells per AC
AC_int <- st_intersection(malla, AC)

cells_perAC <- AC_int %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

names(cells_perAC)[2]="N cells"
cells_perAC$geometry <- NULL

hotspots_perAC <- hotsp_AC_int %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

twofacet_perAC <- twofacet_AC_int %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

onefacet_perAC <- onefacet_AC_int %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

names(hotspots_perAC)[2]="N hotspot cells"
names(twofacet_perAC)[2]="N twofacet cells"
names(onefacet_perAC)[2]="N onefacet cells"

#Percentage of hotspot cells in relation to the number of cells per AC
AC_stats1 <- merge(cells_perAC, hotspots_perAC, by="ccaa.shortname.en")

AC_stats1$Perc_hotp_AC <- AC_stats1$`N hotspot cells`*100/AC_stats1$`N cells`

#Percentage of hotspot cells in relation to the total number of hotspot cells in Spain
AC_stats1$Perc_hotp_Spain <- AC_stats1$`N hotspot cells`*100/240

#Percentage of twofacet cells in relation to the number of cells per AC
AC_stats2 <- merge(cells_perAC, twofacet_perAC, by="ccaa.shortname.en")

AC_stats2$Perc_twof_AC <- AC_stats2$`N twofacet cells`*100/AC_stats2$`N cells`

#Percentage of twofacet cells in relation to the total number of twofacet cells in Spain
AC_stats2$Perc_twof_Spain <- AC_stats2$`N twofacet cells`*100/742

#Percentage of onefacet cells in relation to the number of cells per AC
AC_stats3 <- merge(cells_perAC, onefacet_perAC, by="ccaa.shortname.en")

AC_stats3$Perc_onef_AC <- AC_stats3$`N onefacet cells`*100/AC_stats3$`N cells`

#Percentage of twofacet cells in relation to the total number of onefacet cells in Spain
AC_stats3$Perc_onef_Spain <- AC_stats3$`N onefacet cells`*100/1037

AC_stats0 <- AC_stats1 %>%
  full_join(AC_stats2, by = "ccaa.shortname.en")

AC_stats0$`N cells.x`  <- NULL
AC_stats0$`N cells.y`  <- NULL
AC_stats0$geometry.x  <- NULL
AC_stats0$geometry.y  <- NULL

AC_stats <- AC_stats0 %>%
  full_join(AC_stats3, by = "ccaa.shortname.en")

AC_stats <- merge(cells_perAC, AC_stats, by="ccaa.shortname.en")
AC_stats$`N cells.y`<- NULL
names(AC_stats)[2] <- "N cells"

# I want to know the number of conflict cells inside each AC, and their corresponding percentage in relation to 
# the number of cells per AC and in relation to the total conflict cells in Spain. Also the corresponding areas. 

#Selecting only each category conflict cells

Conflictcells9 <- subset(PVint2, extendedconflicts == "Conflict 9") #There are 55 conflict 9 cells in Spain
Conflictcells8 <- subset(PVint2, extendedconflicts == "Conflict 8") #There are 46 conflict 8 cells in Spain
Conflictcells7 <- subset(PVint2, extendedconflicts == "Conflict 7") #There are 41 conflict 7 cells in Spain
Conflictcells6 <- subset(PVint2, extendedconflicts == "Conflict 6") #There are 141 conflict 7 cells in Spain
Conflictcells5 <- subset(PVint2, extendedconflicts == "Conflict 5") #There are 118 conflict 7 cells in Spain
Conflictcells4 <- subset(PVint2, extendedconflicts == "Conflict 4") #There are 109 conflict 7 cells in Spain
Conflictcells3 <- subset(PVint2, extendedconflicts == "Conflict 3") #There are 177 conflict 7 cells in Spain
Conflictcells2 <- subset(PVint2, extendedconflicts == "Conflict 2") #There are 165 conflict 7 cells in Spain
Conflictcells1 <- subset(PVint2, extendedconflicts == "Conflict 1") #There are 129 conflict 7 cells in Spain

# Spatial intersect between ACs and hotspot cells
conflict_AC_int9 <- st_intersection(Conflictcells9, AC)
conflict_AC_int8 <- st_intersection(Conflictcells8, AC)
conflict_AC_int7 <- st_intersection(Conflictcells7, AC)
conflict_AC_int6 <- st_intersection(Conflictcells6, AC)
conflict_AC_int5 <- st_intersection(Conflictcells5, AC)
conflict_AC_int4 <- st_intersection(Conflictcells4, AC)
conflict_AC_int3 <- st_intersection(Conflictcells3, AC)
conflict_AC_int2 <- st_intersection(Conflictcells2, AC)
conflict_AC_int1 <- st_intersection(Conflictcells1, AC)

#Calculating the area of each intersection section
conflict_AC_int9$areaintconf_AC9 <- st_area(conflict_AC_int9)
conflict_AC_int8$areaintconf_AC8 <- st_area(conflict_AC_int8)
conflict_AC_int7$areaintconf_AC7 <- st_area(conflict_AC_int7)
conflict_AC_int6$areaintconf_AC6 <- st_area(conflict_AC_int6)
conflict_AC_int5$areaintconf_AC5 <- st_area(conflict_AC_int5)
conflict_AC_int4$areaintconf_AC4 <- st_area(conflict_AC_int4)
conflict_AC_int3$areaintconf_AC3 <- st_area(conflict_AC_int3)
conflict_AC_int2$areaintconf_AC2 <- st_area(conflict_AC_int2)
conflict_AC_int1$areaintconf_AC1 <- st_area(conflict_AC_int1)

conflict_perAC9 <- conflict_AC_int9 %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

conflict_perAC8 <- conflict_AC_int8 %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

conflict_perAC7 <- conflict_AC_int7 %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

conflict_perAC6 <- conflict_AC_int6 %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

conflict_perAC5 <- conflict_AC_int5 %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

conflict_perAC4 <- conflict_AC_int4 %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())
conflict_perAC3 <- conflict_AC_int3 %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

conflict_perAC2 <- conflict_AC_int2 %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

conflict_perAC1 <- conflict_AC_int1 %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

names(conflict_perAC9)[2]="N conflict cells 9"
names(conflict_perAC8)[2]="N conflict cells 8"
names(conflict_perAC7)[2]="N conflict cells 7"
names(conflict_perAC6)[2]="N conflict cells 6"
names(conflict_perAC5)[2]="N conflict cells 5"
names(conflict_perAC4)[2]="N conflict cells 4"
names(conflict_perAC3)[2]="N conflict cells 3"
names(conflict_perAC2)[2]="N conflict cells 2"
names(conflict_perAC1)[2]="N conflict cells 1"

#Percentage of conflict cells in relation to the number of cells per AC
AC_confstats9 <- merge(AC_stats, conflict_perAC9, by="ccaa.shortname.en", all.x=TRUE)
AC_confstats8 <- merge(AC_stats, conflict_perAC8, by="ccaa.shortname.en", all.x=TRUE)
AC_confstats7 <- merge(AC_stats, conflict_perAC7, by="ccaa.shortname.en", all.x=TRUE)
AC_confstats6 <- merge(AC_stats, conflict_perAC6, by="ccaa.shortname.en", all.x=TRUE)
AC_confstats5 <- merge(AC_stats, conflict_perAC5, by="ccaa.shortname.en", all.x=TRUE)
AC_confstats4 <- merge(AC_stats, conflict_perAC4, by="ccaa.shortname.en", all.x=TRUE)
AC_confstats3 <- merge(AC_stats, conflict_perAC3, by="ccaa.shortname.en", all.x=TRUE)
AC_confstats2 <- merge(AC_stats, conflict_perAC2, by="ccaa.shortname.en", all.x=TRUE)
AC_confstats1 <- merge(AC_stats, conflict_perAC1, by="ccaa.shortname.en", all.x=TRUE)

#Percentage of conflict cells in relation to the total number of each categoryconflict cells in Spain
AC_confstats9$Perc_conf_Spain9 <- AC_confstats9$`N conflict cells 9`*100/55
AC_confstats8$Perc_conf_Spain8 <- AC_confstats8$`N conflict cells 8`*100/46
AC_confstats7$Perc_conf_Spain7 <- AC_confstats7$`N conflict cells 7`*100/41
AC_confstats6$Perc_conf_Spain6 <- AC_confstats6$`N conflict cells 6`*100/141
AC_confstats5$Perc_conf_Spain5 <- AC_confstats5$`N conflict cells 5`*100/118
AC_confstats4$Perc_conf_Spain4 <- AC_confstats4$`N conflict cells 4`*100/109
AC_confstats3$Perc_conf_Spain3 <- AC_confstats3$`N conflict cells 3`*100/177
AC_confstats2$Perc_conf_Spain2 <- AC_confstats2$`N conflict cells 2`*100/165
AC_confstats1$Perc_conf_Spain1 <- AC_confstats1$`N conflict cells 1`*100/129

AC_confstats9$geometry.x <- NULL
AC_confstats9$geometry.y <- NULL
AC_confstats8$geometry.x <- NULL
AC_confstats8$geometry.y <- NULL
AC_confstats7$geometry.x <- NULL
AC_confstats7$geometry.y <- NULL
AC_confstats6$geometry.x <- NULL
AC_confstats6$geometry.y <- NULL
AC_confstats5$geometry.x <- NULL
AC_confstats5$geometry.y <- NULL
AC_confstats4$geometry.x <- NULL
AC_confstats4$geometry.y <- NULL
AC_confstats3$geometry.x <- NULL
AC_confstats3$geometry.y <- NULL
AC_confstats2$geometry.x <- NULL
AC_confstats2$geometry.y <- NULL
AC_confstats1$geometry.x <- NULL
AC_confstats1$geometry.y <- NULL

#Merging 1 to 9 conflict levels
AC_confstats <- list(AC_confstats1, 
                     AC_confstats2 %>% select(ccaa.shortname.en, `N conflict cells 2`, `Perc_conf_Spain2`), 
                     AC_confstats3 %>% select(ccaa.shortname.en, `N conflict cells 3`, `Perc_conf_Spain3`), 
                     AC_confstats4 %>% select(ccaa.shortname.en, `N conflict cells 4`, `Perc_conf_Spain4`), 
                     AC_confstats5 %>% select(ccaa.shortname.en, `N conflict cells 5`, `Perc_conf_Spain5`), 
                     AC_confstats6 %>% select(ccaa.shortname.en, `N conflict cells 6`, `Perc_conf_Spain6`), 
                     AC_confstats7 %>% select(ccaa.shortname.en, `N conflict cells 7`, `Perc_conf_Spain7`), 
                     AC_confstats8 %>% select(ccaa.shortname.en, `N conflict cells 8`, `Perc_conf_Spain8`),
                     AC_confstats9 %>% select(ccaa.shortname.en, `N conflict cells 9`, `Perc_conf_Spain9`))

AC_confstats <- reduce(AC_confstats, function(x, y) full_join(x, y, by = "ccaa.shortname.en"))

# No-go areas

#Selecting only no-go cells
Veryhigh_nogocells <- subset(PVint2, nogo_risk == "no go: very high risk") #There are 98 high no-go cells in Spain
High_nogocells <- subset(PVint2, nogo_risk == "no go: high risk") #There are 374 moderate no-go cells in Spain
moderate_nogocells <- subset(PVint2, nogo_risk == "no go: moderate risk") #There are 566 lower no-go cells in Spain

# Spatial intersect between ACs and no-go cells
nogo_AC_int_veryhigh <- st_intersection(Veryhigh_nogocells, AC)
#Calculating the area of each intersection section
nogo_AC_int_veryhigh$areaintnogo_AC_vh <- st_area(nogo_AC_int_veryhigh)

# Spatial intersect between ACs and no-go cells
nogo_AC_int_high <- st_intersection(High_nogocells, AC)
#Calculating the area of each intersection section
nogo_AC_int_high$areaintnogo_AC_h <- st_area(nogo_AC_int_high)

# Spatial intersect between ACs and no-go cells
nogo_AC_int_moderate <- st_intersection(moderate_nogocells, AC)
#Calculating the area of each intersection section
nogo_AC_int_moderate$areaintnogo_AC_m <- st_area(nogo_AC_int_moderate)

nogo_perAC_vh <- nogo_AC_int_veryhigh %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

names(nogo_perAC_vh)[2]="very high no-go cells"

nogo_perAC_h <- nogo_AC_int_high %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

names(nogo_perAC_h)[2]="high no-go cells"

nogo_perAC_m <- nogo_AC_int_moderate %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

names(nogo_perAC_m)[2]="moderate no-go cells"

#Percentage of no-go cells in relation to the number of cells per AC and number of cels in Spain
ACnogo_stats3_1 <- merge(AC_confstats9, nogo_perAC_vh, by="ccaa.shortname.en", all.x=TRUE)
ACnogo_stats3 <- ACnogo_stats3_1 %>%
  left_join(cells_perAC %>% select(ccaa.shortname.en, `N cells`), by = "ccaa.shortname.en") %>%
  mutate(`N cells` = coalesce(`N cells.x`,`N cells.y`)) %>%
  select(-`N cells.x`, -`N cells.y`)
ACnogo_stats3$Perc_nogo_AC_vh <- ACnogo_stats3$`very high no-go cells`*100/ACnogo_stats3$`N cells`
ACnogo_stats3$Perc_nogo_Spain_vh <- ACnogo_stats3$`very high no-go cells`*100/98

ACnogo_stats3$geometry <- NULL

ACnogo_stats2_1 <- merge(AC_confstats8, nogo_perAC_h, by="ccaa.shortname.en", all.y=TRUE)
ACnogo_stats2 <- ACnogo_stats2_1 %>%
  left_join(cells_perAC %>% select(ccaa.shortname.en, `N cells`), by = "ccaa.shortname.en") %>%
  mutate(`N cells` = coalesce(`N cells.x`,`N cells.y`)) %>%
  select(-`N cells.x`, -`N cells.y`)
ACnogo_stats2$Perc_nogo_AC_h <- ACnogo_stats2$`high no-go cells`*100/ACnogo_stats2$`N cells`
ACnogo_stats2$Perc_nogo_Spain_h <- ACnogo_stats2$`high no-go cells`*100/374

ACnogo_stats2$geometry <- NULL

ACnogo_stats1_1 <- merge(AC_confstats7, nogo_perAC_m, by="ccaa.shortname.en", all.y=TRUE)
ACnogo_stats1 <- ACnogo_stats1_1 %>%
  left_join(cells_perAC %>% select(ccaa.shortname.en, `N cells`), by = "ccaa.shortname.en") %>%
  mutate(`N cells` = coalesce(`N cells.x`,`N cells.y`)) %>%
  select(-`N cells.x`, -`N cells.y`)
ACnogo_stats1$Perc_nogo_AC_m <- ACnogo_stats1$`moderate no-go cells`*100/ACnogo_stats1$`N cells`
ACnogo_stats1$Perc_nogo_Spain_m <- ACnogo_stats1$`moderate no-go cells`*100/566

ACnogo_stats1$geometry <- NULL

#Merging 1 to 9 conflict levels
ACnogo_stats <- list(ACnogo_stats1 %>% select(ccaa.shortname.en, `moderate no-go cells`, `Perc_nogo_AC_m`, `Perc_nogo_Spain_m`), 
                     ACnogo_stats2 %>% select(ccaa.shortname.en, `high no-go cells`, `Perc_nogo_AC_h`, `Perc_nogo_Spain_h`), 
                     ACnogo_stats3 %>% select(ccaa.shortname.en, `very high no-go cells`, `Perc_nogo_AC_vh`, `Perc_nogo_Spain_vh`))

ACnogo_stats <- reduce(ACnogo_stats, function(x, y) full_join(x, y, by = "ccaa.shortname.en"))

Finalstats  <- merge(AC_confstats, ACnogo_stats, by="ccaa.shortname.en")

# Merging two columns in one to show N hotspot/two facet and one facet cells and percentage in relation to the number of cells in each AC
Finalstats$`N hotspot cells` <- paste(Finalstats$`N hotspot cells`, "[",round(Finalstats$Perc_hotp_AC,1),"%]", sep ="")
Finalstats$`N twofacet cells` <- paste(Finalstats$`N twofacet cells`, "[",round(Finalstats$Perc_twof_AC,1),"%]", sep ="")
Finalstats$`N onefacet cells` <- paste(Finalstats$`N onefacet cells`, "[",round(Finalstats$Perc_onef_AC,1),"%]", sep ="")

Finalstats$Perc_hotp_AC <- NULL
Finalstats$Perc_twof_AC <- NULL
Finalstats$Perc_onef_AC <- NULL

# Merging two columns in one to show N no go cells and percentage in relation to the number of cells in each AC
Finalstats$`moderate no-go cells` <- paste(Finalstats$`moderate no-go cells`, "[",round(Finalstats$Perc_nogo_AC_m,1),"%]", sep ="")
Finalstats$`high no-go cells` <- paste(Finalstats$`high no-go cells`, "[",round(Finalstats$Perc_nogo_AC_h,1),"%]", sep ="")
Finalstats$`very high no-go cells` <- paste(Finalstats$`very high no-go cells`, "[",round(Finalstats$Perc_nogo_AC_vh,1),"%]", sep ="")

Finalstats$Perc_nogo_AC_m <- NULL
Finalstats$Perc_nogo_AC_h <- NULL
Finalstats$Perc_nogo_AC_vh <- NULL

# Calculating the % very high conflict cells (7, 8 and 9 conflict levels) in relation to the number of cells in each AC
Finalstats$`Sum veryhigh conflict cells` <- rowSums(Finalstats[, c("N conflict cells 7", "N conflict cells 8", "N conflict cells 9")], na.rm = TRUE)
Finalstats$Perc_veryhighconf <- (Finalstats$`Sum veryhigh conflict cells` / Finalstats$`N cells`) * 100
Finalstats$`N veryhigh conflict cells` <- paste(Finalstats$`Sum veryhigh conflict cells`, "[",round(Finalstats$Perc_veryhighconf,1),"%]", sep ="")
Finalstats$veryhigh_conf_Spain <- (Finalstats$`Sum veryhigh conflict cells` / 142) * 100

# Calculating the % high conflict cells (4,5, and 6 conflict levels) in relation to the number of cells in each AC
Finalstats$`Sum high conflict cells` <- rowSums(Finalstats[, c("N conflict cells 4", "N conflict cells 5", "N conflict cells 6")], na.rm = TRUE)
Finalstats$Perc_highconf <- (Finalstats$`Sum high conflict cells` / Finalstats$`N cells`) * 100
Finalstats$`N high conflict cells` <- paste(Finalstats$`Sum high conflict cells`, "[",round(Finalstats$Perc_highconf,1),"%]", sep ="")
Finalstats$High_conf_Spain <- (Finalstats$`Sum high conflict cells` / 368) * 100

# Calculating the % moderate conflict cells (1,2, and 3 conflict levels) in relation to the number of cells in each AC
Finalstats$`Sum moderate conflict cells` <- rowSums(Finalstats[, c("N conflict cells 1", "N conflict cells 2", "N conflict cells 3")], na.rm = TRUE)
Finalstats$Perc_moderateconf <- (Finalstats$`Sum moderate conflict cells` / Finalstats$`N cells`) * 100
Finalstats$`N moderate conflict cells` <- paste(Finalstats$`Sum moderate conflict cells`, "[",round(Finalstats$Perc_moderateconf,1),"%]", sep ="")
Finalstats$moderate_conf_Spain <- (Finalstats$`Sum moderate conflict cells` / 471) * 100

Finalstats$Perc_conf_Spain1<- NULL
Finalstats$Perc_conf_Spain2<- NULL
Finalstats$Perc_conf_Spain3<- NULL
Finalstats$Perc_conf_Spain4<- NULL
Finalstats$Perc_conf_Spain5<- NULL
Finalstats$Perc_conf_Spain6<- NULL
Finalstats$Perc_conf_Spain7<- NULL
Finalstats$Perc_conf_Spain8<- NULL
Finalstats$Perc_conf_Spain9<- NULL
Finalstats$`Sum veryhigh conflict cells`<- NULL
Finalstats$Perc_veryhighconf<- NULL
Finalstats$`Sum high conflict cells`<- NULL
Finalstats$Perc_highconf<- NULL
Finalstats$`Sum moderate conflict cells`<- NULL
Finalstats$Perc_moderateconf<- NULL

write_xlsx(Finalstats, 'Data/AC_stats_6November.xlsx') 

#Now I want to scale TD, FD, and PD in a range from 0-1 TO USE ZONATION SOFTWARE. Also 
basic_function <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}
PVint2_1<-PVint2
PVint2_1$SESFRic <- basic_function(PVint2_1$SESFRic)
PVint2_1$SESPD <- basic_function(PVint2_1$SESPD)
PVint2_1$species <- basic_function(PVint2_1$species)
# Replacing NAs by 0 in PV occupancy cells
PVint2_1$PercPV_cell[is.na(PVint2_1$PercPV_cell)] <- 0
#Now also scaling PV occupancy from 0-1 TO USE ZONATION SOFTWARE. However, maximum occupancy % per cell is around 12%
PVint2_1$PercPV_cell <- basic_function(PVint2_1$PercPV_cell)

write_sf(PVint2_1, "Spatial_Data/3facets_scaled.shp")

#Using maps generated by Zonation 5
#First using the map tht combines TD, FD, andPD
Hotsp_Zonation <- raster("Spatial_Data/Zonation/rankmap_hotsp.tif")

Hotsp_percentil <- quantile(Hotsp_Zonation, 0.7)

# Creating a binary raster wth 1 for cells where the value is the same or greater than percentil 70
Hotsp_top30 <- calc(Hotsp_Zonation, fun = function(x) { ifelse(x >= Hotsp_percentil, 1, 0) })
Hotsp_top30_count <- freq(Hotsp_top30)

print(Hotsp_top30_count)
# Guardar el mapa binario como PNG sin título, ejes ni leyenda
png("hotsp_top30_map_clean.png", width = 800, height = 800)

# Plotear el mapa binario sin ejes, sin título, sin marco y sin leyenda
plot(Hotsp_top30, axes = FALSE, box = FALSE, main = "", legend = FALSE, col = c("white", "red"))

# Cerrar el dispositivo gráfico
dev.off()


hotsp_top30perAC <- extract(Hotsp_top30, AC, fun = sum, na.rm = TRUE)

# nEW FIELD TO HAVE THE NUMBER OF TOP30 HOTSPOT CELLS PER AC
AC$top30_hotsp <- hotsp_top30perAC

#Now the map tht combines TD, FD
TDFD_Zonation <- raster("Spatial_Data/Zonation/rankmap_TDFD.tif")

TDFD_percentil <- quantile(TDFD_Zonation, 0.7)

# Creating a binary raster wth 1 for cells where the value is the same or greater than percentil 70
TDFD_top30 <- calc(TDFD_Zonation, fun = function(x) { ifelse(x >= TDFD_percentil, 1, 0) })
TDFD_top30_count <- freq(TDFD_top30)

print(TDFD_top30_count)
plot(TDFD_top30)

TDFD_top30perAC <- extract(TDFD_top30, AC, fun = sum, na.rm = TRUE)

# nEW FIELD TO HAVE THE NUMBER OF TOP30 TDFD CELLS PER AC
AC$top30_TDFD <- TDFD_top30perAC

#Now the map tht combines TD, PD
TDPD_Zonation <- raster("Spatial_Data/Zonation/rankmap_TDPD.tif")

TDPD_percentil <- quantile(TDPD_Zonation, 0.7)

# Creating a binary raster wth 1 for cells where the value is the same or greater than percentil 70
TDPD_top30 <- calc(TDPD_Zonation, fun = function(x) { ifelse(x >= TDPD_percentil, 1, 0) })

plot(TDPD_top30)

TDPD_top30perAC <- extract(TDPD_top30, AC, fun = sum, na.rm = TRUE)

# nEW FIELD TO HAVE THE NUMBER OF TOP30 TDPD CELLS PER AC
AC$top30_TDPD <- TDPD_top30perAC

#Now the map tht combines FD, PD
FDPD_Zonation <- raster("Spatial_Data/Zonation/rankmap_FDPD.tif")

FDPD_percentil <- quantile(FDPD_Zonation, 0.7)

# Creating a binary raster wth 1 for cells where the value is the same or greater than percentil 70
FDPD_top30 <- calc(FDPD_Zonation, fun = function(x) { ifelse(x >= FDPD_percentil, 1, 0) })

plot(FDPD_top30)

FDPD_top30perAC <- extract(FDPD_top30, AC, fun = sum, na.rm = TRUE)

# nEW FIELD TO HAVE THE NUMBER OF TOP30 FDPD CELLS PER AC
AC$top30_FDPD <- FDPD_top30perAC

#Now the map that contains TD
TD_Zonation <- raster("Spatial_Data/Zonation/rankmap_TD.tif")

TD_percentil <- quantile(TD_Zonation, 0.7)

# Creating a binary raster wth 1 for cells where the value is the same or greater than percentil 70
TD_top30 <- calc(TD_Zonation, fun = function(x) { ifelse(x >= TD_percentil, 1, 0) })

# Guardar el mapa binario como PNG sin título, ejes ni leyenda
png("TD_top30_map_clean.png", width = 800, height = 800)

# Plotear el mapa binario sin ejes, sin título, sin marco y sin leyenda
plot(TD_top30, axes = FALSE, box = FALSE, main = "", legend = FALSE, col = c("white", "turquoise"))

# Cerrar el dispositivo gráfico
dev.off()

TD_top30perAC <- extract(TD_top30, AC, fun = sum, na.rm = TRUE)

# nEW FIELD TO HAVE THE NUMBER OF TOP30 TD CELLS PER AC
AC$top30_TD <- TD_top30perAC

#Now the map that contains FD
FD_Zonation <- raster("Spatial_Data/Zonation/rankmap_FD.tif")

FD_percentil <- quantile(FD_Zonation, 0.7)

# Creating a binary raster wth 1 for cells where the value is the same or greater than percentil 70
FD_top30 <- calc(FD_Zonation, fun = function(x) { ifelse(x >= FD_percentil, 1, 0) })
FD_top30_count <- freq(FD_top30)

print(FD_top30_count)

# Guardar el mapa binario como PNG sin título, ejes ni leyenda
png("FD_top30_map_clean.png", width = 800, height = 800)

# Plotear el mapa binario sin ejes, sin título, sin marco y sin leyenda
plot(FD_top30, axes = FALSE, box = FALSE, main = "", legend = FALSE, col = c("white", "paleturquoise3"))

# Cerrar el dispositivo gráfico
dev.off()

FD_top30perAC <- extract(FD_top30, AC, fun = sum, na.rm = TRUE)

# nEW FIELD TO HAVE THE NUMBER OF TOP30 FD CELLS PER AC
AC$top30_FD <- FD_top30perAC

#Now the map that contains PD
PD_Zonation <- raster("Spatial_Data/Zonation/rankmap_PD.tif")

PD_percentil <- quantile(PD_Zonation, 0.7)

# Creating a binary raster wth 1 for cells where the value is the same or greater than percentil 70
PD_top30 <- calc(PD_Zonation, fun = function(x) { ifelse(x >= PD_percentil, 1, 0) })
PD_top30_count <- freq(PD_top30)

print(PD_top30_count)
# Guardar el mapa binario como PNG sin título, ejes ni leyenda
png("PD_top30_map_clean.png", width = 800, height = 800)

# Plotear el mapa binario sin ejes, sin título, sin marco y sin leyenda
plot(PD_top30, axes = FALSE, box = FALSE, main = "", legend = FALSE, col = c("white", "grey"))

# Cerrar el dispositivo gráfico
dev.off()

PD_top30perAC <- extract(PD_top30, AC, fun = sum, na.rm = TRUE)

# nEW FIELD TO HAVE THE NUMBER OF TOP30 PD CELLS PER AC
AC$top30_PD <- PD_top30perAC


# Crear los mapas binarios
TD_top30 <- calc(TD_Zonation, fun = function(x) { ifelse(x >= TD_percentil, 1, 0) })
FD_top30 <- calc(FD_Zonation, fun = function(x) { ifelse(x >= FD_percentil, 1, 0) })
PD_top30 <- calc(PD_Zonation, fun = function(x) { ifelse(x >= PD_percentil, 1, 0) })
Hotsp_top30 <- calc(Hotsp_Zonation, fun = function(x) { ifelse(x >= Hotsp_percentil, 1, 0) })

# Convertir los raster a data frames para ggplot2
TD_df <- as.data.frame(as(TD_top30, "SpatialPixelsDataFrame"))
FD_df <- as.data.frame(as(FD_top30, "SpatialPixelsDataFrame"))
PD_df <- as.data.frame(as(PD_top30, "SpatialPixelsDataFrame"))
Hotsp_df <- as.data.frame(as(Hotsp_top30, "SpatialPixelsDataFrame"))

# Asegurar que los valores sean tratados como factores
TD_df$layer <- factor(TD_df$layer)
FD_df$layer <- factor(FD_df$layer)
PD_df$layer <- factor(PD_df$layer)
Hotsp_df$layer <- factor(Hotsp_df$layer)

# Crear un plot ggplot2 para cada mapa
p1 <- ggplot() +
  geom_tile(data = TD_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_manual(values = c("white", "turquoise"), guide = "none") +
  geom_sf(data = AC, fill = NA, color = "black") +
  ggtitle("Prioritization TD") +
  theme_void() +
  theme(
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)  # Centrar el título
  )

p2 <- ggplot() +
  geom_tile(data = FD_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_manual(values = c("white", "paleturquoise3"), guide = "none") +
  geom_sf(data = AC, fill = NA, color = "black") +
  ggtitle("Prioritization FD") +
  theme_void() +
  theme(
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)  # Centrar el título
  )

p3 <- ggplot() +
  geom_tile(data = PD_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_manual(values = c("white", "grey"), guide = "none") +
  geom_sf(data = AC, fill = NA, color = "black") +
  ggtitle("Prioritization PD") +
  theme_void() +
  theme(
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)  # Centrar el título
  )

p4 <- ggplot() +
  geom_tile(data = Hotsp_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_manual(values = c("white", "red"), guide = "none") +
  geom_sf(data = AC, fill = NA, color = "black") +
  ggtitle("Prioritization TD-FD-PD") +
  theme_void() +
  theme(
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)  # Centrar el título
  ) +
  annotation_scale(location = "br", width_hint = 0.5, height = unit(0.8, "cm"), 
                   text_cex = 2, bar_cols = c("black", "white"), pad_x = unit(0.5, "cm"), 
                   pad_y = unit(0.5, "cm"), dist = 100, dist_unit = "km")

# Combinar los cuatro plots en una sola imagen
png("four_maps.png", width = 1600, height = 1600)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

pdf("four_maps.pdf", width = 16, height = 16) 
grid.arrange(p1, p2, p3, p4, nrow = 2) 
dev.off ()














AC_tops <- AC[, c("ccaa.shortname.en", "top30_hotsp", "top30_TDFD", "top30_TDPD", "top30_FDPD", "top30_TD", "top30_FD", "top30_PD")]

AC_int1 <- st_intersection(datosMAP, AC)

cells_perAC1 <- AC_int1 %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())
AC_tops$geometry <- NULL

AC_tops<-merge(AC_tops, cells_perAC1, by="ccaa.shortname.en")
AC_tops$geometry <- NULL

row_sum <- colSums(AC_tops[2:9])
new_row <- data.frame(ccaa.shortname.en = "Total", t(as.data.frame(row_sum)))
names(new_row) <- names(AC_tops)
AC_tops <- rbind(AC_tops, new_row)

write.csv(AC_tops, "Spatial_Data/Zonation/AC_tops1.csv", row.names = FALSE)

#Calculating overlaps among the three facets
# Combinación de las capas binarias
combined_raster <- (TD_top30 * 1) + (FD_top30 * 2) + (PD_top30 * 4)

# Visualizar el raster combinado
categories <- c("TD", "FD", "PD", "TD-FD", "TD-PD", "FD-PD", "TD-FD-PD")
# Colores específicos para cada categoría
colores <- c("red", "green", "blue", "yellow", "orange", "purple", "black")

# Crear un mapa del solapamiento
plot(combined_raster, main = "Diversity Overlap", col = colores, legend = FALSE)

# Añadir la leyenda
legend("topright", legend = categories, fill = colores, title = "Categories")

# Calcular frecuencias
frequencies <- freq(combined_raster, useNA = "no")
frequencies_df <- as.data.frame(frequencies)
colnames(frequencies_df) <- c("Value", "Cell_Count")

# Filtrar valores no deseados (si es necesario)
frequencies_df <- frequencies_df[frequencies_df$Value != 0, ]

# Asignar etiquetas de categorías a los valores combinados
# Usamos un mapeo explícito para asegurarnos de que las combinaciones se asignen correctamente
category_map <- c(
  "TD" = 1, 
  "FD" = 2, 
  "PD" = 4, 
  "TD-FD" = 3, 
  "TD-PD" = 5, 
  "PD-FD" = 6, 
  "TD-FD-PD" = 7
)

# Asignar categorías basadas en el valor del raster combinado
frequencies_df$Category <- names(category_map)[match(frequencies_df$Value, category_map)]

# Revisar el resultado
print(frequencies_df)

#I want to know the percentage of each facet (individual facet scenarios) that is covered under all scenarios

# Calcular porcentaje de cobertura para TD
TD_overlap <- TD_top30 * TD_top30  # Crear un raster de intersección (1 = solapamiento)
TD_coverage <- cellStats(TD_overlap, sum) / cellStats(TD_top30, sum) * 100

# Calcular porcentaje de cobertura para FD
FD_overlap <- FD_top30 * TD_top30
FD_coverage <- cellStats(FD_overlap, sum) / cellStats(FD_top30, sum) * 100

# Calcular porcentaje de cobertura para PD
PD_overlap <- PD_top30 * TD_top30
PD_coverage <- cellStats(PD_overlap, sum) / cellStats(PD_top30, sum) * 100

# Imprimir resultados
cat("Cobertura relativa de TD:", TD_coverage, "%\n")
cat("Cobertura relativa de FD:", FD_coverage, "%\n")
cat("Cobertura relativa de PD:", PD_coverage, "%\n")

# Calcular porcentaje de cobertura para TD
TD_overlap <- TD_top30 * FD_top30  # Crear un raster de intersección (1 = solapamiento)
TD_coverage <- cellStats(TD_overlap, sum) / cellStats(TD_top30, sum) * 100

# Calcular porcentaje de cobertura para FD
FD_overlap <- FD_top30 * FD_top30
FD_coverage <- cellStats(FD_overlap, sum) / cellStats(FD_top30, sum) * 100

# Calcular porcentaje de cobertura para PD
PD_overlap <- PD_top30 * FD_top30
PD_coverage <- cellStats(PD_overlap, sum) / cellStats(PD_top30, sum) * 100

# Imprimir resultados
cat("Cobertura relativa de TD:", TD_coverage, "%\n")
cat("Cobertura relativa de FD:", FD_coverage, "%\n")
cat("Cobertura relativa de PD:", PD_coverage, "%\n")

# Calcular porcentaje de cobertura para TD
TD_overlap <- TD_top30 * PD_top30  # Crear un raster de intersección (1 = solapamiento)
TD_coverage <- cellStats(TD_overlap, sum) / cellStats(TD_top30, sum) * 100

# Calcular porcentaje de cobertura para FD
FD_overlap <- FD_top30 * PD_top30
FD_coverage <- cellStats(FD_overlap, sum) / cellStats(FD_top30, sum) * 100

# Calcular porcentaje de cobertura para PD
PD_overlap <- PD_top30 * PD_top30
PD_coverage <- cellStats(PD_overlap, sum) / cellStats(PD_top30, sum) * 100

# Imprimir resultados
cat("Cobertura relativa de TD:", TD_coverage, "%\n")
cat("Cobertura relativa de FD:", FD_coverage, "%\n")
cat("Cobertura relativa de PD:", PD_coverage, "%\n")

# Calcular porcentaje de cobertura para TD
TD_overlap <- TD_top30 * TDFD_top30  # Crear un raster de intersección (1 = solapamiento)
TD_coverage <- cellStats(TD_overlap, sum) / cellStats(TD_top30, sum) * 100

# Calcular porcentaje de cobertura para FD
FD_overlap <- FD_top30 * TDFD_top30
FD_coverage <- cellStats(FD_overlap, sum) / cellStats(FD_top30, sum) * 100

# Calcular porcentaje de cobertura para PD
PD_overlap <- PD_top30 * TDFD_top30
PD_coverage <- cellStats(PD_overlap, sum) / cellStats(PD_top30, sum) * 100

# Imprimir resultados
cat("Cobertura relativa de TD:", TD_coverage, "%\n")
cat("Cobertura relativa de FD:", FD_coverage, "%\n")
cat("Cobertura relativa de PD:", PD_coverage, "%\n")

# Calcular porcentaje de cobertura para TD
TD_overlap <- TD_top30 * TDPD_top30  # Crear un raster de intersección (1 = solapamiento)
TD_coverage <- cellStats(TD_overlap, sum) / cellStats(TD_top30, sum) * 100

# Calcular porcentaje de cobertura para FD
FD_overlap <- FD_top30 * TDPD_top30
FD_coverage <- cellStats(FD_overlap, sum) / cellStats(FD_top30, sum) * 100

# Calcular porcentaje de cobertura para PD
PD_overlap <- PD_top30 * TDPD_top30
PD_coverage <- cellStats(PD_overlap, sum) / cellStats(PD_top30, sum) * 100

# Imprimir resultados
cat("Cobertura relativa de TD:", TD_coverage, "%\n")
cat("Cobertura relativa de FD:", FD_coverage, "%\n")
cat("Cobertura relativa de PD:", PD_coverage, "%\n")

# Calcular porcentaje de cobertura para TD
TD_overlap <- TD_top30 * FDPD_top30  # Crear un raster de intersección (1 = solapamiento)
TD_coverage <- cellStats(TD_overlap, sum) / cellStats(TD_top30, sum) * 100

# Calcular porcentaje de cobertura para FD
FD_overlap <- FD_top30 * FDPD_top30
FD_coverage <- cellStats(FD_overlap, sum) / cellStats(FD_top30, sum) * 100

# Calcular porcentaje de cobertura para PD
PD_overlap <- PD_top30 * FDPD_top30
PD_coverage <- cellStats(PD_overlap, sum) / cellStats(PD_top30, sum) * 100

# Imprimir resultados
cat("Cobertura relativa de TD:", TD_coverage, "%\n")
cat("Cobertura relativa de FD:", FD_coverage, "%\n")
cat("Cobertura relativa de PD:", PD_coverage, "%\n")

# Crear mapas de solapamiento
par(mfrow = c(1, 1))  # Crear una cuadrícula de 1 fila y 3 columnas para los mapas

# Mapas de solapamiento para cada faceta
plot(TD_overlap, main = "Solapamiento TD con Hotspots", col = c("white", "green"), legend = FALSE)
plot(FD_overlap, main = "Solapamiento FD con Hotspots", col = c("white", "blue"), legend = FALSE)
plot(PD_overlap, main = "Solapamiento PD con Hotspots", col = c("white", "red"), legend = FALSE)

plot(TD_top30, main = "TD top30", col = c("white", "green"), legend = FALSE)
plot(FD_top30, main = "FD top30", col = c("white", "blue"), legend = FALSE)
plot(PD_top30, main = "PD top30", col = c("white", "red"), legend = FALSE)

# Asegúrate de que el shapefile esté cargado como un objeto sf
facets <- st_read("Spatial_Data/3facets_scaled.shp")

# Crear una función para clasificar y mapear los valores
plot_diversity <- function(data, column, title, borders, show_legend = TRUE, show_annotations = TRUE) {
  # Ordenar y calcular los cuartiles personalizados
  data <- data %>%
    mutate(
      Rank = rank(-!!sym(column)),  # Rango descendente
      Percentile = Rank / nrow(data) * 100  # Percentil
    ) %>%
    mutate(
      Category = case_when(
        Percentile <= 5 ~ "Top 5%",
        Percentile <= 10 ~ "Top 10%",
        Percentile <= 17 ~ "Top 17%",
        Percentile <= 30 ~ "Top 30%",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Category)) %>%
    mutate(Category = factor(Category, levels = c("Top 5%", "Top 10%", "Top 17%", "Top 30%")))  # Ordenar los niveles
  
  p <- ggplot() +
    geom_sf(data = data, aes(fill = Category), color = NA, size = 0) +  # Eliminar bordes blancos
    geom_sf(data = borders, fill = NA, color = "black", size = 0.3) +  # Bordes de comunidades autónomas
    scale_fill_manual(
      values = c("Top 5%" = "red", 
                 "Top 10%" = "orange", 
                 "Top 17%" = "cyan3", 
                 "Top 30%" = "ivory4"),
      name = "Category"
    ) +
    labs(title = title) +
    theme_minimal() +
    theme(
      legend.position = if (show_legend) c(0.8, 0.1) else "none",  # Ajusta la posición de la leyenda para dar más espacio
      legend.text = element_text(size = 8),  # Aumenta el tamaño del texto de la leyenda
      legend.title = element_text(size = 8),  # Aumenta el tamaño del título de la leyenda
      plot.title = element_text(size = 12, hjust = 0.5),  # Aumenta el tamaño del título del mapa y centra el título
      plot.margin = margin(5, 5, 5, 5),  # Ajusta los márgenes
      panel.grid.major = element_blank(),  # Elimina las líneas de cuadrícula mayores
      panel.grid.minor = element_blank(),  # Elimina las líneas de cuadrícula menores
      axis.text = element_blank(),  # Elimina los números de latitud y longitud
      axis.ticks = element_blank()  # Elimina las marcas de los ejes
    )
  
  if (show_annotations) {
    p <- p + annotation_scale(location = "bl", width_hint = 0.3)
  }
  
  return(p)
}

# Generar mapas para cada diversidad, asegurando que sean del mismo tamaño y ajustando la visibilidad de la leyenda y anotaciones
map_td <- plot_diversity(facets, "species", "Taxonomic Diversity (TD)", borders = AC, show_legend = FALSE, show_annotations = FALSE)
map_fd <- plot_diversity(facets, "SESFRic", "Functional Diversity (FD)", borders = AC, show_legend = FALSE, show_annotations = FALSE)
map_pd <- plot_diversity(facets, "SESPD", "Phylogenetic Diversity (PD)", borders = AC, show_legend = TRUE, show_annotations = TRUE)

# Combinar los tres mapas en una sola gráfica asegurando proporciones iguales
combined_map1 <- grid.arrange(
  map_td, 
  map_fd, 
  ncol = 2
)

combined_map2 <- grid.arrange(
  map_fd, 
  map_pd,
  ncol = 2
)

# Guardar la gráfica combinada
ggsave("Figures/combined_map1.png", combined_map1, width = 20, height = 10, units = "cm", dpi = 300)
ggsave("Figures/combined_map1.pdf", combined_map1, width = 20, height = 10, units = "cm", dpi = 300)
ggsave("Figures/combined_map2.png", combined_map2, width = 20, height = 10, units = "cm", dpi = 300)
ggsave("Figures/combined_map2.pdf", combined_map2, width = 20, height = 10, units = "cm", dpi = 300)










facets1 <- facets %>% select(UTMCODE, species, SESFRic, SESPD, geometry)

# Crear las nuevas columnas top5, top10, top17, top30 para TD, FD y PD
facets1 <- facets1 %>%
  mutate(
    Rank_species = rank(-species),  # Rango descendente para species
    Percentile_species = Rank_species / n() * 100,  # Percentil para species
    top30TD = ifelse(Percentile_species <= 30, "top30", NA_character_),
    top17TD = ifelse(Percentile_species <= 17, "top17", NA_character_),
    top10TD = ifelse(Percentile_species <= 10, "top10", NA_character_),
    top5TD = ifelse(Percentile_species <= 5, "top5", NA_character_),
    
    Rank_SESFRic = rank(-SESFRic),  # Rango descendente para SESFRic
    Percentile_SESFRic = Rank_SESFRic / n() * 100,  # Percentil para SESFRic
    top30FD = ifelse(Percentile_SESFRic <= 30, "top30", NA_character_),
    top17FD = ifelse(Percentile_SESFRic <= 17, "top17", NA_character_),
    top10FD = ifelse(Percentile_SESFRic <= 10, "top10", NA_character_),
    top5FD = ifelse(Percentile_SESFRic <= 5, "top5", NA_character_),
    
    Rank_SESPD = rank(-SESPD),  # Rango descendente para SESPD
    Percentile_SESPD = Rank_SESPD / n() * 100,  # Percentil para SESPD
    top30PD = ifelse(Percentile_SESPD <= 30, "top30", NA_character_),
    top17PD = ifelse(Percentile_SESPD <= 17, "top17", NA_character_),
    top10PD = ifelse(Percentile_SESPD <= 10, "top10", NA_character_),
    top5PD = ifelse(Percentile_SESPD <= 5, "top5", NA_character_)
  ) %>%
  select(-Rank_species, -Percentile_species, -Rank_SESFRic, -Percentile_SESFRic, -Rank_SESPD, -Percentile_SESPD)  # Eliminar columnas intermedias

TDtop30 <- subset(facets1, top30TD == "top30") 
TDtop30 <- st_union(TDtop30) 

TDtop17 <- subset(facets1, top17TD == "top17") 
TDtop17 <- st_union(TDtop17) 

TDtop10 <- subset(facets1, top10TD == "top10") 
TDtop10 <- st_union(TDtop10) 

TDtop5 <- subset(facets1, top5TD == "top5") 
TDtop5 <- st_union(TDtop5) 

FDtop30 <- subset(facets1, top30FD == "top30") 
FDtop30 <- st_union(FDtop30) 

FDtop17 <- subset(facets1, top17FD == "top17") 
FDtop17 <- st_union(FDtop17) 

FDtop10 <- subset(facets1, top10FD == "top10") 
FDtop10 <- st_union(FDtop10) 

FDtop5 <- subset(facets1, top5FD == "top5") 
FDtop5 <- st_union(FDtop5) 

PDtop30 <- subset(facets1, top30PD == "top30") 
PDtop30 <- st_union(PDtop30) 

PDtop17 <- subset(facets1, top17PD == "top17") 
PDtop17 <- st_union(PDtop17) 

PDtop10 <- subset(facets1, top10PD == "top10") 
PDtop10 <- st_union(PDtop10) 

PDtop5 <- subset(facets1, top5PD == "top5") 
PDtop5 <- st_union(PDtop5) 

# Spatial intersect between ACs and hotspot cells
TDtop30_sf <- st_sf(geometry = TDtop30, dummy = 1)
TDtop30_int <- st_intersection(TDtop30_sf, AC)
TDtop30_int$TD30_areaAC <- st_area(TDtop30_int)
TDtop30_int <- as.data.frame(st_drop_geometry(TDtop30_int))

TDtop17_sf <- st_sf(geometry = TDtop17, dummy = 1)
TDtop17_int <- st_intersection(TDtop17_sf, AC)
TDtop17_int$TD17_areaAC <- st_area(TDtop17_int)
TDtop17_int <- as.data.frame(st_drop_geometry(TDtop17_int))

TDtop10_sf <- st_sf(geometry = TDtop10, dummy = 1)
TDtop10_int <- st_intersection(TDtop10_sf, AC)
TDtop10_int$TD10_areaAC <- st_area(TDtop10_int)
TDtop10_int <- as.data.frame(st_drop_geometry(TDtop10_int))

TDtop5_sf <- st_sf(geometry = TDtop5, dummy = 1)
TDtop5_int <- st_intersection(TDtop5_sf, AC)
TDtop5_int$TD5_areaAC <- st_area(TDtop5_int)
TDtop5_int <- as.data.frame(st_drop_geometry(TDtop5_int))

FDtop30_sf <- st_sf(geometry = FDtop30, dummy = 1)
FDtop30_int <- st_intersection(FDtop30_sf, AC)
FDtop30_int$FD30_areaAC <- st_area(FDtop30_int)
FDtop30_int <- as.data.frame(st_drop_geometry(FDtop30_int))

FDtop17_sf <- st_sf(geometry = FDtop17, dummy = 1)
FDtop17_int <- st_intersection(FDtop17_sf, AC)
FDtop17_int$FD17_areaAC <- st_area(FDtop17_int)
FDtop17_int <- as.data.frame(st_drop_geometry(FDtop17_int))

FDtop10_sf <- st_sf(geometry = FDtop10, dummy = 1)
FDtop10_int <- st_intersection(FDtop10_sf, AC)
FDtop10_int$FD10_areaAC <- st_area(FDtop10_int)
FDtop10_int <- as.data.frame(st_drop_geometry(FDtop10_int))

FDtop5_sf <- st_sf(geometry = FDtop5, dummy = 1)
FDtop5_int <- st_intersection(FDtop5_sf, AC)
FDtop5_int$FD5_areaAC <- st_area(FDtop5_int)
FDtop5_int <- as.data.frame(st_drop_geometry(FDtop5_int))

PDtop30_sf <- st_sf(geometry = PDtop30, dummy = 1)
PDtop30_int <- st_intersection(PDtop30_sf, AC)
PDtop30_int$PD30_areaAC <- st_area(PDtop30_int)
PDtop30_int <- as.data.frame(st_drop_geometry(PDtop30_int))

PDtop17_sf <- st_sf(geometry = PDtop17, dummy = 1)
PDtop17_int <- st_intersection(PDtop17_sf, AC)
PDtop17_int$PD17_areaAC <- st_area(PDtop17_int)
PDtop17_int <- as.data.frame(st_drop_geometry(PDtop17_int))

PDtop10_sf <- st_sf(geometry = PDtop10, dummy = 1)
PDtop10_int <- st_intersection(PDtop10_sf, AC)
PDtop10_int$PD10_areaAC <- st_area(PDtop10_int)
PDtop10_int <- as.data.frame(st_drop_geometry(PDtop10_int))

PDtop5_sf <- st_sf(geometry = PDtop5, dummy = 1)
PDtop5_int <- st_intersection(PDtop5_sf, AC)
PDtop5_int$PD5_areaAC <- st_area(PDtop5_int)
PDtop5_int <- as.data.frame(st_drop_geometry(PDtop5_int))

# TD Categorías
TDtop30_area <- TDtop30_int %>%
  group_by(ccaa.shortname.en) %>%
  summarize(TD30_total_area = sum(TD30_areaAC, na.rm = TRUE))

TDtop17_area <- TDtop17_int %>%
  group_by(ccaa.shortname.en) %>%
  summarize(TD17_total_area = sum(TD17_areaAC, na.rm = TRUE))

TDtop10_area <- TDtop10_int %>%
  group_by(ccaa.shortname.en) %>%
  summarize(TD10_total_area = sum(TD10_areaAC, na.rm = TRUE))

TDtop5_area <- TDtop5_int %>%
  group_by(ccaa.shortname.en) %>%
  summarize(TD5_total_area = sum(TD5_areaAC, na.rm = TRUE))

# FD Categorías
FDtop30_area <- FDtop30_int %>%
  group_by(ccaa.shortname.en) %>%
  summarize(FD30_total_area = sum(FD30_areaAC, na.rm = TRUE))

FDtop17_area <- FDtop17_int %>%
  group_by(ccaa.shortname.en) %>%
  summarize(FD17_total_area = sum(FD17_areaAC, na.rm = TRUE))

FDtop10_area <- FDtop10_int %>%
  group_by(ccaa.shortname.en) %>%
  summarize(FD10_total_area = sum(FD10_areaAC, na.rm = TRUE))

FDtop5_area <- FDtop5_int %>%
  group_by(ccaa.shortname.en) %>%
  summarize(FD5_total_area = sum(FD5_areaAC, na.rm = TRUE))

# PD Categorías
PDtop30_area <- PDtop30_int %>%
  group_by(ccaa.shortname.en) %>%
  summarize(PD30_total_area = sum(PD30_areaAC, na.rm = TRUE))

PDtop17_area <- PDtop17_int %>%
  group_by(ccaa.shortname.en) %>%
  summarize(PD17_total_area = sum(PD17_areaAC, na.rm = TRUE))

PDtop10_area <- PDtop10_int %>%
  group_by(ccaa.shortname.en) %>%
  summarize(PD10_total_area = sum(PD10_areaAC, na.rm = TRUE))

PDtop5_area <- PDtop5_int %>%
  group_by(ccaa.shortname.en) %>%
  summarize(PD5_total_area = sum(PD5_areaAC, na.rm = TRUE))

# Unir resultados en un único data frame
area_results <- list(
  TDtop30_area, TDtop17_area, TDtop10_area, TDtop5_area,
  FDtop30_area, FDtop17_area, FDtop10_area, FDtop5_area,
  PDtop30_area, PDtop17_area, PDtop10_area, PDtop5_area
) %>%
  reduce(full_join, by = "ccaa.shortname.en")

area_results <- drop_units(area_results) #all calculations are in m2


# Calcular el área total de estudio en España (capa facets1)
area_total_espana <- sum(as.numeric(st_area(facets1)))

# Calcular los porcentajes para cada grupo en relación al área total de estudio
TDtop30_porcentaje <- TDtop30_area %>% mutate(TD30_porcentaje = (TD30_total_area * 100) / st_area(TDtop30))
TDtop17_porcentaje <- TDtop17_area %>% mutate(TD17_porcentaje = (TD17_total_area * 100) / st_area(TDtop17))
TDtop10_porcentaje <- TDtop10_area %>% mutate(TD10_porcentaje = (TD10_total_area * 100) / st_area(TDtop10))
TDtop5_porcentaje <- TDtop5_area %>% mutate(TD5_porcentaje = (TD5_total_area * 100) / st_area(TDtop5))

FDtop30_porcentaje <- FDtop30_area %>% mutate(FD30_porcentaje = (FD30_total_area * 100) / st_area(FDtop30))
FDtop17_porcentaje <- FDtop17_area %>% mutate(FD17_porcentaje = (FD17_total_area * 100) / st_area(FDtop17))
FDtop10_porcentaje <- FDtop10_area %>% mutate(FD10_porcentaje = (FD10_total_area * 100) / st_area(FDtop10))
FDtop5_porcentaje <- FDtop5_area %>% mutate(FD5_porcentaje = (FD5_total_area * 100) / st_area(FDtop5))

PDtop30_porcentaje <- PDtop30_area %>% mutate(PD30_porcentaje = (PD30_total_area * 100) / st_area(PDtop30))
PDtop17_porcentaje <- PDtop17_area %>% mutate(PD17_porcentaje = (PD17_total_area * 100) / st_area(PDtop17))
PDtop10_porcentaje <- PDtop10_area %>% mutate(PD10_porcentaje = (PD10_total_area * 100) / st_area(PDtop10))
PDtop5_porcentaje <- PDtop5_area %>% mutate(PD5_porcentaje = (PD5_total_area * 100) / st_area(PDtop5))

# Unir resultados en un único dataframe con áreas y porcentajes
area_y_porcentaje_results <- list(
  TDtop30_porcentaje, TDtop17_porcentaje, TDtop10_porcentaje, TDtop5_porcentaje,
  FDtop30_porcentaje, FDtop17_porcentaje, FDtop10_porcentaje, FDtop5_porcentaje,
  PDtop30_porcentaje, PDtop17_porcentaje, PDtop10_porcentaje, PDtop5_porcentaje
) %>%
  reduce(full_join, by = "ccaa.shortname.en")

area_y_porcentaje_results  <- drop_units(area_y_porcentaje_results)

# Guardar los dataframes combinados en archivos CSV
write.csv(area_y_porcentaje_results, "Data/area_y_porcentaje_results.csv", row.names = FALSE) #percentage of each top in relation to the study area in Spain









ACstats1 <- merge(cellsperAC, FDtop10_AC, by="ccaa.shortname.en", all.x=TRUE)

# Spatial intersect between ACs and malla to know the number of cells per AC
AC_cells <- st_intersection(facets1, AC)

# Spatial intersect between ACs and hotspot cells
TDtop30_int <- st_intersection(TDtop30, AC)
TDtop30_int <- as.data.frame(st_drop_geometry(TDtop30_int))


TDtop17_int <- st_intersection(TDtop17, AC)

TDtop10_int <- st_intersection(TDtop10, AC)

TDtop5_int <- st_intersection(TDtop5, AC)

FDtop30_int <- st_intersection(FDtop30, AC)

FDtop17_int <- st_intersection(FDtop17, AC)

FDtop10_int <- st_intersection(FDtop10, AC)

FDtop5_int <- st_intersection(FDtop5, AC)

PDtop30_int <- st_intersection(PDtop30, AC)

PDtop17_int <- st_intersection(PDtop17, AC)

PDtop10_int <- st_intersection(PDtop10, AC)

PDtop5_int <- st_intersection(PDtop5, AC)


cellsperAC <- AC_cells %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

names(cellsperAC)[2]="N cells"
cellsperAC$geometry <- NULL

TDtop30_AC <- TDtop30_int %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())
names(TDtop30_AC)[2]="N_TDtop30"

TDtop17_AC <- TDtop17_int %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())
names(TDtop17_AC)[2]="N_TDtop17"

TDtop10_AC <- TDtop10_int %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())
names(TDtop10_AC)[2]="N_TDtop10"

TDtop5_AC <- TDtop5_int %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())
names(TDtop5_AC)[2]="N_TDtop5"

FDtop30_AC <- FDtop30_int %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())
names(FDtop30_AC)[2]="N_FDtop30"

FDtop17_AC <- FDtop17_int %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())
names(FDtop17_AC)[2]="N_FDtop17"

FDtop10_AC <- FDtop10_int %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())
names(FDtop10_AC)[2]="N_FDtop10"

FDtop5_AC <- FDtop5_int %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())
names(FDtop5_AC)[2]="N_FDtop5"

PDtop30_AC <- PDtop30_int %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())
names(PDtop30_AC)[2]="N_PDtop30"

PDtop17_AC <- PDtop17_int %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())
names(PDtop17_AC)[2]="N_PDtop17"

PDtop10_AC <- PDtop10_int %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())
names(PDtop10_AC)[2]="N_PDtop10"

PDtop5_AC <- PDtop5_int %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())
names(PDtop5_AC)[2]="N_PDtop5"

ACstats1 <- merge(cellsperAC, TDtop30_AC, by="ccaa.shortname.en", all.x=TRUE)
#Percentage of hotspot cells in relation to the total number of hotspot cells in Spain
ACstats1$TD30_per <- ACstats1$N_TDtop30*100/nrow(facets1)

ACstats1 <- merge(ACstats1, TDtop17_AC, by="ccaa.shortname.en", all.x=TRUE)
ACstats1$TD17_per <- ACstats1$N_TDtop17*100/nrow(facets1)

ACstats1 <- merge(ACstats1, TDtop10_AC, by="ccaa.shortname.en", all.x=TRUE)
ACstats1$TD10_per <- ACstats1$N_TDtop10*100/nrow(facets1)

ACstats1 <- merge(ACstats1, TDtop5_AC, by="ccaa.shortname.en", all.x=TRUE)
ACstats1$TD5_per <- ACstats1$N_TDtop5*100/nrow(facets1)

ACstats1 <- merge(ACstats1, FDtop30_AC, by="ccaa.shortname.en", all.x=TRUE)
#Percentage of hotspot cells in relation to the total number of hotspot cells in Spain
ACstats1$FD30_per <- ACstats1$N_FDtop30*100/nrow(facets1)

ACstats1 <- merge(ACstats1, FDtop17_AC, by="ccaa.shortname.en", all.x=TRUE)
ACstats1$FD17_per <- ACstats1$N_FDtop17*100/nrow(facets1)

ACstats1 <- merge(ACstats1, FDtop10_AC, by="ccaa.shortname.en", all.x=TRUE)
ACstats1$FD10_per <- ACstats1$N_FDtop10*100/nrow(facets1)

ACstats1 <- merge(ACstats1, FDtop5_AC, by="ccaa.shortname.en", all.x=TRUE)
ACstats1$FD5_per <- ACstats1$N_FDtop5*100/nrow(facets1)


ACstats1 <- merge(ACstats1, PDtop30_AC, by="ccaa.shortname.en", all.x=TRUE)
#Percentage of hotspot cells in relation to the total number of hotspot cells in Spain
ACstats1$PD30_per <- ACstats1$N_PDtop30*100/nrow(facets1)

ACstats1 <- merge(ACstats1, PDtop17_AC, by="ccaa.shortname.en", all.x=TRUE)
ACstats1$PD17_per <- ACstats1$N_PDtop17*100/nrow(facets1)

ACstats1 <- merge(ACstats1, PDtop10_AC, by="ccaa.shortname.en", all.x=TRUE)
ACstats1$PD10_per <- ACstats1$N_PDtop10*100/nrow(facets1)

ACstats1 <- merge(ACstats1, PDtop5_AC, by="ccaa.shortname.en", all.x=TRUE)
ACstats1$PD5_per <- ACstats1$N_PDtop5*100/nrow(facets1)

ACstats1$geometry.x <- NULL
ACstats1$geometry.y <- NULL
ACstats1$geometry.x <- NULL
ACstats1$geometry.y <- NULL
ACstats1$geometry.x <- NULL
ACstats1$geometry.y <- NULL
ACstats1$geometry.x <- NULL
ACstats1$geometry.y <- NULL
ACstats1$geometry.x <- NULL
ACstats1$geometry.y <- NULL
ACstats1$geometry.x <- NULL

write_xlsx(ACstats1, 'Data/TDFDPD_perAC_22Dec.xlsx') 


library(sf)
library(dplyr)

# Crear un identificador único para cada celda en facets1
facets1 <- facets1 %>%
  mutate(id = row_number())

# Crear las capas de celdas top30 para TD, FD y PD
TDtop30 <- facets1 %>% filter(top30TD == "top30")
FDtop30 <- facets1 %>% filter(top30FD == "top30")
PDtop30 <- facets1 %>% filter(top30PD == "top30")

# Calcular las intersecciones entre estas capas
interseccion_TD_FD <- st_intersection(TDtop30, FDtop30) %>% st_drop_geometry() %>% distinct(id)
interseccion_FD_PD <- st_intersection(FDtop30, PDtop30) %>% st_drop_geometry() %>% distinct(id)
interseccion_TD_PD <- st_intersection(TDtop30, PDtop30) %>% st_drop_geometry() %>% distinct(id)

# Calcular el número de celdas en cada intersección sin duplicados
num_celdas_TD <- nrow(TDtop30)
num_celdas_FD <- nrow(FDtop30)
num_celdas_PD <- nrow(PDtop30)

num_celdas_interseccion_TD_FD <- nrow(interseccion_TD_FD)
num_celdas_interseccion_FD_PD <- nrow(interseccion_FD_PD)
num_celdas_interseccion_TD_PD <- nrow(interseccion_TD_PD)

# Calcular los porcentajes de solapamiento correctamente
porcentaje_solapamiento_TD_FD <- (num_celdas_interseccion_TD_FD / num_celdas_TD) * 100
porcentaje_solapamiento_FD_PD <- (num_celdas_interseccion_FD_PD / num_celdas_FD) * 100
porcentaje_solapamiento_TD_PD <- (num_celdas_interseccion_TD_PD / num_celdas_PD) * 100

# Mostrar los resultados
porcentaje_solapamiento_TD_FD
porcentaje_solapamiento_FD_PD
porcentaje_solapamiento_TD_PD

library(dplyr)

# Contar las coincidencias entre top30TD y top30FD
coincidencias_TD_FD <- facets1 %>%
  filter(top30TD == "top30" & top30FD == "top30") %>%
  tally()

# Contar las coincidencias entre top30FD y top30PD
coincidencias_FD_PD <- facets1 %>%
  filter(top30FD == "top30" & top30PD == "top30") %>%
  tally()

# Contar las coincidencias entre top30TD y top30PD
coincidencias_TD_PD <- facets1 %>%
  filter(top30TD == "top30" & top30PD == "top30") %>%
  tally()

# Contar las coincidencias entre top30TD, top30FD y top30PD
coincidencias_TD_FD_PD <- facets1 %>%
  filter(top30TD == "top30" & top30FD == "top30" & top30PD == "top30") %>%
  tally()

# Mostrar los resultados
coincidencias_TD_FD$n
coincidencias_FD_PD$n
coincidencias_TD_PD$n
coincidencias_TD_FD_PD$n


# Instalar y cargar el paquete VennDiagram
if (!require("VennDiagram")) {
  install.packages("VennDiagram")
}
library(VennDiagram)

# Calcular las coincidencias entre las combinaciones
coincidencias_TD_FD <- nrow(facets1 %>% filter(top30TD == "top30" & top30FD == "top30"))
coincidencias_FD_PD <- nrow(facets1 %>% filter(top30FD == "top30" & top30PD == "top30"))
coincidencias_TD_PD <- nrow(facets1 %>% filter(top30TD == "top30" & top30PD == "top30"))
coincidencias_TD_FD_PD <- nrow(facets1 %>% filter(top30TD == "top30" & top30FD == "top30" & top30PD == "top30"))

# Calcular los porcentajes de solapamiento
porcentaje_TD_FD <- (coincidencias_TD_FD / 1157) * 100
porcentaje_FD_PD <- (coincidencias_FD_PD / 1151) * 100
porcentaje_TD_PD <- (coincidencias_TD_PD / 1157) * 100
porcentaje_TD_FD_PD <- (coincidencias_TD_FD_PD / 1157) * 100

# Crear etiquetas de porcentajes
label_TD_FD <- paste0(round(porcentaje_TD_FD, 2), "%")
label_FD_PD <- paste0(round(porcentaje_FD_PD, 2), "%")
label_TD_PD <- paste0(round(porcentaje_TD_PD, 2), "%")
label_TD_FD_PD <- paste0(round(porcentaje_TD_FD_PD, 2), "%")

# Crear un diagrama de Venn con porcentajes de solapamiento
venn.plot <- draw.triple.venn(
  area1 = 1157,
  area2 = 1151,
  area3 = 1151,
  n12 = coincidencias_TD_FD,
  n23 = coincidencias_FD_PD,
  n13 = coincidencias_TD_PD,
  n123 = coincidencias_TD_FD_PD,
  category = c("TD Top 30", "FD Top 30", "PD Top 30"),
  fill = c("red", "green", "blue"),
  alpha = 0.5,
  label.col = "white",
  cat.cex = 1.5,
  cat.col = c("red", "green", "blue"),
  lwd = 2,
  cex = 1.5,
  scaled = TRUE,
  euler.d = TRUE,
  cat.pos = c(-15, 15, 0),
  cat.dist = c(0.05, 0.05, 0.05)
)

# Agregar porcentajes de solapamiento al diagrama de Venn
grid.text(label_TD_FD, x = 0.4, y = 0.6, gp = gpar(fontsize = 12, col = "black"))
grid.text(label_FD_PD, x = 0.7, y = 0.6, gp = gpar(fontsize = 12, col = "black"))
grid.text(label_TD_PD, x = 0.55, y = 0.25, gp = gpar(fontsize = 12, col = "black"))
grid.text(label_TD_FD_PD, x = 0.55, y = 0.5, gp = gpar(fontsize = 12, col = "black"))

# Guardar el diagrama de Venn como un archivo PNG
png(filename = "venn_diagram.png")
grid.draw(venn.plot)
dev.off()

# Filtrar las celdas Top 30 de TD, FD y PD
TD_top30 <- facets1 %>% filter(top30TD == "top30")
FD_top30 <- facets1 %>% filter(top30FD == "top30")
PD_top30 <- facets1 %>% filter(top30PD == "top30")

# Calcular el área total de las celdas Top 30 de cada categoría
TD_top30_area <- sum(st_area(TD_top30))
FD_top30_area <- sum(st_area(FD_top30))
PD_top30_area <- sum(st_area(PD_top30))

# Calcular intersecciones
TD_FD_intersection <- st_intersection(TD_top30, FD_top30)
FD_PD_intersection <- st_intersection(FD_top30, PD_top30)
TD_PD_intersection <- st_intersection(TD_top30, PD_top30)

# Calcular el área de las intersecciones
TD_FD_area <- sum(st_area(TD_FD_intersection))
FD_PD_area <- sum(st_area(FD_PD_intersection))
TD_PD_area <- sum(st_area(TD_PD_intersection))

# Calcular porcentajes de solapamiento
TD_FD_overlap <- (TD_FD_area / TD_top30_area) * 100  # Porcentaje de TD solapado con FD
FD_PD_overlap <- (FD_PD_area / FD_top30_area) * 100  # Porcentaje de FD solapado con PD
TD_PD_overlap <- (TD_PD_area / TD_top30_area) * 100  # Porcentaje de TD solapado con PD

# Resultados
cat("Porcentaje de solapamiento de TD con FD (Top 30):", TD_FD_overlap, "%\n")
cat("Porcentaje de solapamiento de FD con PD (Top 30):", FD_PD_overlap, "%\n")
cat("Porcentaje de solapamiento de TD con PD (Top 30):", TD_PD_overlap, "%\n")

# Instalar y cargar el paquete VennDiagram
if (!require("VennDiagram")) {
  install.packages("VennDiagram")
}
library(VennDiagram)

# Definir las etiquetas de porcentajes
label_TD_FD <- "TD-FD: 49.87%"
label_TD_PD <- "TD-PD: 33.79%"
label_FD_PD <- "FD-PD: 53.95%"
label_TD_FD_PD <- "TD-FD-PD: 24.20%"

# Crear un diagrama de Venn con los porcentajes de solapamiento
venn.plot <- draw.triple.venn(
  area1 = 1157,
  area2 = 1151,
  area3 = 1151,
  n12 = 49.87,
  n23 = 53.95,
  n13 = 33.79,
  n123 = 24.20,
  category = c("TD Top 30", "FD Top 30", "PD Top 30"),
  fill = c("red", "green", "blue"),
  alpha = 0.5,
  label.col = "white",
  cat.cex = 1.5,
  cat.col = c("red", "green", "blue"),
  lwd = 2,
  cex = 1.5,
  scaled = FALSE,  # Cambiar a FALSE para mostrar los porcentajes como texto
  euler.d = TRUE,
  cat.pos = c(-15, 15, 0),
  cat.dist = c(0.05, 0.05, 0.05)
)

# Agregar porcentajes de solapamiento al diagrama de Venn como texto
grid.text(label_TD_FD, x = 0.3, y = 0.65, gp = gpar(fontsize = 12, col = "black"))
grid.text(label_FD_PD, x = 0.7, y = 0.65, gp = gpar(fontsize = 12, col = "black"))
grid.text(label_TD_PD, x = 0.5, y = 0.25, gp = gpar(fontsize = 12, col = "black"))
grid.text(label_TD_FD_PD, x = 0.5, y = 0.5, gp = gpar(fontsize = 12, col = "black"))

# Guardar el diagrama de Venn como un archivo PNG
png(filename = "Figures/venn_diagram.png")
grid.draw(venn.plot)
dev.off()

library(dplyr)

# Crear un subconjunto eliminando filas donde top30TD, top30FD y top30PD son todos NA, esto es para calcular porcentajes para diagramas venn
facets1_subset <- facets1 %>%
  filter(!(is.na(top30TD) & is.na(top30FD) & is.na(top30PD)))

# Verificar el resultado
summary(facets1_subset)

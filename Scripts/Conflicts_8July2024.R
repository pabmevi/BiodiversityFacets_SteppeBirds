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

tmap_save(mapaSES, filename = "Figures/3facetsSEShotsp.png")

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

Hotspotcells <- subset(PVint2, overlapSES == "Hotspot") #There are 240 hotspot cells in Spain

# Spatial intersect between ACs and hotspot cells
hotsp_AC_int <- st_intersection(Hotspotcells, AC)

#Calculating the area of each intersection section
hotsp_AC_int$areainthots_AC <- st_area(hotsp_AC_int)

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

names(hotspots_perAC)[2]="N hotspot cells"

#Percentage of hotspot cells in relation to the number of cells per AC
AC_stats <- merge(cells_perAC, hotspots_perAC, by="ccaa.shortname.en")

AC_stats$Perc_hotp_AC <- AC_stats$`N hotspot cells`*100/AC_stats$`N cells`

#Percentage of hotspot cells in relation to the total number of hotspot cells in Spain
AC_stats$Perc_hotp_Spain <- AC_stats$`N hotspot cells`*100/240

#As one hotspot cell may fall in different ACs, calculating areas could be more precise. 

# Calculating the total area from hotspots cells
areaHots <- st_area(Hotspotcells)
sum(areaHots) #The 240 hotspot cells equal to 23844622279 [m^2] or 23844.62 km2

# Calculating the area from each AC 
AC$areaAC <- st_area(AC)

#Summing the total hotspot area per AC 
areas_AC_hotsp <- hotsp_AC_int %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(Hotsarea_pAC = sum(areainthots_AC)) 

areas_AC_hotsp$geometry <- NULL

areas_AC_hotsp1 <- merge(areas_AC_hotsp, AC[,c("ccaa.shortname.en", "areaAC")], by="ccaa.shortname.en" ) 

areas_AC_hotsp1 <- drop_units(areas_AC_hotsp1)

AC_stats <- merge(AC_stats, areas_AC_hotsp1, by="ccaa.shortname.en" ) 
AC_stats <- drop_units(AC_stats)
AC_stats$geometry.x <- NULL

AC_stats$per_hotareaAC <- AC_stats$Hotsarea_pAC*100/AC_stats$areaAC
AC_stats$per_hotareaSpain <- AC_stats$Hotsarea_pAC*100/23844622279 

# I want to know the number of conflict cells inside each AC, and their corresponding percentage in relation to 
# the number of cells per AC and in relation to the total conflict cells in Spain. Also the corresponding areas. 

#Selecting only strong conflict cells

Conflictcells9 <- subset(PVint2, extendedconflicts == "Conflict 9") #There are 55 conflict 9 cells in Spain
Conflictcells8 <- subset(PVint2, extendedconflicts == "Conflict 8") #There are 46 conflict 8 cells in Spain
Conflictcells7 <- subset(PVint2, extendedconflicts == "Conflict 7") #There are 41 conflict 7 cells in Spain

# Spatial intersect between ACs/AC and hotspot cells
conflict_AC_int9 <- st_intersection(Conflictcells9, AC)
conflict_AC_int8 <- st_intersection(Conflictcells8, AC)
conflict_AC_int7 <- st_intersection(Conflictcells7, AC)

#Calculating the area of each intersection section
conflict_AC_int9$areaintconf_AC9 <- st_area(conflict_AC_int9)
conflict_AC_int8$areaintconf_AC8 <- st_area(conflict_AC_int8)
conflict_AC_int7$areaintconf_AC7 <- st_area(conflict_AC_int7)

conflict_perAC9 <- conflict_AC_int9 %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

conflict_perAC8 <- conflict_AC_int8 %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

conflict_perAC7 <- conflict_AC_int7 %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(UTMCODE = n())

names(conflict_perAC9)[2]="N conflict cells 9"
names(conflict_perAC8)[2]="N conflict cells 8"
names(conflict_perAC7)[2]="N conflict cells 7"

#Percentage of strong conflict cells (3, 2 or 1 level) in relation to the number of cells per AC
AC_confstats9 <- merge(AC_stats, conflict_perAC9, by="ccaa.shortname.en", all.x=TRUE)
AC_confstats8 <- merge(AC_stats, conflict_perAC8, by="ccaa.shortname.en", all.x=TRUE)
AC_confstats7 <- merge(AC_stats, conflict_perAC7, by="ccaa.shortname.en", all.x=TRUE)

#Percentage of strong conflict cells (3, 2 or 1 level) in relation to the total number of strong conflict cells (3, 2 or 1 level) in Spain
AC_confstats9$Perc_conf_Spain9 <- AC_confstats9$`N conflict cells 9`*100/55
AC_confstats8$Perc_conf_Spain8 <- AC_confstats8$`N conflict cells 8`*100/46
AC_confstats7$Perc_conf_Spain7 <- AC_confstats7$`N conflict cells 7`*100/41

AC_confstats9$geometry.x <- NULL
AC_confstats9$geometry.y <- NULL

AC_confstats8$geometry.x <- NULL
AC_confstats8$geometry.y <- NULL

AC_confstats7$geometry.x <- NULL
AC_confstats7$geometry.y <- NULL

#As one conflict cell may fall in different ACs, calculating areas could be more precise. 

# Calculating the total area from hotspots cells
areaconflict9 <- st_area(Conflictcells9)
areaconflict8 <- st_area(Conflictcells8)
areaconflict7 <- st_area(Conflictcells7)

areaconflict9sum  <- sum(areaconflict9) 
areaconflict8sum  <- sum(areaconflict8) 
areaconflict7sum  <- sum(areaconflict7) 

totalconfarea  <- (areaconflict9sum +areaconflict8sum +areaconflict7sum) 
#The sum of the 142 strong conflict cells (9, 8 and 7) equal to 16454744256 [m^2] or 16454.74 km2

#Summing the total hotspot area per AC 
areas_AC_confl9 <- conflict_AC_int9 %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(Confarea_pAC9 = sum(areaintconf_AC9)) 

areas_AC_confl8 <- conflict_AC_int8 %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(Confarea_pAC8 = sum(areaintconf_AC8))

areas_AC_confl7 <- conflict_AC_int7 %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(Confarea_pAC7 = sum(areaintconf_AC7))

areas_AC_confl9$geometry <- NULL
areas_AC_confl8$geometry <- NULL
areas_AC_confl7$geometry <- NULL

areas_AC_confl9_1 <- merge(areas_AC_confl9, AC[,c("ccaa.shortname.en", "areaAC")], by="ccaa.shortname.en" ) 
areas_AC_confl9_1 <- drop_units(areas_AC_confl9_1)

AC_confstats9 <- merge(AC_confstats9, areas_AC_confl9_1[,c("ccaa.shortname.en", "Confarea_pAC9")], by="ccaa.shortname.en", all.x=TRUE ) 
AC_confstats9 <- drop_units(AC_confstats9)

AC_confstats9$per_confareaAC9 <- AC_confstats9$Confarea_pAC*100/AC_confstats9$areaAC
AC_confstats9$per_confareaSpain9 <- AC_confstats9$Confarea_pAC*100/5484914752 

areas_AC_confl8_1 <- merge(areas_AC_confl8, AC[,c("ccaa.shortname.en", "areaAC")], by="ccaa.shortname.en" ) 
areas_AC_confl8_1 <- drop_units(areas_AC_confl8_1)

AC_confstats8 <- merge(AC_confstats8, areas_AC_confl8_1[,c("ccaa.shortname.en", "Confarea_pAC8")], by="ccaa.shortname.en", all.x=TRUE ) 
AC_confstats8 <- drop_units(AC_confstats8)

AC_confstats8$per_confareaAC8 <- AC_confstats8$Confarea_pAC*100/AC_confstats8$areaAC
AC_confstats8$per_confareaSpain8 <- AC_confstats8$Confarea_pAC*100/4600617564 

areas_AC_confl7_1 <- merge(areas_AC_confl7, AC[,c("ccaa.shortname.en", "areaAC")], by="ccaa.shortname.en" ) 
areas_AC_confl7_1 <- drop_units(areas_AC_confl7_1)

AC_confstats7 <- merge(AC_confstats7, areas_AC_confl7_1[,c("ccaa.shortname.en", "Confarea_pAC7")], by="ccaa.shortname.en", all.x=TRUE ) 
AC_confstats7 <- drop_units(AC_confstats7)

AC_confstats7$per_confareaAC7 <- AC_confstats7$Confarea_pAC*100/AC_confstats7$areaAC
AC_confstats7$per_confareaSpain7 <- AC_confstats7$Confarea_pAC*100/4100189172 

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

#As one no-go cell may fall in different ACs, calculating areas could be more precise. 

# Calculating the total area from hotspots cells
areanogo_vh <- st_area(Veryhigh_nogocells)
sum(areanogo_vh) #The 98 no-go cells equal to 9658900791 [m^2] 

areanogo_h <- st_area(High_nogocells)
sum(areanogo_h) #The 374 no-go cells equal to 36987363492 [m^2] 

areanogo_m <- st_area(moderate_nogocells)
sum(areanogo_m) #The 566 no-go cells equal to 55753073269 [m^2]

#Summing the total nogo area per AC 
areas_AC_nogo_vh <- nogo_AC_int_veryhigh %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(nogoarea_pAC_vh = sum(areaintnogo_AC_vh)) 

areas_AC_nogo_vh$geometry <- NULL

#Summing the total nogo area per AC 
areas_AC_nogo_h <- nogo_AC_int_high %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(nogoarea_pAC_h = sum(areaintnogo_AC_h)) 

areas_AC_nogo_h$geometry <- NULL

#Summing the total nogo area per AC 
areas_AC_nogo_m <- nogo_AC_int_moderate %>% 
  group_by(ccaa.shortname.en) %>% 
  summarize(nogoarea_pAC_m = sum(areaintnogo_AC_m)) 

areas_AC_nogo_m$geometry <- NULL

areas_AC_nogo_vh <- merge(areas_AC_nogo_vh, AC[,c("ccaa.shortname.en", "areaAC")], by="ccaa.shortname.en") 
areas_AC_nogo_vh <- drop_units(areas_AC_nogo_vh)

areas_AC_nogo_h <- merge(areas_AC_nogo_h, AC[,c("ccaa.shortname.en", "areaAC")], by="ccaa.shortname.en") 
areas_AC_nogo_h <- drop_units(areas_AC_nogo_h)

areas_AC_nogo_m <- merge(areas_AC_nogo_m, AC[,c("ccaa.shortname.en", "areaAC")], by="ccaa.shortname.en") 
areas_AC_nogo_m <- drop_units(areas_AC_nogo_m)

AC_stats_total3 <- merge(ACnogo_stats3, areas_AC_nogo_vh[,c("ccaa.shortname.en", "nogoarea_pAC_vh")], by="ccaa.shortname.en", all.x=TRUE ) 
AC_stats_total3 <- drop_units(AC_stats_total3)

AC_stats_total3$per_nogoareaAC_vh <- AC_stats_total3$nogoarea_pAC*100/AC_stats_total3$areaAC
AC_stats_total3$per_nogoareaSpain_vh <- AC_stats_total3$nogoarea_pAC*100/9658900791 

AC_stats_total2 <- merge(ACnogo_stats2, areas_AC_nogo_h[,c("ccaa.shortname.en", "nogoarea_pAC_h")], by="ccaa.shortname.en", all.x=TRUE ) 
AC_stats_total2 <- drop_units(AC_stats_total2)

AC_stats_total2$per_nogoareaAC_h <- AC_stats_total2$nogoarea_pAC*100/AC_stats_total2$areaAC
AC_stats_total2$per_nogoareaSpain_h <- AC_stats_total2$nogoarea_pAC*100/36987363492 

AC_stats_total1 <- merge(ACnogo_stats1, areas_AC_nogo_m[,c("ccaa.shortname.en", "nogoarea_pAC_m")], by="ccaa.shortname.en", all.x=TRUE ) 
AC_stats_total1 <- drop_units(AC_stats_total1)

AC_stats_total1$per_nogoareaAC_m <- AC_stats_total1$nogoarea_pAC*100/AC_stats_total1$areaAC
AC_stats_total1$per_nogoareaSpain_m <- AC_stats_total1$nogoarea_pAC*100/55753073269 

# Merging AC_stats_total3, AC_stats_total2, and AC_stats_total1
AC_stats_merged <- AC_stats_total1 %>%
  full_join(AC_stats_total2, by = "ccaa.shortname.en")

AC_stats_total <- AC_stats_merged %>%
  full_join(AC_stats_total3, by = "ccaa.shortname.en")

# Merging two columns in one to show N hotspot cells and percentage in relation to the number of cells in each AC
AC_stats_total$`N_hotspot cells` <- paste(AC_stats_total$`N hotspot cells.x`, "[",round(AC_stats_total$Perc_hotp_AC.x,1),"%]", sep ="")

# Calculating the % highest conflict cells (7, 8 and 9 conflict levels) in relation to the number of cells in each AC
AC_stats_total$`Sum highest conflict cells` <- rowSums(AC_stats_total[, c("N conflict cells 7", "N conflict cells 8", "N conflict cells 9")], na.rm = TRUE)
AC_stats_total$Perc_highestconf <- (AC_stats_total$`Sum highest conflict cells` / AC_stats_total$`N cells.x`) * 100
AC_stats_total$`N highest conflict cells` <- paste(AC_stats_total$`Sum highest conflict cells`, "[",round(AC_stats_total$Perc_highestconf,1),"%]", sep ="")
AC_stats_total$Highest_conf_Spain <- (AC_stats_total$`Sum highest conflict cells` / 142) * 100

AC_stats_total$m_nogocells <- paste(AC_stats_total$`moderate no-go cells`, "[",round(AC_stats_total$Perc_nogo_AC_m,1),"%]", sep ="")
AC_stats_total$h_nogocells <- paste(AC_stats_total$`high no-go cells`, "[",round(AC_stats_total$Perc_nogo_AC_h,1),"%]", sep ="")
AC_stats_total$vh_nogocells <- paste(AC_stats_total$`very high no-go cells`, "[",round(AC_stats_total$Perc_nogo_AC_vh,1),"%]", sep ="")

AC_stats_total$geometry.x.x <- NULL
AC_stats_total$geometry.y.x <- NULL
AC_stats_total$geometry.x.y <- NULL
AC_stats_total$geometry.y.y <- NULL
AC_stats_total$geometry.x <- NULL
AC_stats_total$geometry.y <- NULL
AC_stats_total$`N hotspot cells.y` <- NULL
AC_stats_total$Perc_hotp_AC.y <- NULL
AC_stats_total$Perc_hotp_Spain.y <- NULL
AC_stats_total$Hotsarea_pAC.y <- NULL
AC_stats_total$areaAC.y <- NULL
AC_stats_total$per_hotareaAC.y <- NULL
AC_stats_total$per_hotareaSpain.y <- NULL
AC_stats_total$Hotsarea_pAC.y <- NULL
AC_stats_total$Hotsarea_pAC.y <- NULL
AC_stats_total$Hotsarea_pAC.y <- NULL
AC_stats_total$Hotsarea_pAC.y <- NULL
AC_stats_total$`N cells` <- NULL
AC_stats_total$`N hotspot cells` <- NULL
AC_stats_total$Perc_hotp_AC <- NULL
AC_stats_total$Perc_hotp_Spain <- NULL
AC_stats_total$Hotsarea_pAC <- NULL
AC_stats_total$areaAC <- NULL
AC_stats_total$per_hotareaAC <- NULL
AC_stats_total$per_hotareaSpain <- NULL
AC_stats_total$`N hotspot cells.x`<- NULL
AC_stats_total$Perc_hotp_AC.x<- NULL
AC_stats_total$`Sum highest conflict cells`<- NULL
AC_stats_total$Perc_highestconf<- NULL
AC_stats_total$`moderate no-go cells`<- NULL
AC_stats_total$Perc_nogo_AC_m<- NULL
AC_stats_total$`high no-go cells`<- NULL
AC_stats_total$Perc_nogo_AC_h<- NULL
AC_stats_total$`very high no-go cells`<- NULL
AC_stats_total$Perc_nogo_AC_vh<- NULL
write_xlsx(AC_stats_total, 'Data/AC_stats_5November.xlsx') 

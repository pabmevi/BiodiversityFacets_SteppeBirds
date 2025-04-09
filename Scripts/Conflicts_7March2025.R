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
library(biscale)
library(terra)
library(cowplot)
library(ggvenn)

setwd("~/GitHub/BiodiversityFacets")

####################################################################
############### Overlapping biodiversity facets############################

# clean environment
rm(list = ls())

# Adding Autonomous communities from Spain
AC <-esp_get_ccaa()
AC <- AC[!AC$iso2.ccaa.name.es %in% c("Canarias"),]
AC <- st_transform(AC, 25830) 

# load Spanish UTM grid
malla <- st_read("Spatial_Data/Malla_municipios/Malla10x10_Ter_p.shp")
malla <- st_transform(malla, 25830) 

###########Percentage of cells occupied by PV #################################
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

length(unique(PVint$UTMCODE))

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

PVint2 <- PVint11[, c("UTMCODE","species", "SESFRic", "SESPD", "PercPV_cell")]

#Next, will calculate percentile 30, to overlap biodiversity facets 
# Calculating percentil 70 for each index
cutoff_sp_richn <- quantile(PVint2$species, probs = 0.70, na.rm = TRUE)
cutoff_SESfric <- quantile(PVint2$SESFRic, probs = 0.70, na.rm = TRUE)
cutoff_SESPD <- quantile(PVint2$SESPD, probs = 0.70, na.rm = TRUE)

# New variable for each index. It is TRUE if the value is inside the top 30%
PVint2$top_sp_richn <- PVint2$species > cutoff_sp_richn
PVint2$top_SESfric <- PVint2$SESFRic > cutoff_SESfric
PVint2$top_SESPD <- PVint2$SESPD > cutoff_SESPD

# Creating a new variable to represent overlapping using SES values
PVint2$overlapSES <- with(PVint2, ifelse(top_sp_richn & top_SESfric & top_SESPD, "Three-facet diverse",
                                          ifelse(top_sp_richn & top_SESfric, "Top 30% TD and FD",
                                                 ifelse(top_sp_richn & top_SESPD, "Top 30% TD and PD",
                                                        ifelse(top_SESfric & top_SESPD, "Top 30% FD and PD",
                                                               ifelse(top_sp_richn, "Top 30% TD",
                                                                      ifelse(top_SESfric, "Top 30% FD",
                                                                             ifelse(top_SESPD, "Top 30% PD", "Outside top 30%"))))))))
# Sorting overlapSES levels
PVint2$overlapSES <- factor(PVint2$overlap, levels = c("Three-facet diverse", "Top 30% TD and FD", 
                                                    "Top 30% TD and PD", "Top 30% FD and PD", 
                                                    "Top 30% TD", "Top 30% FD", "Top 30% PD", 
                                                    "Outside top 30%"))

# Colours for each overlapping category 
colores <- c("Three-facet diverse" = "darkred", "Top 30% TD and FD" = "#332288", "Top 30% TD and PD" = "#DDCC77", 
             "Top 30% FD and PD" = "#CC79A7", "Top 30% TD" = "grey", "Top 30% FD" = "goldenrod3",
             "Top 30% PD" = "slategray3", "Outside top 30%" = "white")

PVint2$overlapSES[is.na(PVint2$overlapSES)] <- "Outside top 30%"

# Final map with TD, SESFRIC and SESPD
mapaSES <- tm_shape(PVint2) +
  tm_fill(col = "overlapSES", palette = colores, title = "Index overlapping", style = "cat")

# Adding Spain borders to final map
mapaSES <- tm_shape(PVint2) +
  tm_fill(col = "overlapSES", palette = colores, title = "Index overlapping", style = "cat") +
  tm_shape(AC) +
  tm_borders() +
  tm_scale_bar(breaks = c(0, 50, 100, 150, 200), position = c("left", "bottom")) +  
  tm_compass(type = "arrow", position = c("right", "top"))

tmap_save(mapaSES, filename = "Figures/3facetsSES_22decemb.png")
tmap_save(mapaSES, filename = "Figures/3facetsSES_22decemb.pdf")

# Generating a Venn diagramm overlapping three biodiversity facets
# First countinng the quantity of overlapping cells 
summary(PVint2$overlapSES) #Total cells 3838, however, 1819 are outside any top 30%, then cells in any top 30% count 2019
#This means that to generate a venn diagramm with overlapping facets, I need to consider 2019 cells only.  
# Top 30% TD, FD and PD count 240, this is 11.89% of 2019 cells, 
# Top 30% TD and FD count 276, this is 13.67%
# Top 30% TD and PD count 85, this is 4.21%
# Top 30% FD and PD count 381, this is 18.87%
# Top 30% TD count 337, this is 16.69%
# Top 30% FD count 255, this is 12.63%
# Top 30% PD count 445, this is 22.04%

# Selecting only top cells for each facet
venn_data <- list(
  "TD" = PVint2$UTMCODE[PVint2$top_sp_richn == TRUE],
  "FD" = PVint2$UTMCODE[PVint2$top_SESfric == TRUE],
  "PD" = PVint2$UTMCODE[PVint2$top_SESPD == TRUE]
)

# Generating the venn diagram
venn_plot <-ggvenn(venn_data, 
       fill_color = c("grey", "goldenrod3", "slategray3"), 
       stroke_size = 0.5, 
       text_size = 4)

ggsave("Figures/venn_diagram.png", plot = venn_plot, width = 8, height = 6)

# I want to know the number of 3, 2 and 1 facet diverse cells in the study area. Just for curiosity intersecting with AC to know the number of diverse cells per AC. 

Hotspotcells <- subset(PVint2, overlapSES == "Three-facet diverse") #240 hotspot cells in Spain
Twofacetdiv <- subset(PVint2, overlapSES %in% c("Top 30% TD and FD", "Top 30% TD and PD", "Top 30% FD and PD"))#742 Twofacetdiv cells in Spain
Onefacetdiv <- subset(PVint2, overlapSES %in% c("Top 30% TD", "Top 30% FD", "Top 30% PD")) #1037 Onefacetdiv cells in Spain
Outsidetop30 <- subset(PVint2, overlapSES %in% c("Outside top 30%")) #1819 outside top 30% cells in Spain

# Counting the number of 3, 2 and 1 facet diverse cells where PV have been installed. 
hotspotswPV <- PVint2[!is.na(PVint2$PercPV_cell) & PVint2$overlap == "Three-facet diverse", ] #142 out of 142 Three-facet diverse cells contain PV
TwofacetwPV <- PVint2[!is.na(PVint2$PercPV_cell) & (PVint2$overlap == "Top 30% FD and PD"| PVint2$overlap == "Top 30% TD and PD"| PVint2$overlap == "Top 30% TD and FD"),]
##368 out of 742 Three-facet diverse cells contain PV
OnefacetwPV <- PVint2[!is.na(PVint2$PercPV_cell) & (PVint2$overlap == "Top 30% TD"| PVint2$overlap == "Top 30% FD"| PVint2$overlap == "Top 30% PD"),]
##471 out of 1037 Three-facet diverse cells contain PV

# Spatial intersect between ACs and hotspot cells
hotsp_AC_int <- st_intersection(Hotspotcells, AC)
twofacet_AC_int <- st_intersection(Twofacetdiv, AC)
onefacet_AC_int <- st_intersection(Onefacetdiv, AC)

#Calculating the area of each intersection section
hotsp_AC_int$areainthots_AC <- st_area(hotsp_AC_int)
twofacet_AC_int$areainttwof_AC <- st_area(twofacet_AC_int)
onefacet_AC_int$areaintonef_AC <- st_area(onefacet_AC_int)

#Now I want to scale TD, FD, and PD in a range from 0-1 TO USE ZONATION SOFTWARE.
basic_function <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}
PVint2_1<-PVint2
PVint2_1$SESFRic <- basic_function(PVint2_1$SESFRic)
PVint2_1$SESPD <- basic_function(PVint2_1$SESPD)
PVint2_1$species <- basic_function(PVint2_1$species)
# Replacing NAs by 0 in PV occupancy cells
PVint2_1$PercPV_cell[is.na(PVint2_1$PercPV_cell)] <- 0
#Also scaling PV occupancy from 0-1. However, maximum occupancy % per cell is around 12%
PVint2_1$PercPV_cell0  <- PVint2_1$PercPV_cell #Keeping the PV percentage occupancy per cell before scaling 0-1
PVint2_1$PercPV_cell <- basic_function(PVint2_1$PercPV_cell)

# PV occupancy in tertiles
terciles <- quantile(PVint2$PercPV_cell, probs = c(0, 0.33, 0.67, 1), na.rm = TRUE)

# Creating the new column to show high, medium and low occupancy of PV plants inside each cell based on tertiles.
PVint2_1 <- PVint2_1 %>%
  mutate(PV_occupancy = cut(PercPV_cell, breaks = terciles, labels = c("low", "medium", "high"), include.lowest = TRUE))

write_sf(PVint2_1, "Spatial_Data/3facets_scaled_1.shp")

facets <- st_read("Spatial_Data/3facets_scaled_1.shp")

# Now categorizing each biodiversity facet into 5, 10, 17 and 30% top values to facilitate the comparison of spatial patterns of facets
# Function to classify and map values 
plot_diversity <- function(data, column, title, borders, show_legend = TRUE, show_annotations = TRUE) {
  # Sorting
  data <- data %>%
    mutate(
      Rank = rank(-!!sym(column)),  
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
    mutate(Category = factor(Category, levels = c("Top 5%", "Top 10%", "Top 17%", "Top 30%")))  # Sorting
  
  p <- ggplot() +
    geom_sf(data = data, aes(fill = Category), color = NA, size = 0) +  # Removing white lies
    geom_sf(data = borders, fill = NA, color = "black", size = 0.3) +  # AC borders
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
      legend.position = if (show_legend) c(0.8, 0.1) else "none",  
      legend.text = element_text(size = 8),  
      legend.title = element_text(size = 8), 
      plot.title = element_text(size = 12, hjust = 0.5), 
      plot.margin = margin(5, 5, 5, 5), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      axis.text = element_blank(), 
      axis.ticks = element_blank() 
    )
  
  if (show_annotations) {
    p <- p + annotation_scale(location = "bl", width_hint = 0.3)
  }
  
  return(p)
}

# Maps for each facet
map_td <- plot_diversity(facets, "species", "Taxonomic Diversity (TD)", borders = AC, show_legend = FALSE, show_annotations = FALSE)
map_fd <- plot_diversity(facets, "SESFRic", "Functional Diversity (FD)", borders = AC, show_legend = FALSE, show_annotations = FALSE)
map_pd <- plot_diversity(facets, "SESPD", "Phylogenetic Diversity (PD)", borders = AC, show_legend = TRUE, show_annotations = TRUE)

# Combining maps in a single graph 
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

ggsave("Figures/combined_map1_1.png", combined_map1, width = 20, height = 10, units = "cm", dpi = 300)
ggsave("Figures/combined_map1_1.pdf", combined_map1, width = 20, height = 10, units = "cm", dpi = 300)
ggsave("Figures/combined_map2_1.png", combined_map2, width = 20, height = 10, units = "cm", dpi = 300)
ggsave("Figures/combined_map2_1.pdf", combined_map2, width = 20, height = 10, units = "cm", dpi = 300)

#################################################################
#####I rasterized TD, FD and PD in ArcMap, as it is slighltly more precise than R. Then each rasterized facet was used in Zonation 
#####to generate prioritized maps. For Zonation process I used the whole values not subsets of top values for each facet.
################################################################

#Using maps generated by Zonation 5 to generate Top 30% maps of each prioritization scenario (TD, FD, PD, TD-FD, FD-PD, TD-FD-PD)
#First using the map tht combines TD, FD, and PD
Hotsp_Zonation <- raster("Spatial_Data/Zonation/rankmap_hotsp.tif")

Hotsp_percentil30 <- quantile(Hotsp_Zonation, 0.7)
Hotsp_percentil17 <- quantile(Hotsp_Zonation, 0.83)
Hotsp_percentil10 <- quantile(Hotsp_Zonation, 0.9)
Hotsp_percentil5 <- quantile(Hotsp_Zonation, 0.95)

Hotsp_classes <- calc(Hotsp_Zonation, fun = function(x) {
  ifelse(x >= Hotsp_percentil5, 4,  # Top 5%
         ifelse(x >= Hotsp_percentil10, 3,  # Top 10%
                ifelse(x >= Hotsp_percentil17, 2,  # Top 17%
                       ifelse(x >= Hotsp_percentil30, 1, 0))))  # Top 30% y fondo
})

# Definir nombres de clases
custom_colors <- c("ivory2", "yellow", "pink", "orange", "red") 
categories <- c("Below 30%", "Top 30%", "Top 17%", "Top 10%", "Top 5%")

# Convertir el raster a data frame para ggplot
df <- as.data.frame(rasterToPoints(Hotsp_classes))
colnames(df) <- c("x", "y", "category")

#Loading ACs
comm <-esp_get_ccaa()
comm <- comm[!comm$iso2.ccaa.name.es %in% c("Canarias"),]
comm <- st_transform(comm, st_crs(Hotsp_Zonation))


# Graficar con ggplot2
TDFDPD <- ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = factor(category))) +
  geom_sf(data = comm, fill = NA, color = "gray", linewidth = 0.5) +  # Agregar límites de comunidades en gris claro
  scale_fill_manual(values = custom_colors, labels = categories, name = "Hotspot Level") +
  theme_minimal() +
  labs(title = "Hotspot Prioritization Levels",
       x = "Longitude",
       y = "Latitude") +
  theme(legend.position = "right")


ggsave("Figures/TDFDPD_tops.tiff", plot = TDFDPD, width = 10, height = 8, dpi = 300, compression = "lzw")

dev.off()

#Now the map that contains TD
TD_Zonation <- raster("Spatial_Data/Zonation/rankmap_TD.tif")

TD_percentil30 <- quantile(TD_Zonation, 0.7)
TD_percentil17 <- quantile(TD_Zonation, 0.83)
TD_percentil10 <- quantile(TD_Zonation, 0.9)
TD_percentil5 <- quantile(TD_Zonation, 0.95)

TD_classes <- calc(TD_Zonation, fun = function(x) {
  ifelse(x >= TD_percentil5, 4,  # Top 5%
         ifelse(x >= TD_percentil10, 3,  # Top 10%
                ifelse(x >= TD_percentil17, 2,  # Top 17%
                       ifelse(x >= TD_percentil30, 1, 0))))  # Top 30% y fondo
})

# Definir nombres de clases
custom_colors <- c("gray91", "yellow", "pink", "orange", "red") 
categories <- c("Below 30%", "Top 30%", "Top 17%", "Top 10%", "Top 5%")

# Convertir el raster a data frame para ggplot
df1 <- as.data.frame(rasterToPoints(TD_classes))
colnames(df1) <- c("x", "y", "category")

# Graficar con ggplot2
TD <- ggplot() +
  geom_raster(data = df1, aes(x = x, y = y, fill = factor(category))) +
  geom_sf(data = comm, fill = NA, color = "gray", linewidth = 0.5) +
  scale_fill_manual(values = custom_colors, labels = categories, name = "") +
  theme_minimal() +
  labs(title = "Taxonomic diversity",
       x = "Longitude",
       y = "Latitude") +
  theme(legend.position = "left")

ggsave("Figures/TD_tops.tiff", plot = TD, width = 10, height = 8, dpi = 300, compression = "lzw")

#Now the map that contains FD
FD_Zonation <- raster("Spatial_Data/Zonation/rankmap_FD.tif")

FD_percentil30 <- quantile(FD_Zonation, 0.7)
FD_percentil17 <- quantile(FD_Zonation, 0.83)
FD_percentil10 <- quantile(FD_Zonation, 0.9)
FD_percentil5 <- quantile(FD_Zonation, 0.95)

FD_classes <- calc(FD_Zonation, fun = function(x) {
  ifelse(x >= FD_percentil5, 4,  # Top 5%
         ifelse(x >= FD_percentil10, 3,  # Top 10%
                ifelse(x >= FD_percentil17, 2,  # Top 17%
                       ifelse(x >= FD_percentil30, 1, 0))))  # Top 30% y fondo
})

# Definir nombres de clases
custom_colors <- c("gray91", "yellow", "pink", "orange", "red") 
categories <- c("Below 30%", "Top 30%", "Top 17%", "Top 10%", "Top 5%")

# Convertir el raster a data frame para ggplot
df2 <- as.data.frame(rasterToPoints(FD_classes))
colnames(df2) <- c("x", "y", "category")

# Graficar con ggplot2
FD <- ggplot() +
  geom_raster(data = df2, aes(x = x, y = y, fill = factor(category))) +
  geom_sf(data = comm, fill = NA, color = "gray", linewidth = 0.5) +
  scale_fill_manual(values = custom_colors, labels = categories, name = "") +
  theme_minimal() +
  labs(title = "Functional diversity",
       x = "Longitude",
       y = "Latitude") +
  theme(legend.position = "left")
ggsave("Figures/FD_tops.tiff", plot = FD, width = 10, height = 8, dpi = 300, compression = "lzw")

#Now the map that contains PD
PD_Zonation <- raster("Spatial_Data/Zonation/rankmap_PD.tif")

PD_percentil30 <- quantile(PD_Zonation, 0.7)
PD_percentil17 <- quantile(PD_Zonation, 0.83)
PD_percentil10 <- quantile(PD_Zonation, 0.9)
PD_percentil5 <- quantile(PD_Zonation, 0.95)

PD_classes <- calc(PD_Zonation, fun = function(x) {
  ifelse(x >= PD_percentil5, 4,  # Top 5%
         ifelse(x >= PD_percentil10, 3,  # Top 10%
                ifelse(x >= PD_percentil17, 2,  # Top 17%
                       ifelse(x >= PD_percentil30, 1, 0))))  # Top 30% y fondo
})

# Definir nombres de clases
custom_colors <- c("gray91", "yellow", "pink", "orange", "red") 
categories <- c("Below 30%", "Top 30%", "Top 17%", "Top 10%", "Top 5%")

# Convertir el raster a data frame para ggplot
df3 <- as.data.frame(rasterToPoints(PD_classes))
colnames(df3) <- c("x", "y", "category")

# Graficar con ggplot2
PD <- ggplot() +
  geom_raster(data = df3, aes(x = x, y = y, fill = factor(category))) +
  geom_sf(data = comm, fill = NA, color = "gray", linewidth = 0.5) +
  scale_fill_manual(values = custom_colors, labels = categories, name = "") +
  theme_minimal() +
  labs(title = "Phylogenetic diversity",
       x = "Longitude",
       y = "Latitude") +
  theme(legend.position = "left")
ggsave("Figures/PD_tops.tiff", plot = PD, width = 10, height = 8, dpi = 300, compression = "lzw")


#Now the map tht combines TD, PD
TDPD_Zonation <- raster("Spatial_Data/Zonation/rankmap_TDPD.tif")

TDFD_percentil <- quantile(TDFD_Zonation, 0.7)

# Creating a binary raster wth 1 for cells where the value is the same or greater than percentil 70
TDFD_top30 <- calc(TDFD_Zonation, fun = function(x) { ifelse(x >= TDFD_percentil, 1, 0) })
TDFD_top30_count <- freq(TDFD_top30)

#Now the map tht combines TD, PD
TDPD_Zonation <- raster("Spatial_Data/Zonation/rankmap_TDPD.tif")

TDPD_percentil <- quantile(TDPD_Zonation, 0.7)

# Creating a binary raster wth 1 for cells where the value is the same or greater than percentil 70
TDPD_top30 <- calc(TDPD_Zonation, fun = function(x) { ifelse(x >= TDPD_percentil, 1, 0) })

#Now the map tht combines FD, PD
FDPD_Zonation <- raster("Spatial_Data/Zonation/rankmap_FDPD.tif")

FDPD_percentil <- quantile(FDPD_Zonation, 0.7)

# Creating a binary raster wth 1 for cells where the value is the same or greater than percentil 70
FDPD_top30 <- calc(FDPD_Zonation, fun = function(x) { ifelse(x >= FDPD_percentil, 1, 0) })

# Raster into data frames for ggplot2
TD_df <- as.data.frame(as(TD_top30, "SpatialPixelsDataFrame"))
FD_df <- as.data.frame(as(FD_top30, "SpatialPixelsDataFrame"))
PD_df <- as.data.frame(as(PD_top30, "SpatialPixelsDataFrame"))
TDFD_df <- as.data.frame(as(TDFD_top30, "SpatialPixelsDataFrame"))
TDPD_df <- as.data.frame(as(TDPD_top30, "SpatialPixelsDataFrame"))
FDPD_df <- as.data.frame(as(FDPD_top30, "SpatialPixelsDataFrame"))
Hotsp_df <- as.data.frame(as(Hotsp_top30, "SpatialPixelsDataFrame"))

# As factors
TD_df$layer <- factor(TD_df$layer)
FD_df$layer <- factor(FD_df$layer)
PD_df$layer <- factor(PD_df$layer)
TDFD_df$layer <- factor(TDFD_df$layer)
TDPD_df$layer <- factor(TDPD_df$layer)
FDPD_df$layer <- factor(FDPD_df$layer)
Hotsp_df$layer <- factor(Hotsp_df$layer)

# One plot per map
p1 <- ggplot() +
  geom_tile(data = TD_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_manual(values = c("white", "grey"), guide = "none") +
  geom_sf(data = AC, fill = NA, color = "black") +
  ggtitle("Prioritization TD") +
  theme_void() +
  theme(
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)  # Centrar el título
  )

p2 <- ggplot() +
  geom_tile(data = FD_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_manual(values = c("white", "goldenrod3"), guide = "none") +
  geom_sf(data = AC, fill = NA, color = "black") +
  ggtitle("Prioritization FD") +
  theme_void() +
  theme(
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)  # Centrar el título
  )

p3 <- ggplot() +
  geom_tile(data = PD_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_manual(values = c("white", "slategray3"), guide = "none") +
  geom_sf(data = AC, fill = NA, color = "black") +
  ggtitle("Prioritization PD") +
  theme_void() +
  theme(
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)  # Centrar el título
  )

p4 <- ggplot() +
  geom_tile(data = TDFD_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_manual(values = c("white", "#332288"), guide = "none") +
  geom_sf(data = AC, fill = NA, color = "black") +
  ggtitle("Prioritization TD-FD") +
  theme_void() +
  theme(
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)  # Centrar el título
  ) 

p5 <- ggplot() +
  geom_tile(data = TDPD_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_manual(values = c("white", "#DDCC77"), guide = "none") +
  geom_sf(data = AC, fill = NA, color = "black") +
  ggtitle("Prioritization TD-PD") +
  theme_void() +
  theme(
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)  # Centrar el título
  ) 

p6 <- ggplot() +
  geom_tile(data = FDPD_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_manual(values = c("white", "#CC79A7"), guide = "none") +
  geom_sf(data = AC, fill = NA, color = "black") +
  ggtitle("Prioritization FD-PD") +
  theme_void() +
  theme(
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)  # Centrar el título
  ) 

p7 <- ggplot() +
  geom_tile(data = Hotsp_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_manual(values = c("white", "darkred"), guide = "none") +
  geom_sf(data = AC, fill = NA, color = "black") +
  ggtitle("Prioritization TD-FD-PD") +
  theme_void() +
  theme(
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)  # Centrar el título
  ) +
  annotation_scale(location = "br", width_hint = 0.4, height = unit(0.5, "cm"), 
                   text_cex = 2, bar_cols = c("black", "white"), pad_x = unit(0.5, "cm"), 
                   pad_y = unit(0.5, "cm"), dist = 50, dist_unit = "km")

# Prioritization maps for TD, FD, PD, and TD-FD-PD scenarios 
png("Figures/Prioritization_scenarios.png", width = 1000, height = 1600)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, nrow = 4)
dev.off()

pdf("Figures/Prioritization_scenarios.pdf", width = 10, height = 16) 
grid.arrange(p1, p2, p3, p4, p5, p6, p7, nrow = 4) 
dev.off ()

#I want to know the percentage of individual facet scenarios covered under all seven scenarios

# First we calculate for the TD scenario

# Obviously 100% of TD scenario is covered under the TD scenario
TD_overlap <- TD_top30 * TD_top30  # This creates an overlappin raster (1 = solapamiento)
TD_coverage <- cellStats(TD_overlap, sum) / cellStats(TD_top30, sum) * 100

# Cover % for FD under the TD scenario
FD_overlap <- FD_top30 * TD_top30
FD_coverage <- cellStats(FD_overlap, sum) / cellStats(FD_top30, sum) * 100

# Cover % for PD under the TD scenario 
PD_overlap <- PD_top30 * TD_top30
PD_coverage <- cellStats(PD_overlap, sum) / cellStats(PD_top30, sum) * 100

# Results
cat("Cobertura relativa de TD:", TD_coverage, "%\n") # 100%
cat("Cobertura relativa de FD:", FD_coverage, "%\n") # 50.1%
cat("Cobertura relativa de PD:", PD_coverage, "%\n") # 33.9%

# Now we calculate for the FD scenario

# Cover % of TD under the FD scenario
TD_overlap <- TD_top30 * FD_top30  
TD_coverage <- cellStats(TD_overlap, sum) / cellStats(TD_top30, sum) * 100

# Cover % of FD under the FD scenario
FD_overlap <- FD_top30 * FD_top30
FD_coverage <- cellStats(FD_overlap, sum) / cellStats(FD_top30, sum) * 100

# Cover % of PD under the FD scenario
PD_overlap <- PD_top30 * FD_top30
PD_coverage <- cellStats(PD_overlap, sum) / cellStats(PD_top30, sum) * 100

# Results
cat("Cobertura relativa de TD:", TD_coverage, "%\n") # 50.1%
cat("Cobertura relativa de FD:", FD_coverage, "%\n") # 100%
cat("Cobertura relativa de PD:", PD_coverage, "%\n") # 54.2%

# Now under the PD scenario
# Cover % of TD under the PD scenario
TD_overlap <- TD_top30 * PD_top30  
TD_coverage <- cellStats(TD_overlap, sum) / cellStats(TD_top30, sum) * 100

# Cover % of FD under the PD scenario
FD_overlap <- FD_top30 * PD_top30
FD_coverage <- cellStats(FD_overlap, sum) / cellStats(FD_top30, sum) * 100

# Cover % of PD under the PD scenario
PD_overlap <- PD_top30 * PD_top30
PD_coverage <- cellStats(PD_overlap, sum) / cellStats(PD_top30, sum) * 100

# Results
cat("Cobertura relativa de TD:", TD_coverage, "%\n") # 33.9%
cat("Cobertura relativa de FD:", FD_coverage, "%\n") # 54.2%
cat("Cobertura relativa de PD:", PD_coverage, "%\n") # 100%

# Now under the TDFD scenario
# TD coverage% under the TDFD scenario
TD_overlap <- TD_top30 * TDFD_top30  
TD_coverage <- cellStats(TD_overlap, sum) / cellStats(TD_top30, sum) * 100

# FD coverage% under the TDFD scenario
FD_overlap <- FD_top30 * TDFD_top30
FD_coverage <- cellStats(FD_overlap, sum) / cellStats(FD_top30, sum) * 100

# PD coverage% under the TDFD scenario
PD_overlap <- PD_top30 * TDFD_top30
PD_coverage <- cellStats(PD_overlap, sum) / cellStats(PD_top30, sum) * 100

# Results
cat("Cobertura relativa de TD:", TD_coverage, "%\n") # 77.9%
cat("Cobertura relativa de FD:", FD_coverage, "%\n") # 72.2%
cat("Cobertura relativa de PD:", PD_coverage, "%\n") # 44.8%

# Now under the TDPD scenario
# TD coverage% under the TDPD scenario
TD_overlap <- TD_top30 * TDPD_top30  
TD_coverage <- cellStats(TD_overlap, sum) / cellStats(TD_top30, sum) * 100

# FD coverage% under the TDPD scenario
FD_overlap <- FD_top30 * TDPD_top30
FD_coverage <- cellStats(FD_overlap, sum) / cellStats(FD_top30, sum) * 100

# PD coverage% under the TDPD scenario
PD_overlap <- PD_top30 * TDPD_top30
PD_coverage <- cellStats(PD_overlap, sum) / cellStats(PD_top30, sum) * 100

# Results
cat("Cobertura relativa de TD:", TD_coverage, "%\n") # 74.3%
cat("Cobertura relativa de FD:", FD_coverage, "%\n") # 58.9%
cat("Cobertura relativa de PD:", PD_coverage, "%\n") # 59.6%

# Now under the FDPD scenario
# TD coverage% under the FDPD scenario
TD_overlap <- TD_top30 * FDPD_top30 
TD_coverage <- cellStats(TD_overlap, sum) / cellStats(TD_top30, sum) * 100

# FD coverage% under the FDPD scenario
FD_overlap <- FD_top30 * FDPD_top30
FD_coverage <- cellStats(FD_overlap, sum) / cellStats(FD_top30, sum) * 100

# PD coverage% under the FDPD scenario
PD_overlap <- PD_top30 * FDPD_top30
PD_coverage <- cellStats(PD_overlap, sum) / cellStats(PD_top30, sum) * 100

# Results
cat("Cobertura relativa de TD:", TD_coverage, "%\n") # 43.2%
cat("Cobertura relativa de FD:", FD_coverage, "%\n") # 81.5%
cat("Cobertura relativa de PD:", PD_coverage, "%\n") # 72.7%

# Now under the TD-FD-PD scenario
# TD coverage% under the TD-FD-PD scenario
TD_overlap <- TD_top30 * Hotsp_top30  
TD_coverage <- cellStats(TD_overlap, sum) / cellStats(TD_top30, sum) * 100

# FD coverage% under the TD-FD-PD scenario
FD_overlap <- FD_top30 * Hotsp_top30
FD_coverage <- cellStats(FD_overlap, sum) / cellStats(FD_top30, sum) * 100

# PD coverage% under the TD-FD-PD scenario
PD_overlap <- PD_top30 * Hotsp_top30
PD_coverage <- cellStats(PD_overlap, sum) / cellStats(PD_top30, sum) * 100

# Results
cat("Cobertura relativa de TD:", TD_coverage, "%\n") # 68.5%
cat("Cobertura relativa de FD:", FD_coverage, "%\n") # 75.5%
cat("Cobertura relativa de PD:", PD_coverage, "%\n") # 59.7%

# I used ArcMap to generate bivariate maps showing the combination of Biodiversity value and PV occupancy.
# In ArcMap, I loaded the two tif files: 1) rankmap_hotsp.tif (i.e., Biodiversity value represented by the TD-FD-PD prioritization scenario - generated in Zonation 5)
# and 2) PV_occupancy.tif (rasterized in ArcMap from 3facets_scaled_1.shp). 
# Next I used the reclassify tool (from Spatial Analyst) to categorize data from both tif files into 3 categories using quantiles. 
# Therefore, for each raster I had cells in 3categories: 1, 2 and 3, where 3 represents the highest values and 1 the lowest.
# Later, I used the Raster Calculator tool to generate a new raster that combines the categories from both rasters:
# For this, I used the expression: "Biodvalue" * 10 + "PV", where "Biodvalue" is the rankmap_hotsp.tif, and "PV" is the PV occupancy.
# This generated the categories 11, 12, 13, 21, 22, 23, 31, 32, 33, representing low (1), medium (2), and high (3)
# values of Biodvalue and PV occupancy. For example, 13 means low biodvalue and high PV occupancy, while 31 means high biodvalue and low PV occupancy.


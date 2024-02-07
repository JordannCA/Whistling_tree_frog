##L. verreauxii tadpole paper code 

##data set up##
pacman::p_load(DHARMa, stats, MASS, sf, raster, rasterVis, reshape2, gridExtra, lme4, rptR, pscl, glmmTMB, ggplot2, dplyr, readr, tidyverse, patchwork, epitools, cowplot, here)

#raw data files 
site_data <- read_csv("final_site_data.csv") #aggregated site level data 
data <- read_csv(here("formatted_data.csv")) #all data 
zoop <- read_csv(here("updated_zooplankton_data_site_sampled.csv")) #microfauna summary
temp <- read_csv(here("modified_ibutton_data.csv"))
temp_summary <- read_csv(here("temperature_summary.csv")) #ibutton data 
community <- read_csv(here("community_data_clean.csv")) #frog community make up 

# Create a subset of data 
tadpole_data <- subset(data, record_type == "Tadpole")
adult_data <- subset(data, record_type == "Adult")

#all age groups model Bd prevelance 
model1 <- glmer(Bd_binary ~ factor(record_type) + scale(raw_count) + scale(shannon_index) + scale(mean_temp) + scale(elevation) + (1 | site), family = binomial, data = data)

summary(model1)

#check fit 
confint(model1)
simulationOutput1 <- simulateResiduals(fittedModel= model1)
plot(simulationOutput1)

# Tadpole Bd prevelance only model 
model2 <- glmer(Bd_binary ~ scale(survey_period) + scale(elevation) + scale(mean_temp) + scale(surface_temp) + scale(raw_count) + scale(species_richness) + (1 | site), 
                family = binomial(link = "logit"), 
                data = tadpole_data,
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

summary(model2)

#check fit of model 
confint(model2)
simulationOutput2 <- simulateResiduals(fittedModel= model2)
plot(simulationOutput2)

#adult Bd prevelance only 
model3 <- glmer(Bd_binary ~ scale(amb_temp) + scale(SVL_mm) + scale(elevation) + (1 | site), 
                          family = binomial(link = "logit"), data = adult_data)
summary(model3)

#check fit of model 
confint(model3)

simulationOutput3 <- simulateResiduals(fittedModel= model3)

plot(simulationOutput3)
  
##intensity models###

#BD GEload model 
model4 <- lmer(GE_load~ factor(record_type) + scale(survey_period) + scale(raw_count) + scale(species_richness) + scale(mean_temp) + scale(elevation) + (1 | site), data = data)

summary(model4)

#check fit
confint(model4)
simulationOutput4 <- simulateResiduals(fittedModel= model4)

plot(simulationOutput4)

##very overdispersed... need a different model structure
#zero inflated model 
model4_glmmTMB <- glmmTMB(log10(GE_load +1) ~ factor(record_type) + scale(elevation) + scale(raw_count) + scale(shannon_index) + scale(mean_temp) + (1 | site), 
                          zi = ~1,  #this assumes 0 inflation is constant across alls obvservations
                          family = poisson, 
                          data = data)
summary(model4_glmmTMB)

simulationOutput5 <- simulateResiduals(fittedModel= model4_glmmTMB)

plot(simulationOutput5)

# tadpole Bd intensity only model - NOT VERY WELL FIT!!! 
model5 <- lmer(log10(GE_load +1)~ scale(survey_period) + scale(elevation) + scale(mean_temp) + scale(raw_count) + scale(species_richness) + scale(surface_temp) + (1 | site), data = tadpole_data) 

#check fit of model
confint(model5)

simulationOutput5 <- simulateResiduals(fittedModel= model5)

plot(simulationOutput5)

#adults only Bd GE load 
model6 <- lmer(log10(GE_load +1)~ scale(amb_temp) + scale(SVL_mm) + scale(elevation) + (1 | site), data = adult_data)

summary(model6)

#check fit of model
confint(model5)

simulationOutput6 <- simulateResiduals(fittedModel= model6)

plot(simulationOutput6)

# frog community analysis - 
#Crinia signifera is present at all sites, so we need to disclude. cannot use in the regression. 

#Bd proportion per site vs presence of other frog species, accounting for adult and tadpole (record type)
glmm_community <- glmer(log(Proportion_Infected +1) ~ record_type + cp + limt + limd + limp + ul + lp + (1 | site), family = binomial(link = "logit"), data = community)
summary(glmm_community)

#check fit of model 
simulationOutput <- simulateResiduals(fittedModel= glmm_community)
simulationOutput
plot(simulationOutput)

#Bd intensity (GE load) per site vs presence of other frog species, accounting for adult and tadpole (record type)
lmm_intensity <- lmer(log(Mean_Bd +1) ~ record_type + cs + cp + limt + limd + limp + ul + lp + (1 | site), data = community)

summary (lmm_intensity)

simulationOutput <- simulateResiduals(fittedModel= lmm_intensity)
simulationOutput
plot(simulationOutput)

##Glm for body size 
glm_model_body <- glm(Bd_binary ~ length_mm, family = poisson, data = tadpole_data)

lm_model_body <- lm(GE_load ~ length_mm, data = tadpole_data)

summary(glm_model_body)
summary(lm_model_body)

##abundance and diversity vs elevation for microfauna  
tadpole_data_cleaned_shannon <- na.omit(tadpole_data[c('elevation', 'shannon_index', 'site')])
glm_shannon_cleaned <- glm(shannon_index ~ elevation, data = tadpole_data_cleaned_shannon, family = gaussian())

summary(glm_shannon_cleaned)

#Negative Binomial Model- abundance 
#Log transformation, remove NAs 
tadpole_data_cleaned <- na.omit(tadpole_data[c('elevation', 'raw_count', 'site')])
tadpole_data_cleaned$log_raw_count <- log(tadpole_data_cleaned$raw_count + 1)

#Negative Binomial Model
nb_model_lme4 <- lmer(log_raw_count ~ elevation + (1 | site), data = tadpole_data_cleaned)  # Example with a random effect for 'site'

summary(nb_model_lme4)

# CHECKING FIT OF MODELS 
plot(glm_shannon_cleaned)
plot(nb_model_lme4)

##secchi data 

#Negative Binomial model
aggregated_data <- read_csv(here("aggregated_site_data.csv"))

#Gaussian GLM for Shannon Index
model_shannon <- lm(Shannon_index ~ secchi_avg, data = aggregated_data)
summary(model_shannon)

plot(model_shannon)

#Negative Binomial Model for abundance 
# Log-transform
aggregated_data$log_Total_count <- log(aggregated_data$Total_count)
model_total_count <- lm(log_Total_count ~ secchi_avg, data = aggregated_data)
summary(model_total_count)

plot(model_total_count)


###temperature data####
tdata <- read_csv(here(extracted_temp_data.csv))

# Adding a column for hours above 28 degrees
data$hours_above_28 <- ifelse(tdata$temperature > 28, 1, 0)

# Calculate the required summaries
summaries <- tdata %>%
  group_by(`depth_category`) %>%
  summarise(
    `Hours within 17-25°C` = sum(temperature >= 17 & temperature <= 25),
    `Hours below 17°C` = sum(temperature < 17),
    `Hours above 25°C` = sum(temperature > 25),
    `Hours above 28°C` = sum(hours_above_28),
    `Most Common Temperature` = as.numeric(names(sort(table(temperature), decreasing = TRUE)[1])),
    `Temperature Range` = max(temperature) - min(temperature)
  )


##visualisations### 

##BD prevelance and intensity graph ##
##summary data
categories <- c('ADULT', 'TADPOLE (EARLY)', 'TADPOLE (LATE)')
prevalence <- c(50.53, 0.82, 1.86)
lower_ci <- c(43.19, -0.1, 0.39)
upper_ci <- c(57.84, 1.75, 3.33)
errors <- data.frame(
  category = categories,
  prevalence = prevalence,
  lower = prevalence - lower_ci,
  upper = upper_ci - prevalence
)

sample_sizes <- c(190, 350, 314)
errors$sample_size <- sample_sizes

# Intensity data for boxplots
set.seed(0)
intensity <- data.frame(
  category = rep(categories, each = 100),
  intensity = c(runif(100, min = 0, max = 6870),
                runif(100, min = 0, max = 531),
                runif(100, min = 0, max = 659.7))
)

#bar graph for Bd prevalence 
p1 <- ggplot(errors, aes(x = category, y = prevalence)) +
  geom_bar(stat = 'identity', fill = 'lightgrey') +
  geom_errorbar(aes(ymin = prevalence - lower, ymax = prevalence + upper), width = 0.05) +
  labs(y = 'Prevalence (%)', x = '') +
  theme_classic() 

#boxplot for intensity 
p2 <- ggplot(intensity, aes(x = category, y = intensity)) +
  geom_boxplot(fill = 'lightgrey') +
  labs(y = 'Intensity (GE Range)', x = '') +
  theme_classic() +
  scale_y_log10()

BdPlot <- grid.arrange(p1, p2, ncol = 2)

ggsave(filename = "BdPlot.png",  
       plot = BdPlot,         
       width = 9,               
       height = 6,                  
       dpi = 300)                      


##microfauna abundance graph ##

#Aggregate total counts, per order, across all sites and survey periods
aggregate_abundance <- aggregate(Total_Count ~ Order, data = zoop, sum)

#bar graph 
microfauna <- ggplot(aggregate_abundance, aes(x = Order, y = Total_Count)) +
  geom_bar(stat = "identity", fill = "lightgrey") +
  theme_minimal() +
  labs(x = "Microfauna Order",
       y = "Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
          theme_classic()

microfauna

ggsave(filename = "MicrofaunaFigure.png",  
       plot = microfauna,           
       width = 9,                      
       height = 6,                     
       dpi = 300)                     


###per site abundance graph for microfauna ###

#Average total count across sites 
average_data <- zoop %>%
  group_by(Order, site_sampled) %>%
  summarise(Average_Count = mean(Total_Count))

#unique orders
orders <- unique(average_data$Order)

#graph colours
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


total_microfauna <- ggplot(average_data, aes(x = factor(site_sampled), y = Average_Count, fill = Order)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Site",
       y = "Average Total Count") +
  scale_fill_manual(values=cbPalette)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") # Rotate x-axis labels and adjust legend position

total_microfauna 

ggsave(filename = "TotalMicrofaunaFigure.png",  
       plot = total_microfauna,           
       width = 9,                      
       height = 6,                     
       dpi = 300)                      

##temperature per site heat map ####

#site_sampled to a factor
temp_summary$site_sampled <- as.factor(temp_summary$site_sampled)

#Reshape
df_melted <- melt(temp_summary, id.vars = c("Site", "Depth Category", "site_sampled"))

df_melted <- na.omit(df_melted)

df_melted$Category_Variable <- with(df_melted, paste(`Depth Category`, variable))

#Order factor levels
df_melted$Category_Variable <- factor(df_melted$Category_Variable, 
                                      levels = unique(df_melted$Category_Variable[order(df_melted$`Depth Category`)]))

#Graph colours 
color_palette <- brewer.pal(9, "Blues") 

heat_map <- ggplot(df_melted, aes(x = Category_Variable, y = site_sampled, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 1)), color = "black", size = 3) +
  scale_fill_gradientn(colors = color_palette, na.value = "white", name = "Total Hours") +
  theme_minimal() +
  labs(x = "Depth and Temperature Range", y = "Site") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

heat_map

ggsave(filename = "HeatMapFigure.png",  
       plot = heat_map,          
       width = 8,                      
       height = 6,                     
       dpi = 300)                     


#####Map code ######
gps_data <- read.csv(here("sites_for_map.csv"))

#Transform GPS data to spatial data
gps_points <- st_as_sf(gps_data, coords = c("lon", "lat"), crs = 4326)

#Australia SA2 boundary data
australia_sa2_shape <- st_read('~/Library/CloudStorage/OneDrive-AustralianNationalUniversity/Jordann PhD Files/whistling tree frog/tadpole_chytrid/SA2_2021_AUST_SHP_GDA2020')

#ACT regions
act_regions <- australia_sa2_shape %>% 
  filter(STE_NAME21 == "Australian Capital Territory")

#Dissolve all ACT regions into a single boundary
act_boundary <- act_regions %>%
  summarize(geometry = st_union(geometry)) %>%
  st_cast("POLYGON")

elevation_data <- raster("~/Library/CloudStorage/OneDrive-AustralianNationalUniversity/Jordann PhD Files/whistling tree frog/tadpole_chytrid/5m_DEM.tif")

#Check CRS of all layers
crs_act_boundary <- st_crs(act_boundary)
crs_gps_points <- st_crs(gps_points)
crs_elevation_data <- crs(elevation_data)

crs_act_boundary
crs_gps_points
crs_elevation_data

#Transform ACT boundary and GPS points to match elevation data CRS
act_boundary_transformed <- st_transform(act_boundary, crs_elevation_data)
gps_points_transformed <- st_transform(gps_points, crs_elevation_data)

#Define extent for ACT
extent_of_interest <- extent(c(148.7, 149.4, -36, -35.1))

#Clip raster to extent of interest
elev_raster_clipped <- crop(elevation_data, act_boundary_transformed)

#aggregation factor
aggregation_factor <- 5

#Aggregate raster
elev_raster_aggregated <- aggregate(elev_raster_clipped, fact = aggregation_factor, fun = mean)

#Convert elevation data to same CRS as the ACT boundary
elevation_data_transformed <- projectRaster(elev_raster_aggregated, crs = crs(act_boundary))

#Mask
elevation_data_masked <- mask(elevation_data_transformed, act_boundary)

#convert to dataframe
elevation_df <- raster::rasterToPoints(elevation_data_masked)
elevation_df <- data.frame(elevation_df)

gps_points_transformed$Sampling_Status <- as.factor(gps_points_transformed$Sampling_Status)

### final map ####
#boundary limits
x_min <- 148.7
x_max <- 149.4
y_min <- -36
y_max <- -35.1

#spread out map
padding <- 0.01 


main_map <- ggplot() +
  geom_raster(data = elevation_df, aes(x = x, y = y, fill = X5m_DEM), alpha = 0.7) +
  scale_fill_viridis_c(name = "Elevation", guide = guide_colourbar(direction = "horizontal", title.position = "top", title.hjust = 0.5, position = "bottom", hjust = 0)) +
  geom_sf(data = gps_points_transformed, aes(color = Sampling_Status), shape = 20, size = 2, stroke = 0.5) +
  geom_sf(data = act_boundary, fill = NA, color = "black") +
  theme_bw() +
  labs(
    x = "Longitude",  
    y = "Latitude",   
    color = "Age Group Sampled"
  ) +
  coord_sf(crs = crs_elevation_data, xlim = c(x_min - padding, x_max + padding), ylim = c(y_min - padding, y_max + padding)) +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )

main_map

### additions to the map ###
library(ggplot2)
library(sf)
library(rnaturalearth)
library(ggspatial)

main_map_with_scale <- main_map +
  annotation_scale(location = "bl", width_hint = 0.3) + 
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering)

main_map_with_scale

#inlay map 
world <- ne_countries(scale = "medium", returnclass = "sf")
australia <- world[world$admin == "Australia", ]

#check projection 
crs_australia <- st_crs(australia)
crs_act_boundary <- st_crs(act_boundary)

##transform CRS
australia_transformed <- st_transform(australia, crs = crs_act_boundary)

australia_act_map <- ggplot() +
  geom_sf(data = australia_transformed, fill = NA, color = "black") +
  geom_sf(data = act_boundary_transformed, fill = "black", color = "black") +
  theme_void()

australia_act_map

##bounding box to show ACT 
australia_map <- australia_act_map +
  geom_rect(
    aes(
      xmin = 150, xmax = 148, 
      ymin = -36, ymax = -35),
    color = "red", fill = NA, size = 1) +
  theme_void()

australia_map

library(cowplot)

final_map_figure <- ggdraw() +
  draw_plot(main_map_with_scale) +
  draw_plot(australia_map, x = 0.40, y = 0.01, width = 0.3, height = 0.3)

final_map_figure

#### NEED TO CENTRE RIGHT FOR EXPORT!!!! 
####
ggsave(filename = "MapFigure.png",  # or .pdf, .svg, .eps for vectorized graphics
       plot = final_map_figure,           # your ggplot object
       width = 7.5,                      # Width of the plot in inches
       height = 10,                     # Height of the plot in inches
       dpi = 300)                      # Resolution in dots per inch, for raster formats



###Temperature Graph ####

#temp already imported
ibutton_data <- temp

# Convert the date/time variable to a POSIXct format
ibutton_data$date <- as.POSIXct(ibutton_data$date, format = "%d/%m/%Y %H:%M")

#separate date and time
ibutton_data$date <- as.Date(ibutton_data$date)
ibutton_data$time <- format(ibutton_data$date, "%H:%M:%S")

#daily mean temperatures for each site
daily_mean_temp <- ibutton_data %>%
  group_by(date, site_sampled) %>%
  summarize(mean_temp = mean(temperature, na.rm = TRUE))

#site to a factor for labelling
daily_mean_temp$site_sampled <- factor(daily_mean_temp$site_sampled)

##remove the outliers if needed##
daily_mean_temp <- daily_mean_temp %>%
  filter(mean_temp <= 35)

dailymeantemp_plot <- ggplot() +
  geom_line(data = daily_mean_temp, aes(x = date, y = mean_temp, color = site_sampled)) +
  scale_color_discrete(name = "Site") +
  labs(x = "Date", y = "Temperature (°C)") +
  ggtitle("Daily Mean Temperatures (°C) by Site") +
  geom_hline(yintercept = 25, color = "red", linetype = 2, linewidth = 1.0) +
  geom_hline(yintercept = 17, color = "red", linetype = 2, linewidth = 1.0) +
  geom_hline(yintercept = 28, color = "blue", linetype = 2, linewidth = 1.0) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right") + theme(axis.title = element_text(size = 10), axis.text = element_text(size = 10))

dailymeantemp_plot

# Filter Bd_records for records with record_type "Tadpole"
tadpole_records <- tadpole_data %>%
  filter(record_type == "Tadpole", !is.na(ID))

#tadpole temps plot 
tadpole_plot <- ggplot(tadpole_records, aes(x = date, y = surface_temp, color = Bd)) +
  geom_point()+
  geom_hline(yintercept = 25, color = "red", linetype = 2, linewidth = 1.5) +
  geom_hline(yintercept = 17, color = "red", linetype = 2, linewidth = 1.5) +
  geom_hline(yintercept = 28, color = "blue", linetype = 2, linewidth = 1.0) +
  labs(x = "Date", y = "Temperature (°C)") +
  ggtitle("Temperatures of Tadpoles at Time of Capture")

tadpole_plot
```

#########FINAL TEMP PLOT #######
#combined plot 

#temperature column in both data frames needs to have same name
daily_mean_temp$temperature <- daily_mean_temp$mean_temp
tadpole_records$temperature <- tadpole_records$surface_temp

#remove the original columns
daily_mean_temp$mean_temp <- NULL
tadpole_records$surface_temp <- NULL

#ensure all same date type
tadpole_records$date <- as.Date(tadpole_records$date, format = "%d/%m/%Y")

#check factors 
tadpole_records$site_sampled <- as.factor(tadpole_records$site_sampled)
tadpole_records$elevation <- as.factor(tadpole_records$elevation)

#identify source of data then combine two dataframes 
daily_mean_temp$source <- "MeanTemp"
tadpole_records$source <- "Tadpole"

combined_data <- rbind(daily_mean_temp, tadpole_records)

#sort tadpole data into postive and negative to see all positive data points  
tadpole_negatives <- subset(combined_data, source == "Tadpole" & Bd == "Negative")
tadpole_positives <- subset(combined_data, source == "Tadpole" & Bd == "Positive")

#make sure tadpoles are above mean temp 
combined_data$source <- factor(combined_data$source, levels = c("Tadpole", "MeanTemp"))

#graph colours
colors_for_sites <- c(
  "#e6194B", # Red
  "#3cb44b", # Green
  "#ffe119", # Yellow
  "#4363d8", # Blue
  "#f58231", # Orange
  "#911eb4", # Purple
  "#46f0f0", # Cyan
  "#f032e6", # Magenta
  "#bcf60c", # Lime
  "#fabebe", # Pink
  "#008080", # Teal
  "#e6beff", # Lavender
  "#9a6324", # Brown
  "#fffac8", # Beige
  "#800000", # Maroon
  "#aaffc3", # Mint
  "#808000"  # Olive
)

colors_for_Bd <- c("darkgrey", "red2")


combined_plot <- ggplot(combined_data, aes(x = date, y = temperature)) +
  geom_line(data = subset(combined_data, source == "MeanTemp"), aes(color = site_sampled)) +
  geom_point(data = tadpole_negatives, aes(color = Bd)) +
  geom_point(data = tadpole_positives, aes(color = Bd)) +
  geom_hline(yintercept = 25, color = "blue", linetype = 2, linewidth = 0.5) +
  geom_hline(yintercept = 17, color = "blue", linetype = 2, linewidth = 0.5) +
  geom_hline(yintercept = 28, color = "red", linetype = 2, linewidth = 0.5) +
  facet_grid(source ~ ., scales = "free_y") +
  scale_color_manual(name= "Sites and Bd status", values = c(colors_for_sites, colors_for_Bd)) +
  labs(x = "Date", y = "Temperature (°C)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) 


combined_plot

ggsave(filename = "TempFigure.png",  
       plot = combined_plot,           
       width = 8,                     
       height = 6,                    
       dpi = 300)                      

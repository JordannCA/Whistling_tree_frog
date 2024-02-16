##L. verreauxii tadpole paper code 

##data set up##
pacman::p_load(DHARMa, stats, MASS, rnaturalearth, ggspatial, sf, raster, rasterVis, reshape2, gridExtra, lme4, rptR, pscl, glmmTMB, ggplot2, dplyr, readr, tidyverse, patchwork, epitools, cowplot, here)

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

## Bd Prevelance Models ## 
#all age groups Bd prevelance 
model1 <- glmer(Bd_binary ~ factor(record_type) + scale(raw_count) + scale(shannon_index) + scale(mean_temp) + scale(elevation) + (1 | site), family = binomial(link = "logit"), data = data)

summary(model1)

#check fit 
confint(model1)
simulationOutput1 <- simulateResiduals(fittedModel= model1)
plot(simulationOutput1)

# Tadpole Bd prevalence only model 
model2 <- glmer(Bd_binary ~ scale(survey_period) + scale(elevation) + scale(mean_temp) + scale(length_mm) + scale(raw_count) + scale(species_richness) + (1 | site), 
                family = binomial(link = "logit"), 
                data = tadpole_data,
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

summary(model2)

#check fit of model 
confint(model2)
simulationOutput2 <- simulateResiduals(fittedModel= model2)
plot(simulationOutput2)

# Adult Bd prevalence only 
model3 <- glmer(Bd_binary ~ scale(amb_temp) + scale(SVL_mm) + scale(elevation) + (1 | site), 
                          family = binomial(link = "logit"), data = adult_data)
summary(model3)

#check fit of model 
confint(model3)

simulationOutput3 <- simulateResiduals(fittedModel= model3)

plot(simulationOutput3)

## body size for tadpoles ## 

Model4 <- glmer(Bd_binary ~ scale(survey_period) + scale(length_mm) + (1 | site), 
                family = binomial(link = "logit"), 
                data = tadpole_data,
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

summary(Model4)

simulationOutput4 <- simulateResiduals(fittedModel= Model4)
plot(simulationOutput4)

  
## adult intensity model ###

#adults Bd GE load 

model5_glmmTMB <- glmmTMB(GE_load ~ scale(amb_temp) + scale(SVL_mm) + scale(elevation) + (1 | site), 
                          zi = ~1, 
                          family = nbinom2(), 
                          data = adult_data)
summary(model5_glmmTMB)

#check fit of model
simulationOutput5 <- simulateResiduals(fittedModel= model5_glmmTMB)

plot(simulationOutput5)

## frog community analysis ## 

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

#check fit
simulationOutput <- simulateResiduals(fittedModel= lmm_intensity)
simulationOutput
plot(simulationOutput)

##abundance and diversity vs elevation for microfauna ##

tadpole_data_cleaned_shannon <- na.omit(tadpole_data[c('elevation', 'shannon_index', 'site')])
glm_shannon_cleaned <- glm(shannon_index ~ elevation, data = tadpole_data_cleaned_shannon, family = gaussian())

summary(glm_shannon_cleaned)

#Negative Binomial Model- abundance 

#Log transform, remove NAs 
tadpole_data_cleaned <- na.omit(tadpole_data[c('elevation', 'raw_count', 'site')])
tadpole_data_cleaned$log_raw_count <- log(tadpole_data_cleaned$raw_count + 1)

#Negative Binomial Model
nb_model_lme4 <- lmer(log_raw_count ~ elevation + (1 | site), data = tadpole_data_cleaned)  # Example with a random effect for 'site'

summary(nb_model_lme4)

#fit
plot(glm_shannon_cleaned)
plot(nb_model_lme4)

##secchi data ##

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


## temperature data summary ##

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

#**********************************************************************************************************
# Paper: Uncovering correlates of decline and critical refuges for a threatened terrestrial mammal
# Journal: Conservation Biology
# Date: May 27th 2026
# Purpose: Code for core components of analysis
#**********************************************************************************************************

##### Exploration of 3000 random sampling points ------------------------------
# Import DFs and combine --------------------------------------------------------------
ext_alpha_points <- read.csv("your/path/here/extirpated_data.csv")
contemp_alpha_points<-read.csv("your/path/here/contemporary_data.csv")


cols_to_factor <- c("land_clear", "inapp_fire", "geology", "vegetation", "grazing")
set.seed(12345)

contemp <- contemp_alpha_points %>%
  filter(!apply(. == -9999, 1, any)) %>%
  select(-OID_, -CID) %>%
  sample_n(3000) %>%
  mutate(across(all_of(cols_to_factor), as.factor), Status = 0)

ext <- ext_alpha_points %>%
  filter(!apply(. == -9999, 1, any)) %>%
  select(-OID_, -CID) %>%
  sample_n(3000) %>%
  mutate(across(all_of(cols_to_factor), as.factor), Status = 1)

combined_data <- bind_rows(contemp, ext)
rm(contemp,contemp_alpha_points,ext, ext_alpha_points, cols_to_factor)

head(combined_data)
# Comparing variables across time periods --------------------------------------------------------------

library(vcd)
library(dplyr)
library(effsize)
library(rcompanion)

categorical_vars <- c("land_clear", "inapp_fire", "geology", "vegetation", "grazing")
numerical_vars <- c("water", "slope", "rugged", "productivity", "elevation", "min_temp", 
                      "avg_rain", "aspect", "rabbit", "pig", "temp_diff", "lantana", "horse", 
                      "habitat", "goat", "fox_den", "flood", "drought", "dog", "cat_den")

ext_data <- subset(combined_data, Status == 1)
contemp_data <- subset(combined_data, Status == 0)

#****************************************************************************
# Categorical: Chi-square / Fisher + Cramer's V
#****************************************************************************

cat_results <- data.frame()
for (var in categorical_vars) {
  tab <- table(combined_data[[var]], combined_data$Status)
  test <- chisq.test(tab)
  
  cat_results <- rbind(cat_results, data.frame(
    Variable  = var,
    Test = "Chi-sq",
    P_value = round(test$p.value, 4),
    Cramers_V = round(as.numeric(rcompanion::cramerV(tab)), 2)
  ))
}
print(cat_results)

#****************************************************************************
# Numerical: Wilcoxon + Cliff's Delta
#****************************************************************************

calculate_ci <- function(x, conf = 0.95) {
  m  <- median(x)
  se <- 1.253 * (sd(x) / sqrt(length(x)))
  m + c(-1, 1) * qnorm((1 + conf) / 2) * se
}

num_results <- data.frame()

for (var in numerical_vars) {
  ext_var<- na.omit(ext_data[[var]])
  contemp_var <- na.omit(contemp_data[[var]])
  
  ext_ci <- calculate_ci(ext_var)
  contemp_ci <- calculate_ci(contemp_var)
  wilcox <- wilcox.test(ext_var, contemp_var, exact = FALSE)
  delta <- effsize::cliff.delta(ext_var, contemp_var)$estimate
  
  num_results <- rbind(num_results, data.frame(
    Variable  = var,
    Extirpated_Median = round(median(ext_var), 2),
    Extirpated_CI = sprintf("[%.2f - %.2f]", ext_ci[1], ext_ci[2]),
    Contemporary_Median = round(median(contemp_var), 2),
    Contemporary_CI = sprintf("[%.2f - %.2f]", contemp_ci[1], contemp_ci[2]),
    P_value = round(wilcox$p.value, 8),
    Cliffs_Delta = round(delta, 2)
  ))
}

num_results$Variable <- gsub("dog", "dingo", num_results$Variable)
num_results <- num_results %>% arrange(desc(Cliffs_Delta))

print(num_results)
##### GAMM with 400 random sampling points ------------------------------
# Import DFs and combine --------------------------------------------------------------
rm(list=ls())

ext_alpha_points <- read.csv("your/path/here/extirpated_data.csv")
contemp_alpha_points<-read.csv("your/path/here/contemporary_data.csv")


cols_to_factor <- c("land_clear", "inapp_fire", "geology", "vegetation", "grazing")
set.seed(12345)

contemp <- contemp_alpha_points %>%
  filter(!apply(. == -9999, 1, any)) %>%
  select(-OID_, -CID) %>%
  sample_n(400) %>%
  mutate(across(all_of(cols_to_factor), as.factor), Status = 0)

ext <- ext_alpha_points %>%
  filter(!apply(. == -9999, 1, any)) %>%
  select(-OID_, -CID) %>%
  sample_n(400) %>%
  mutate(across(all_of(cols_to_factor), as.factor), Status = 1)

combined_data <- bind_rows(contemp, ext)
rm(contemp,ext, cols_to_factor)

# Final GAMM --------------------------------------------------------------
library(mgcv)
start_time <- Sys.time()

gamm_final <- gamm (Status ~ water + s(min_temp, k = 5) + s(rugged, k = 5) + inapp_fire + 
                             s(goat, k = 5) + s(fox_den, k = 5) + s(dog, k = 5) + s(cat_den,k = 5) + 
                             ti(fox_den, rugged),
                                data = combined_data, 
                                family = binomial(), 
                                method = "REML",
                                niterPQL=40,
                                correlation = corExp(value=c(10000), form = ~ lat + long, 
                                                     nugget = TRUE))

end_time <- Sys.time()
print(end_time -start_time)

summary(gamm_final$gam)
gam.check(gamm_final$gam) #no significant
concurvity(gamm_final$gam) # below 0.8 for worst

# Threshold for 80% probability of persistence -------------------------------------------
#Illustrated using fox density below 

combined_data$dingo <- combined_data$dog

predict_prob <- function(model, combined_data, fox_den_value) {
  combined_data$fox_den <- fox_den_value
  combined_data$water <- median(combined_data$water, na.rm = TRUE)
  combined_data$min_temp <- median(combined_data$min_temp, na.rm = TRUE)
  combined_data$rugged <- median(combined_data$rugged, na.rm = TRUE)
  combined_data$goat <- median(combined_data$goat, na.rm = TRUE)
  combined_data$dingo <- median(combined_data$dingo, na.rm = TRUE)
  combined_data$cat_den <- median(combined_data$cat_den, na.rm = TRUE)
  combined_data$inapp_fire<- as.factor(names(which.max(table(combined_data$inapp_fire))))
  mean(predict(model, newdata = combined_data, type = "response"))
}
fox_den_values <- seq(min(combined_data$fox_den, na.rm = TRUE), 
                      max(combined_data$fox_den, na.rm = TRUE), 
                      length.out = 10000)

probabilities <- sapply(fox_den_values, function(x) predict_prob(gamm_final$gam, combined_data, x))

# Estimated threshold for foxes
closest_index <- which.min(abs(probabilities - 0.2))
fox_den_values[closest_index]



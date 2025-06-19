################################################################################
#Estevão Cardoso - Pset 2#######################################################

#Importing the libraries we need ###############################################

library(MetBrewer)
library(forecast)
library(dplyr)
library(tidyr)
library(tibble)
library(xtable)
library(ggplot2)
library(gridExtra)
library(tseries)
library(quantmod)
library(doParallel)
library(gganimate)
library(transformr)
library(tictoc) 
library(urca)  
library(lmtest) 
library(car)  
library(lpdensity)  
library(kableExtra)
library(stargazer)
library(modelsummary)
library(glue)
library(scales)
library(geomtextpath)

################################################################################
######################Question 2################################################
rm(list = ls())
set.seed(666)

################################################################################
##############################Item a############################################
# I did more than the 1 and 5 Degrees of Freedoms. But if you want to do only
# them, it is just delete or run the part of the code that has the others DF.

#parameters
alpha = 0  # intercept
delta = 1  # coefficient of time trend
size = .1  # test size


# Adjusting the MC  
sample_size = 10000  # sample size
number_sim = 10000  # number of simulations

data = data.frame(t = seq(1, sample_size, by = 1))

# Creating a XtXt' matrix and its inverse to manually execute the test

XtXt = matrix(c(sample_size, 
                1/2 * sample_size * (sample_size + 1),
                1/2 * sample_size * (sample_size + 1),
                1/6 * sample_size * (sample_size + 1) * (2*sample_size + 1)
), nrow = 2, ncol = 2)

XtXt_inv = solve(XtXt)


#Setting the storage to place the simulations----------------------------------- 

df_mc_1t_results <- data.frame(
  "estimate" = rep(NA, number_sim),
  "bias" = rep(NA, number_sim),
  "t_stat" = rep(NA, number_sim),
  "p_value" = rep(NA, number_sim),
  "reject" = rep(NA, number_sim)
)

df_mc_2t_results <- data.frame(
  "estimate" = rep(NA, number_sim),
  "bias" = rep(NA, number_sim),
  "t_stat" = rep(NA, number_sim),
  "p_value" = rep(NA, number_sim),
  "reject" = rep(NA, number_sim)
)


df_mc_5t_results <- data.frame(
  "estimate" = rep(NA, number_sim),
  "bias" = rep(NA, number_sim),
  "t_stat" = rep(NA, number_sim),
  "p_value" = rep(NA, number_sim),
  "reject" = rep(NA, number_sim)
)

df_mc_3t_results <- data.frame(
  "estimate" = rep(NA, number_sim),
  "bias" = rep(NA, number_sim),
  "t_stat" = rep(NA, number_sim),
  "p_value" = rep(NA, number_sim),
  "reject" = rep(NA, number_sim)
)

df_mc_4t_results <- data.frame(
  "estimate" = rep(NA, number_sim),
  "bias" = rep(NA, number_sim),
  "t_stat" = rep(NA, number_sim),
  "p_value" = rep(NA, number_sim),
  "reject" = rep(NA, number_sim)
)

df_mc_10t_results <- data.frame(
  "estimate" = rep(NA, number_sim),
  "bias" = rep(NA, number_sim),
  "t_stat" = rep(NA, number_sim),
  "p_value" = rep(NA, number_sim),
  "reject" = rep(NA, number_sim)
)


## Monte Carlo -----------------------------------------------------------------
tic()
for (sim in 1:number_sim){
  # Sampling errors from the desired t-distributions
  data %<>% 
    mutate(
      t_1df = rt(sample_size, df = 1),
      t_2df = rt(sample_size, df = 2),
      t_3df = rt(sample_size, df = 3),
      t_4df = rt(sample_size, df = 4),
      t_5df = rt(sample_size, df = 5),
      t_10df = rt(sample_size, df = 10),
      data_1df = alpha + delta*t + t_1df,
      data_2df = alpha + delta*t + t_2df,
      data_3df = alpha + delta*t + t_3df,
      data_4df = alpha + delta*t + t_4df,
      data_5df = alpha + delta*t + t_5df,
      data_10df = alpha + delta*t + t_10df
    )
  
  # Estimating with OLS
  model_1df = lm(data_1df ~ t, data = data)
  model_2df = lm(data_2df ~ t, data = data)
  model_3df = lm(data_3df ~ t, data = data)
  model_4df = lm(data_4df ~ t, data = data)
  model_5df = lm(data_5df ~ t, data = data)
  model_10df = lm(data_10df ~ t, data = data)
  
  # Estimate of error variance
  sT2_1df = sum(model_1df$residuals^2) / (sample_size - 2)
  sT2_2df = sum(model_2df$residuals^2) / (sample_size - 2)
  sT2_3df = sum(model_3df$residuals^2) / (sample_size - 2)
  sT2_4df = sum(model_4df$residuals^2) / (sample_size - 2)
  sT2_5df = sum(model_5df$residuals^2) / (sample_size - 2)
  sT2_10df = sum(model_10df$residuals^2) / (sample_size - 2)
  
  # t statistics
  t_1df = (model_1df$coefficients[[2]] - delta) / sqrt(sT2_1df * (c(0,1) %*% XtXt_inv %*% c(0,1))[1,1])
  t_2df = (model_2df$coefficients[[2]] - delta) / sqrt(sT2_2df * (c(0,1) %*% XtXt_inv %*% c(0,1))[1,1])
  t_3df = (model_3df$coefficients[[2]] - delta) / sqrt(sT2_3df * (c(0,1) %*% XtXt_inv %*% c(0,1))[1,1])
  t_4df = (model_4df$coefficients[[2]] - delta) / sqrt(sT2_4df * (c(0,1) %*% XtXt_inv %*% c(0,1))[1,1])
  t_5df = (model_5df$coefficients[[2]] - delta) / sqrt(sT2_5df * (c(0,1) %*% XtXt_inv %*% c(0,1))[1,1])
  t_10df = (model_10df$coefficients[[2]] - delta) / sqrt(sT2_10df * (c(0,1) %*% XtXt_inv %*% c(0,1))[1,1])
  
  # Two-sided P-values from N(0, 1)
  p_1df = 2 * min(pnorm(t_1df, lower.tail = TRUE), pnorm(t_1df, lower.tail = FALSE))
  p_2df = 2 * min(pnorm(t_2df, lower.tail = TRUE), pnorm(t_2df, lower.tail = FALSE))
  p_3df = 2 * min(pnorm(t_3df, lower.tail = TRUE), pnorm(t_3df, lower.tail = FALSE))
  p_4df = 2 * min(pnorm(t_4df, lower.tail = TRUE), pnorm(t_4df, lower.tail = FALSE))
  p_5df = 2 * min(pnorm(t_5df, lower.tail = TRUE), pnorm(t_5df, lower.tail = FALSE))
  p_10df = 2 * min(pnorm(t_10df, lower.tail = TRUE), pnorm(t_10df, lower.tail = FALSE))
  
  # Storing results
  df_mc_1t_results$estimate[sim] = model_1df$coefficients[[2]]
  df_mc_1t_results$bias[sim] = model_1df$coefficients[[2]] - delta
  df_mc_1t_results$t_stat[sim] = t_1df
  df_mc_1t_results$p_value[sim] = p_1df
  df_mc_1t_results$reject[sim] = as.numeric(p_1df < size)
  
  df_mc_2t_results$estimate[sim] = model_2df$coefficients[[2]]
  df_mc_2t_results$bias[sim] = model_2df$coefficients[[2]] - delta
  df_mc_2t_results$t_stat[sim] = t_2df
  df_mc_2t_results$p_value[sim] = p_2df
  df_mc_2t_results$reject[sim] = as.numeric(p_2df < size)
  
  df_mc_3t_results$estimate[sim] = model_3df$coefficients[[2]]
  df_mc_3t_results$bias[sim] = model_3df$coefficients[[2]] - delta
  df_mc_3t_results$t_stat[sim] = t_3df
  df_mc_3t_results$p_value[sim] = p_3df
  df_mc_3t_results$reject[sim] = as.numeric(p_3df < size)
  
  df_mc_4t_results$estimate[sim] = model_4df$coefficients[[2]]
  df_mc_4t_results$bias[sim] = model_4df$coefficients[[2]] - delta
  df_mc_4t_results$t_stat[sim] = t_4df
  df_mc_4t_results$p_value[sim] = p_4df
  df_mc_4t_results$reject[sim] = as.numeric(p_4df < size)
  
  df_mc_5t_results$estimate[sim] = model_5df$coefficients[[2]]
  df_mc_5t_results$bias[sim] = model_5df$coefficients[[2]] - delta
  df_mc_5t_results$t_stat[sim] = t_5df
  df_mc_5t_results$p_value[sim] = p_5df
  df_mc_5t_results$reject[sim] = as.numeric(p_5df < size)
  
  df_mc_10t_results$estimate[sim] = model_10df$coefficients[[2]]
  df_mc_10t_results$bias[sim] = model_10df$coefficients[[2]] - delta
  df_mc_10t_results$t_stat[sim] = t_10df
  df_mc_10t_results$p_value[sim] = p_10df
  df_mc_10t_results$reject[sim] = as.numeric(p_10df < size)
}
toc()

df_results_list <- list(
  "1" = df_mc_1t_results,
  "2" = df_mc_2t_results,
  "3" = df_mc_3t_results,
  "4" = df_mc_4t_results,
  "5" = df_mc_5t_results,
  "10" = df_mc_10t_results
)
colors <- met.brewer("Homer2")[c(1, 2, 3, 4, 5, 10)]
dfs <- c(1, 2, 3, 4, 5, 10)


################################Rejection rates#################################



rejection_rates <- list()

for (df_name in names(df_results_list)) {
  
  
  df_value <- gsub("df_", "", df_name)
  df_data <- df_results_list[[df_name]]
  
  
  rejection_rate <- round(100 * mean(df_data$reject, na.rm = TRUE), 2)
  
  
  rejection_rates[[df_name]] <- data.frame(
    df = df_value,
    rejection_rate = rejection_rate
  )
}


rejection_rates_df <- bind_rows(rejection_rates)

print(
  xtable(rejection_rates_df, 
         caption = "Rejection rates by degrees of freedom", 
         label = "tab:rejection_rates"),
  include.rownames = FALSE
)

##########Comparing rejec rate nrmal dist with a uniform items a and b##########

size <- 0.1
density_df_list <- list()
for (df_name in names(df_results_list)) {
  
  df_data <- df_results_list[[df_name]]
  
  # Estimate p-value density
  p_density <- lpdensity(
    data = df_data$p_value,
    grid = quantile(df_data$p_value, probs = seq(0.01, 0.99, by = 0.01), na.rm = TRUE)
  )
  
  # Rejection rate
  rejection_rate <- round(100 * mean(df_data$reject, na.rm = TRUE), 2)
  
  # Build dataframe for plotting
  p_density_plot <- data.frame(
    grid = p_density$Estimate[, 1],
    p_density = p_density$Estimate[, 5],
    df = paste0("t(", df_name, ")"),
    rejection_rate = rejection_rate
  )
  
  # Append to list
  density_df_list[[df_name]] <- p_density_plot
}

dfs <- names(df_results_list)
df_labels <- paste0("t(", dfs, ")") 
uniform_density_df <- data.frame(
  x = rep(seq(0, 1, length.out = 1000), times = length(dfs)),
  y = 1,
  df = rep(df_labels, each = 1000)
)

#Combine all into one data frame

full_density_df <- bind_rows(density_df_list)

decision_rate_all <- 
  ggplot(full_density_df, aes(x = grid, y = p_density, color = df)) +
  geom_line(size = 1.2) +
  geom_line(
    data = uniform_density_df,
    aes(x = x, y = y, color = "Uniform(0,1)"),
    linetype = "dashed", size = 0.7, inherit.aes = FALSE
  ) +
  geom_vline(xintercept = size, linetype = "dotted", size = 0.8) +
  facet_wrap(~ df, ncol = 3) +
  labs(x = "P-Values", y = "Density", color = "") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual(
    values = c(setNames(colors, df_labels), "Uniform(0,1)" = "black")
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 13),
    legend.text = element_text(size = 13),
    strip.background = element_blank()
  )

ggsave("decision_rate_all.png", plot = decision_rate_all, height = 7, width = 10, dpi = 300)



#############builting the t distribution with DF items a and b##################

densidades_t <- list()

for (df in dfs) {
  result_df <- results_list[[paste0("df_", df)]]
  
  grid_df <- quantile(result_df$t_stat, probs = seq(0.01, 0.99, 0.01), na.rm = TRUE)
  
  densidade <- lpdensity(data = result_df$t_stat, grid = grid_df)
  
  densidade_df <- data.frame(
    grid = densidade$Estimate[, 1],
    t_density = densidade$Estimate[, 5],
    df = paste0("df = ", df)
  )
  
  densidades_t[[as.character(df)]] <- densidade_df
}

df_densidades_plot <- bind_rows(densidades_t)



# Creating the plots for all densities
t_density_all_dfs <- ggplot(df_densidades_plot, aes(x = grid, y = t_density, color = df)) +
  geom_line(size = 1.2) +
  # Normal(0,1) density line for reference — add separately for each facet
  geom_line(data = data.frame(x = seq(-5, 5, length.out = 1000), 
                              y = dnorm(seq(-5, 5, length.out = 1000)),
                              df = rep(unique(df_densidades_plot$df), each = 1000)),
            aes(x = x, y = y, color = "Normal(0,1)"), 
            linetype = "dashed", size = 0.5, inherit.aes = FALSE) +
  scale_color_manual(values = c(setNames(colors, paste0("df = ", dfs)), "Normal(0,1)" = "black")) +
  labs(x = "T-Statistics", y = "Density", color = "") +
  facet_wrap(~ df, scales = "free_y", ncol = 3) +  
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    strip.background = element_blank()
  )


ggsave("t_density_all_df.png", plot = t_density_all_dfs, height = 8, width = 12, dpi = 300)





################################Question 4######################################
################################################################################

#ploting the graphs for prod and land-------------------------------------------

df_us <- df %>%
  filter(Entity == "United States", Year >= 1950, Year <= 2021) %>%
  select(Year, production = `Corn production (tonnes)`,
         area = `Corn, area harvested (hectares)`) %>%
  drop_na()

head(df_us)

p_production <- ggplot(df_us, aes(x = Year, y = production / 1e6)) +
  geom_line(color = "#1f78b4", size = 1) +
  labs(
    #title = "Corn Production in the U.S. (1950–2021)",
    x = "Year",
    y = "Production (millions of tonnes)"
  ) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

ggsave("us_corn_production.png", plot = p_production, height = 6, width = 10, dpi = 300)

p_area <- ggplot(df_us, aes(x = Year, y = area / 1e6)) +
  geom_line(color = "#FF4040", size = 1) +
  labs(
    #title = "Corn Cultivated Area in the U.S. (1950–2021)",
    x = "Year",
    y = "Area (millions of hectares)"
  ) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

ggsave("us_corn_area.png", plot = p_area, height = 6, width = 10, dpi = 300)




#running the model and testing the optimal number of lags in the AR-------------

serie_ts <- ts(df_us$production, start = 1950, end = 2021, frequency = 1)
teste_adf_auto <- ur.df(serie_ts, type = "trend", selectlags = "BIC")

summary(teste_adf_auto)


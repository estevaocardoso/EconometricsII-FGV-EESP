################################################################################
#Estevão Cardoso - Pset 1#######################################################

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


################################################################################
######################Question 1################################################

#We clean and import the dataset################################################

pset_1_gdp <- read.csv("data_gdp_brazil.csv")
colnames(pset_1_gdp) <- c("year", "gdp")

head(pset_1_gdp)

gdp_ts <- ts(pset_1_gdp$gdp, frequency = 1)

###########################Question 1 items a) and b)###########################

###########################Setting the models###################################

#Setting a dataframe to store the models########################################
model_list <- list(
  "AR(1)"     = c(1, 0, 0),
  "AR(2)"     = c(2, 0, 0),
  "MA(1)"     = c(0, 0, 1),
  "MA(2)"     = c(0, 0, 2),
  "ARMA(1,1)" = c(1, 0, 1),
  "ARMA(2,1)" = c(2, 0, 1),
  "ARMA(1,2)" = c(1, 0, 2),
  "ARMA(2,2)" = c(2, 0, 2)
)


results_list <- list()


model_params <- function(model, name) {
  coefs <- coef(model)
  se <- sqrt(diag(vcov(model)))
  zval <- coefs / se
  pval <- 2 * (1 - pnorm(abs(zval)))
  
  info <- list()
  for (param in names(coefs)) {
    info[[param]] <- coefs[param]
    info[[paste0("std_error_", param)]] <- se[param]
    info[[paste0("pval_", param)]] <- pval[param]
  }
  
  info$AIC <- AIC(model)
  info$BIC <- BIC(model)
  
  return(info)
}

#Estimate all models
for (model_name in names(model_list)) {
  order <- model_list[[model_name]]
  fit <- arima(gdp_ts, order = order)
  results_list[[model_name]] <- model_params(fit, model_name)
}

# Get all possible coefficients names
all_keys <- unique(unlist(lapply(results_list, names)))
print(all_keys)

# Build the final data frame
results_df <- as.data.frame(matrix(NA, nrow = length(all_keys), ncol = length(results_list)))
rownames(results_df) <- all_keys
colnames(results_df) <- names(results_list)

# Fill the data frame
for (model_name in names(results_list)) {
  for (key in names(results_list[[model_name]])) {
    results_df[key, model_name] <- round(results_list[[model_name]][[key]], 4)
  }
}

#Just setting correct the cols ans rows names

row.names(results_df) <- c("AR1", "Std Error AR1", "Pval AR1", "Intercept", 
                           "Std Error Intercept", "Pval Intercept", "AIC", "BIC", 
                           "AR2", "Std Error AR2", "Pval AR2", "MA1", 
                           "Std Error MA1", "Pval MA1", "MA2", 
                           "Std Error MA2", "Pval MA2")
colnames(results_df) <- c("AR(1)", "AR(2)", "MA(1)", "MA(2)", "ARMA(1,1)", "ARMA(2,1)", 
                          "ARMA(1,2)", "ARMA(2,2)")

results_df <- results_df[c("Intercept", "Std Error Intercept",  "Pval Intercept",  
                           "AR1", "Std Error AR1", "Pval AR1", "AR2", "Std Error AR2", 
                           "Pval AR2", "MA1",  "Std Error MA1",  "Pval MA1", "MA2",  
                           "Std Error MA2", "Pval MA2", "AIC", "BIC"), ]
xtable(results_df, 
       type = "latex", 
       caption = "Model Estimates for Brazilian GDP Growth (AR/MA/ARMA)", 
       digits = 4)


################################################################################
##################################Item c########################################

gdp <- pset_1_gdp$gdp
years <- pset_1_gdp$year
gdp_ts <- ts(gdp, start = min(years), frequency = 1)
models <- list(
  AR1 = arima(gdp_ts, order = c(1, 0, 0)),
  AR2 = arima(gdp_ts, order = c(2, 0, 0)),
  MA1 = arima(gdp_ts, order = c(0, 0, 1)),
  MA2 = arima(gdp_ts, order = c(0, 0, 2)),
  ARMA11 = arima(gdp_ts, order = c(1, 0, 1)),
  ARMA21 = arima(gdp_ts, order = c(2, 0, 1)),
  ARMA12 = arima(gdp_ts, order = c(1, 0, 2)),
  ARMA22 = arima(gdp_ts, order = c(2, 0, 2))
)

model_names <- c("AR(1)", "AR(2)", "MA(1)", "MA(2)",
                 "ARMA(1,1)", "ARMA(2,1)", "ARMA(1,2)", "ARMA(2,2)")




# Forecast for 10 periods ahead
h <- 10
plots <- list()
last_row_indices <- c(7, 8)


# Loop through the models to generate forecasts and plots
for (i in seq_along(models)) {
  fit <- models[[i]]
  model_name <- model_names[i]
  fc <- forecast(fit, h = h)
  
 
  last_hist_year <- tail(time(gdp_ts), 1)
  cutoff_year <- as.integer(last_hist_year) - 29  # last 30 years before the forecast
  
  
  gdp_ts_cut <- window(gdp_ts, start = cutoff_year)
  fitted_cut <- fitted(fit)
  fitted_cut <- tail(fitted_cut, length(gdp_ts_cut)) 
  
  hist_df <- data.frame(
    Year = as.integer(time(gdp_ts_cut)),
    Value = as.numeric(gdp_ts_cut),
    Fitted = as.numeric(fitted_cut)
  )
  
  forecast_df <- data.frame(
    Year = as.integer(time(fc$mean)),
    Forecast = as.numeric(fc$mean),
    Lower = as.numeric(fc$lower[,2]),
    Upper = as.numeric(fc$upper[,2])
  )
  
  # setting total range: last 30 years + 10 years forecast
  x_min <- cutoff_year
  x_max <- max(forecast_df$Year)
  
  x_breaks <- unique(c(seq(min(hist_df$Year), max(forecast_df$Year), by = 10), max(forecast_df$Year)))
  
  
  x_axis_theme <- if (i %in% last_row_indices) {
    list(
      scale_x_continuous(breaks = x_breaks, limits = c(x_min, x_max)),
      theme(
        axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks.x = element_line()
      )
    )
  } else {
    list(
      scale_x_continuous(limits = c(x_min, x_max)),
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    )
  }
  
  forecast_plot <- ggplot() +
    geom_line(data = hist_df, aes(x = Year, y = Value), color = "black") +
    geom_line(data = hist_df, aes(x = Year, y = Fitted), color = "#C5240E", linetype = "dashed") +
    geom_line(data = forecast_df, aes(x = Year, y = Forecast), color = "#559F52") +
    geom_ribbon(data = forecast_df, aes(x = Year, ymin = Lower, ymax = Upper), fill = "#926590", alpha = 0.2) +
    labs(title = paste("Forecast -", model_name),
         x = "Year",
         y = "GDP") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 8, face = "bold"),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8)
    ) +
    x_axis_theme
  
  plots[[i]] <- forecast_plot
}

ggsave("forecast_10y_ahead.png", 
       plot = grid.arrange(grobs = plots, ncol = 2), 
       width = 10, height = 7, dpi = 300)


################################################################################
##################################Item d########################################

model_metrics <- list()
residual_plots <- list()


for (model_name in names(model_list)) {
  order <- model_list[[model_name]]
  
  
  model_fit <- arima(gdp_ts, order = order)
  
  
  #fc <- forecast(model_fit, h = 10, level = 95)
  
  actual <- gdp_ts
  fitted_vals <- fitted(model_fit)
  
  # getting the residuals
  resid <- residuals(model_fit)
  
  # calculating the metrics
  rmse <- sqrt(mean(resid^2, na.rm = TRUE))
  mae <- mean(abs(resid), na.rm = TRUE)
  
  # saving all the metrics
  model_metrics[[model_name]] <- tibble(
    Model = model_name,
    AIC = AIC(model_fit),
    BIC = BIC(model_fit),
    RMSE = rmse,
    MAE = mae,
  )
  # creating the residuals graph
  residual_df <- data.frame(
    Time = time(residuals(model_fit)),
    Residuals = as.numeric(residuals(model_fit))
  )
  
  residual_plot <- ggplot(residual_df, aes(x = Time, y = Residuals)) +
    geom_line(color = "#FF5733") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = paste("Residuals -", model_name),
         x = "Year",
         y = "Residuals") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 8, face = "bold"),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      axis.text.x = element_blank(),  
      axis.ticks.x = element_blank() 
    )
  
  residual_plots[[model_name]] <- residual_plot
}

ggsave("forecast_10y_ahead_residuals.png", 
       plot = grid.arrange(grobs = residual_plots, ncol = 2), 
       width = 10, height = 5, dpi = 300) 

metrics_df <- bind_rows(model_metrics)

gridExtra::grid.arrange(grobs = residual_plots, ncol = 2)

metrics_df <- as.data.frame(metrics_df)
metrics_df <- metrics_df %>%
  arrange(desc(AIC)) 

xtable(metrics_df, type = "latex",
       caption = "Models Analysis",
       digits = 4,
       include.rownames = FALSE)




################################################################################
##################################Question 2####################################

################################################################################
###################################MA(1)########################################

set.seed(123)

# Define sample sizes
sample_sizes <- 30:500

# Simulation parameters
n_sim <- 1000  # Number of simulations per sample size
true_theta <- 0.5
true_c <- 0

# Function to run simulations for a given error distribution
run_simulations <- function(error_dist = "normal") {
  results <- data.frame()
  
  for (n in sample_sizes) {
    c_estimates_normal <- numeric(n_sim)
    theta_estimates_normal <- numeric(n_sim)
    p_values <- numeric(n_sim)
    
    for (i in 1:n_sim) {
      # Generate errors
      if (error_dist == "normal") {
        eps <- rnorm(n, mean = 0, sd = 1)
      } else if (error_dist == "exponential") {
        eps <- rexp(n, rate = 1) - 1  # Center exponential(1) to have mean 0
      }
      
      # Generate MA(1) process
      y <- numeric(n)
      y[1] <- true_c + eps[1]
      for (t in 2:n) {
        y[t] <- true_c + true_theta * eps[t-1] + eps[t]
      }
      
      # Fit MA(1) model
      fit <- try(Arima(y, order = c(0, 0, 1), include.mean = TRUE, method = "ML"), silent = TRUE)
      
      if (!inherits(fit, "try-error")) {
        c_estimates_normal[i] <- coef(fit)["intercept"]
        theta_estimates_normal[i] <- coef(fit)["ma1"]
        
        # Test H0: c = 0
        se <- sqrt(diag(fit$var.coef))["intercept"]
        test_stat <- c_estimates_normal[i]/se
        p_values[i] <- 2 * (1 - pnorm(abs(test_stat)))
      } else {
        c_estimates_normal[i] <- NA
        theta_estimates_normal[i] <- NA
        p_values[i] <- NA
      }
    }
    
    # Calculate summary statistics
    mean_c <- mean(c_estimates_normal, na.rm = TRUE)
    bias_c <- mean_c - true_c
    var_c <- var(c_estimates_normal, na.rm = TRUE)
    mse_c <- bias_c^2 + var_c
    
    # Test size (proportion of rejections when H0 is true)
    test_size <- mean(p_values < 0.05, na.rm = TRUE)
    
    results <- rbind(results, data.frame(
      n = n,
      mean_c = mean_c,
      bias_c = bias_c,
      var_c = var_c,
      mse_c = mse_c,
      test_size = test_size,
      error_dist = error_dist
    ))
  }
  
  return(results)
}

# Run simulations for both error distributions
results_normal <- run_simulations("normal")
results_exp <- run_simulations("exponential")

# Combine results
all_results <- rbind(results_normal, results_exp)

## a. Convergence in probability (bias and MSE plots)
# Bias plot
p_bias <- ggplot(all_results, aes(x = n, y = bias_c, color = error_dist)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Bias of Intercept Estimator",
       x = "Sample Size", y = "Bias",
       color = "Error Distribution") +
  theme_minimal()

# MSE plot
p_mse <- ggplot(all_results, aes(x = n, y = mse_c, color = error_dist)) +
  geom_line() +
  geom_point() +
  labs(title = "MSE of Intercept Estimator",
       x = "Sample Size", y = "MSE",
       color = "Error Distribution") +
  theme_minimal()

## b. Convergence in distribution (QQ plots for largest sample size)
# Get estimates for largest sample size (n=500)
get_large_sample_estimates <- function(error_dist) {
  n_large <- 500
  c_estimates <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # Generate errors
    if (error_dist == "normal") {
      eps <- rnorm(n_large, mean = 0, sd = 1)
    } else if (error_dist == "exponential") {
      eps <- rexp(n_large, rate = 1) - 1
    }
    
    # Generate MA(1) process
    y <- numeric(n_large)
    y[1] <- true_c + eps[1]
    for (t in 2:n_large) {
      y[t] <- true_c + true_theta * eps[t-1] + eps[t]
    }
    
    # Fit MA(1) model
    fit <- try(Arima(y, order = c(0, 0, 1), include.mean = TRUE, method = "ML"), silent = TRUE)
    
    if (!inherits(fit, "try-error")) {
      c_estimates[i] <- coef(fit)["intercept"]
    } else {
      c_estimates[i] <- NA
    }
  }
  
  return(na.omit(c_estimates))
}

c_normal <- get_large_sample_estimates("normal")
c_exp <- get_large_sample_estimates("exponential")

# Standardize the estimates
c_normal_std <- (c_normal - mean(c_normal))/sd(c_normal)
c_exp_std <- (c_exp - mean(c_exp))/sd(c_exp)

# Create QQ plot data
qq_normal <- data.frame(
  theoretical = qnorm(ppoints(length(c_normal_std))),
  sample = sort(c_normal_std),
  dist = "Normal Errors"
)

qq_exp <- data.frame(
  theoretical = qnorm(ppoints(length(c_exp_std))),
  sample = sort(c_exp_std),
  dist = "Exponential Errors"
)

qq_data <- rbind(qq_normal, qq_exp)

# QQ plot
p_qq <- ggplot(qq_data, aes(x = theoretical, y = sample)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  facet_wrap(~dist) +
  labs(title = "QQ Plots of Standardized Intercept Estimates (n=500)",
       x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()

## c. Test size control
p_test_size <- ggplot(all_results, aes(x = n, y = test_size, color = error_dist)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  labs(title = "Test Size for H0: c = 0 (α = 0.05)",
       x = "Sample Size", y = "Rejection Rate",
       color = "Error Distribution") +
  ylim(0, 0.15) +
  theme_minimal()

# Display all plots
print(p_bias)
print(p_mse)
print(p_qq)
print(p_test_size)





################################################################################
############################Question 2##########################################



################Rascunhos dont run anything bellow##############################
convergence_in_dist <- function(results, sample_sizes_vec){
  
  
  if (length(results) != length(sample_sizes_vec)){
    estimated_samples <- c()
    for (i in 1:length(results)){
      estimated_samples[i] <- results[[i]]$Normal$sample_size[1]
    }
    sample_sizes_vec <- estimated_samples
  }
  
  # Create DataFrames to store the results
  ## Normal error
  temp_normal <- data.frame(
    "sample_size" = rep(NA, 6 * 101),
    "Fy" = rep(NA, 6 * 101),
    "Qy" = rep(NA, 6 * 101)
  )
  
  ## Exponential error
  temp_exp <- data.frame(
    "sample_size" = rep(NA, 6 * 101),
    "Fy" = rep(NA, 6 * 101),
    "Qy" = rep(NA, 6 * 101)
  )
  
  # Find the relevant sample size indexes
  i_vec <- which(sample_sizes_vec %in% c(30, 40, 50, 100, 250, 500))
  
  # Loop over the sample sizes (for the normal)
  for (i in i_vec) {
    # Index within i_vec
    j <- which(i == i_vec)
    
    # Write the sample size
    temp_normal$sample_size[(1 + (j - 1) * 101):(j * 101)] <- 
      results[[i]]$Normal$sample_size[1]
    temp_exp$sample_size[(1 + (j - 1) * 101):(j * 101)] <- 
      results[[i]]$Exponential$sample_size[1]
    
    # Write down the probabilities
    temp_normal$Fy[(1 + (j - 1) * 101):(j * 101)] <- seq(0, 1, 0.01)
    temp_exp$Fy[(1 + (j - 1) * 101):(j * 101)] <- seq(0, 1, 0.01)
    
    # Write down the quantiles
    # Para Normal
    coefs_normal <- results[[i]]$Normal$normalized_coef
    if (!is.null(coefs_normal) && any(!is.na(coefs_normal))) {
      temp_normal$Qy[(1 + (j - 1) * 101):(j * 101)] <- quantile(
        coefs_normal, probs = seq(0, 1, 0.01), na.rm = TRUE
      )
    } else {
      temp_normal$Qy[(1 + (j - 1) * 101):(j * 101)] <- rep(NA, 101)
    }
    
    # Para Exponential
    coefs_exp <- results[[i]]$Exponential$normalized_coef
    if (!is.null(coefs_exp) && any(!is.na(coefs_exp))) {
      temp_exp$Qy[(1 + (j - 1) * 101):(j * 101)] <- quantile(
        coefs_exp, probs = seq(0, 1, 0.01), na.rm = TRUE
      )
    } else {
      temp_exp$Qy[(1 + (j - 1) * 101):(j * 101)] <- rep(NA, 101)
    }
  }
  
  # Remove os valores NA de Qy e Fy
  temp_normal <- temp_normal %>% na.omit()
  temp_exp <- temp_exp %>% na.omit()
  
  # Write sample size as a factor to enforce the ordering
  temp_normal$sample_size <- factor(temp_normal$sample_size)
  temp_exp$sample_size <- factor(temp_exp$sample_size)
  
  # Create ggplots
  ## Normal
  gg_normal <- ggplot(temp_normal, aes(x = Qy, y = Fy)) +
    theme_bw(base_size = 25) +
    theme(plot.margin = unit(c(5, 7, 2, 2), "mm")) +
    xlab("Normalized Coefficient") + ylab("CDF") +
    ggtitle("Convergence in Distribution - Normal Errors") + 
    geom_line(aes(colour = sample_size), size = 1) +
    scale_colour_manual(name = "", values = met.brewer("Demuth")) +
    stat_function(fun = pnorm, size = .5, linetype = "dashed", colour = '#171110')  +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, title = "Sample Size")) +
    theme(legend.position = "bottom",
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.margin = margin(c(1,5,5,5)),
          plot.title = element_text(size = 16, face="bold"))
  
  ## Exponential
  gg_exp <- ggplot(temp_exp, aes(x = Qy, y = Fy)) +
    theme_bw(base_size = 25) +
    theme(plot.margin = unit(c(5, 7, 2, 2), "mm")) +
    xlab("Normalized Coefficient") + ylab("CDF") +
    ggtitle("Convergence in Distribution - Exponential Errors") + 
    geom_line(aes(colour = sample_size), size = 1) +
    scale_colour_manual(name = "", values = met.brewer("Moreau")) +
    stat_function(fun = pnorm, size = .5, linetype = "dashed", colour = '#171110')  +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, title = "Sample Size")) +
    theme(legend.position = "bottom",
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.margin = margin(c(1,5,5,5)),
          plot.title = element_text(size = 16, face="bold"))
  
  # With facets
  ## Normal
  gg_normal_facet <- ggplot(temp_normal, aes(x = Qy, y = Fy)) +
    theme_bw(base_size = 25) +
    theme(plot.margin = unit(c(5, 7, 2, 2), "mm")) +
    xlab("Normalized Coefficient") + ylab("CDF") +
    ggtitle("Convergence in Distribution - Normal Errors") + 
    geom_line(aes(colour = sample_size), size = 1) +
    scale_colour_manual(name = "", values = met.brewer("Demuth")) +
    stat_function(fun = pnorm, size = .5, linetype = "dashed", colour = '#171110')  +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, title = "Sample Size")) +
    theme(legend.position = "bottom",
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.margin = margin(c(1,5,5,5)),
          plot.title = element_text(size = 16, face="bold")) + 
    facet_wrap(vars(sample_size), labeller = label_parsed, ncol = 2, nrow = 3)
  
  ## Exponential
  gg_exp_facet <- ggplot(temp_exp, aes(x = Qy, y = Fy)) +
    theme_bw(base_size = 25) +
    theme(plot.margin = unit(c(5, 7, 2, 2), "mm")) +
    xlab("Normalized Coefficient") + ylab("CDF") +
    ggtitle("Convergence in Distribution - Exponential Errors") + 
    geom_line(aes(colour = sample_size), size = 1) +
    scale_colour_manual(name = "", values = met.brewer("Demuth")) +
    stat_function(fun = pnorm, size = .5, linetype = "dashed", colour = '#171110')  +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, title = "Sample Size")) +
    theme(legend.position = "bottom",
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.margin = margin(c(1,5,5,5)),
          plot.title = element_text(size = 16, face="bold")) + 
    facet_wrap(vars(sample_size), labeller = label_parsed, ncol = 2, nrow = 3)
  
  # Returning
  output_list <- list("gg_normal" = gg_normal,
                      "gg_exp" = gg_exp,
                      "gg_normal_facet" = gg_normal_facet,
                      "gg_exp_facet" = gg_exp_facet)
  return(output_list)
}
monte_carlo <- function(c, theta, phi, 
                        n_sim, sample_sizes_vec, 
                        delta_vec, order_arma, coefficient){
  
  # Determinando o modelo
  if (identical(order_arma, c(1, 0, 0))){ # AR
    model_list <- list(ar = phi)
    index_coefficient <- ifelse(coefficient == "intercept", 2, 1)
    parameter <- ifelse(coefficient == "intercept", c, phi)
    
  } else if (identical(order_arma, c(0, 0, 1))){ # MA
    model_list <- list(ma = theta)
    index_coefficient <- ifelse(coefficient == "intercept", 2, 1)
    parameter <- ifelse(coefficient == "intercept", c, theta)
    
  } else { # ARMA
    model_list <- list(ar = phi, ma = theta)
    index_coefficient <- ifelse(coefficient == "intercept", 3, 
                                ifelse(coefficient == "ma1", 2, 1))
    parameter <- ifelse(coefficient == "intercept", c, 
                        ifelse(coefficient == "ma1", theta, phi))
  }
  
  # Lista de resultados
  results <- list()
  
  for (sample in sample_sizes_vec) {
    
    df_mc_normal_results <- data.frame(
      "sample_size" = rep(sample, n_sim),
      "parameter_estimate" = NA,
      "bias" = NA,
      "normalized_coef" = NA,
      "reject" = NA
    )
    
    df_mc_exp_results <- data.frame(
      "sample_size" = rep(sample, n_sim),
      "parameter_estimate" = NA,
      "bias" = NA,
      "normalized_coef" = NA,
      "reject" = NA
    )
    
    for (sim in 1:n_sim) {
      # Simulação com erro normal
      data_normal <- arima.sim(model = model_list, n = sample, rand.gen = rnorm)
      
      # Simulação com erro exponencial
      data_exp <- arima.sim(model = model_list, n = sample, rand.gen = rexp)
      
      # Estimação dos modelos
      model_normal <- arima(x = data_normal, order = order_arma, include.mean = TRUE)
      model_exp <- arima(x = data_exp, order = order_arma, include.mean = TRUE)
      
      # Estimativa do parâmetro
      df_mc_normal_results$parameter_estimate[sim] <- model_normal$coef[index_coefficient]
      df_mc_exp_results$parameter_estimate[sim] <- model_exp$coef[index_coefficient]
      
      # Viés
      df_mc_normal_results$bias[sim] <- model_normal$coef[index_coefficient] - parameter
      df_mc_exp_results$bias[sim] <- model_exp$coef[index_coefficient] - parameter
      
      # Estatística normalizada
      df_mc_normal_results$normalized_coef[sim] <- 
        (model_normal$coef[index_coefficient] - parameter) /
        sqrt(model_normal$var.coef[index_coefficient, index_coefficient])
      
      df_mc_exp_results$normalized_coef[sim] <- 
        (model_exp$coef[index_coefficient] - parameter) /
        sqrt(model_exp$var.coef[index_coefficient, index_coefficient])
      
      # Teste de hipótese (nível de 5%)
      df_mc_normal_results$reject[sim] <- 
        as.numeric(abs(df_mc_normal_results$normalized_coef[sim]) >= qnorm(0.975))
      
      df_mc_exp_results$reject[sim] <- 
        as.numeric(abs(df_mc_exp_results$normalized_coef[sim]) >= qnorm(0.975))
    }
    
    # Salvar os resultados em lista
    results[[as.character(sample)]] <- list(
      Normal = df_mc_normal_results,
      Exponential = df_mc_exp_results
    )
  }
  
  return(results)
}
##########################You should run form this on###########################



monte_carlo <- function(c, theta, phi, 
                        n_sim, sample_sizes_vec, 
                        delta_vec, n_cores,
                        order_arma, coefficient){
  
  # Determining model
  if (identical(order_arma, c(1, 0, 0))){ # AR
    model_list <- list(ar = phi)
    index_coefficient <- ifelse(coefficient == "intercept", 2, 1)
    parameter <- ifelse(coefficient == "intercept", c, phi)
    
  } else if (identical(order_arma, c(0, 0, 1))){ # MA
    model_list <- list(ma = theta)
    index_coefficient <- ifelse(coefficient == "intercept", 2, 1)
    parameter <- ifelse(coefficient == "intercept", c, theta)
    
  } else { # ARMA
    model_list <- list(ar = phi, ma = theta)
    index_coefficient <- ifelse(coefficient == "intercept", 3, 
                                ifelse(coefficient == "ma1", 2, 1))
    parameter <- ifelse(coefficient == "intercept", c, 
                        ifelse(coefficient == "ma1", theta, phi))
  }
  
  # Setup parallel backend to use all cores (yes, I have 2)
  n_cores <- 16
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Run a parallel loop over sample sizes and keeping time
  tic()
  results <- foreach(
    sample = sample_sizes_vec, .inorder = TRUE, 
    .errorhandling = "remove", .verbose = FALSE
  ) %dopar% {
    # DataFrame to store the result of each MC simulation
    # Each column is a result; each line is a simulation
    df_mc_normal_results <- data.frame(
      "sample_size" = rep(sample, n_sim),
      "parameter_estimate" = rep(NA, n_sim),
      "bias" = rep(NA, n_sim),
      "normalized_coef" = rep(NA, n_sim),
      "reject" = rep(NA, n_sim)
    )
    
    df_mc_exp_results <- data.frame(
      "sample_size" = rep(sample, n_sim),
      "parameter_estimate" = rep(NA, n_sim),
      "bias" = rep(NA, n_sim),
      "normalized_coef" = rep(NA, n_sim),
      "reject" = rep(NA, n_sim)
    )
    
    # Looping over MC repetitions and simulating MA(1) processes
    for (sim in 1:n_sim) {
      # Errors following a Normal distribution (they are iid by assumption)
      data_normal <- arima.sim(model = model_list, 
                               n = sample, rand.gen = rnorm)
      
      # Errors following an Exponential distribution (they are iid by assumption)
      data_exp <- arima.sim(model = model_list, 
                            n = sample, rand.gen = rexp)
      
      # Estimating models
      model_normal <- arima(x = data_normal, order = order_arma, 
                            include.mean = T)
      model_exp <- arima(x = data_exp, order = order_arma, 
                         include.mean = T)
      
      # Storing results
      ## Parameter estimate
      df_mc_normal_results$parameter_estimate[sim] <- 
        model_normal$coef[index_coefficient]
      
      df_mc_exp_results$parameter_estimate[sim] <- 
        model_exp$coef[index_coefficient]
      
      ## Bias
      df_mc_normal_results$bias[sim] <- 
        model_normal$coef[index_coefficient] - parameter
      
      df_mc_exp_results$bias[sim] <- 
        model_exp$coef[index_coefficient] - parameter
      
      ## Normalized coefficients
      df_mc_normal_results$normalized_coef[sim] <- 
        (model_normal$coef[index_coefficient] - parameter) /
        sqrt(model_normal$var.coef[index_coefficient, index_coefficient])
      
      df_mc_exp_results$normalized_coef[sim] <- 
        (model_exp$coef[index_coefficient] - parameter) /
        sqrt(model_exp$var.coef[index_coefficient, index_coefficient])
      
      ## Test decision
      df_mc_normal_results$reject[sim] <- as.numeric(
        abs(df_mc_normal_results$normalized_coef[sim]) >= qnorm(0.975)
      )
      
      df_mc_exp_results$reject[sim] <- as.numeric(
        abs(df_mc_exp_results$normalized_coef[sim]) >= qnorm(0.975)
      )
    }
    
    # Return the results
    ## List
    output <- list("Normal" = df_mc_normal_results,
                   "Exponential" = df_mc_exp_results)
    ## Output
    return(output)
  }
  
  # Stop parallel backend
  stopCluster(cl)
  toc()
  
  # Returning
  return(results)
}

convergence_in_prob <- function(results, delta_vec, sample_sizes_vec){
  
  # In ARMA, some sample sizes are run with error (message says something
  # about non-stationary AR with KSS test), so `results` is of a smaller length
  # This fixes it: I take only the sample sizes thar throw no error
  if (length(results) != length(sample_sizes_vec)){
    estimated_samples <- c()
    for (i in 1:length(results)){
      estimated_samples[i] <- results[[i]]$Normal$sample_size[1]
    }
    sample_sizes_vec <- estimated_samples
  }
  
  # DataFrames that stores the results
  ## Normal
  df_conv_prob_normal <- data.frame(
    "sample_size" = sample_sizes_vec,
    "delta1" = rep(NA, length(sample_sizes_vec)),
    "delta2" = rep(NA, length(sample_sizes_vec)),
    "delta3" = rep(NA, length(sample_sizes_vec)),
    "delta4" = rep(NA, length(sample_sizes_vec))
  )
  
  ## Exponential
  df_conv_prob_exp <- df_conv_prob_normal
  
  # Looping over samples
  for (sample_size in 1:length(sample_sizes_vec)){
    # For each one of them, loop across the value of delta
    for (delta in delta_vec){
      # Compare results against delta
      ## Is it significant?
      significant_normal <- abs(results[[sample_size]]$Normal$bias) > delta
      significant_exp <- abs(results[[sample_size]]$Exponential$bias) > delta
      
      ## Storing (just copied, don't know what is going on)
      df_conv_prob_normal[sample_size, which(delta == delta_vec) + 1] <-
        mean(significant_normal)
      df_conv_prob_exp[sample_size, which(delta == delta_vec) + 1] <-
        mean(significant_exp)
      
    }
  }
  
  # Create a plot with the results
  ## Normal
  gg_normal <- ggplot(data = df_conv_prob_normal, aes(x = sample_size)) +
    theme_bw(base_size = 25) +
    theme(plot.margin = unit(c(5, 7, 2, 2), "mm")) +
    xlab("Sample Size") + ylab("Probability") +
    ggtitle("Convergence in Probability - Normal Errors") + 
    geom_line(aes(y = delta1, color = "d = 0.25"), size = 1.5) +
    geom_line(aes(y = delta2, color = "d = 0.20"), size = 1.5) +
    geom_line(aes(y = delta3, color = "d = 0.15"), size = 1.5) +
    geom_line(aes(y = delta4, color = "d = 0.10"), size = 1.5) +
    scale_colour_manual(name = "", values = met.brewer("Demuth")
    ) +
    theme(
      legend.title = element_blank(),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      plot.title = element_text(size = 16, face="bold"),
      legend.position = "bottom"
    )
  
  ## Exponential
  gg_exp <- ggplot(data = df_conv_prob_exp, aes(x = sample_size)) +
    theme_bw(base_size = 25) +
    theme(plot.margin = unit(c(5, 7, 2, 2), "mm")) +
    xlab("Sample Size") + ylab("Probability") +
    ggtitle("Convergence in Probability - Exponential Errors") + 
    geom_line(aes(y = delta1, color = "d = 0.25"), size = 1.5) +
    geom_line(aes(y = delta2, color = "d = 0.20"), size = 1.5) +
    geom_line(aes(y = delta3, color = "d = 0.15"), size = 1.5) +
    geom_line(aes(y = delta4, color = "d = 0.10"), size = 1.5) +
    scale_colour_manual(values = met.brewer("Demuth")[c(1:4)]) +
    theme(
      legend.title = element_blank(),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      plot.title = element_text(size = 16, face="bold"),
      legend.position = "bottom"
    )
  
  # Returning
  output_list <- list("data_normal" = df_conv_prob_normal,
                      "data_exp" = df_conv_prob_exp,
                      "gg_normal" = gg_normal,
                      "gg_exp" = gg_exp)
  return(output_list)
}




convergence_in_dist <- function(results, sample_sizes_vec){
  
  if (length(results) != length(sample_sizes_vec)) {
    estimated_samples <- c()
    for (i in 1:length(results)) {
      estimated_samples[i] <- results[[i]]$Normal$sample_size[1]
    }
    sample_sizes_vec <- estimated_samples
  }
  
  # Índices dos tamanhos de amostra desejados
  i_vec <- which(sample_sizes_vec %in% c(30, 40, 50, 100, 250, 500))
  
  # Número real de tamanhos de amostra válidos
  n_samples <- length(i_vec)
  
  # Criação de dataframes com o tamanho certo
  temp_normal <- data.frame(
    "sample_size" = rep(NA, n_samples * 101),
    "Fy" = rep(NA, n_samples * 101),
    "Qy" = rep(NA, n_samples * 101)
  )
  
  temp_exp <- data.frame(
    "sample_size" = rep(NA, n_samples * 101),
    "Fy" = rep(NA, n_samples * 101),
    "Qy" = rep(NA, n_samples * 101)
  )
  
  for (j in seq_along(i_vec)) {
    i <- i_vec[j]
    
    # Preenche sample_size
    temp_normal$sample_size[(1 + (j - 1) * 101):(j * 101)] <- results[[i]]$Normal$sample_size[1]
    temp_exp$sample_size[(1 + (j - 1) * 101):(j * 101)] <- results[[i]]$Exponential$sample_size[1]
    
    # Probabilidades
    temp_normal$Fy[(1 + (j - 1) * 101):(j * 101)] <- seq(0, 1, 0.01)
    temp_exp$Fy[(1 + (j - 1) * 101):(j * 101)] <- seq(0, 1, 0.01)
    
    # Quantis Normal
    coefs_normal <- results[[i]]$Normal$normalized_coef
    if (!is.null(coefs_normal) && any(!is.na(coefs_normal))) {
      temp_normal$Qy[(1 + (j - 1) * 101):(j * 101)] <- quantile(
        coefs_normal, probs = seq(0, 1, 0.01), na.rm = TRUE
      )
    } else {
      temp_normal$Qy[(1 + (j - 1) * 101):(j * 101)] <- rep(NA, 101)
    }
    
    # Quantis Exponential
    coefs_exp <- results[[i]]$Exponential$normalized_coef
    if (!is.null(coefs_exp) && any(!is.na(coefs_exp))) {
      temp_exp$Qy[(1 + (j - 1) * 101):(j * 101)] <- quantile(
        coefs_exp, probs = seq(0, 1, 0.01), na.rm = TRUE
      )
    } else {
      temp_exp$Qy[(1 + (j - 1) * 101):(j * 101)] <- rep(NA, 101)
    }
  }
  
  # Remove NAs
  temp_normal <- temp_normal %>% na.omit()
  temp_exp <- temp_exp %>% na.omit()
  
  # Fatores
  temp_normal$sample_size <- factor(temp_normal$sample_size)
  temp_exp$sample_size <- factor(temp_exp$sample_size)
  
  
  # Create ggplots
  ## Normal
  gg_normal <- ggplot(temp_normal, aes(x = Qy, y = Fy)) +
    theme_bw(base_size = 25) +
    theme(plot.margin = unit(c(5, 7, 2, 2), "mm")) +
    xlab("Normalized Coefficient") + ylab("CDF") +
    ggtitle("Convergence in Distribution - Normal Errors") + 
    geom_line(aes(colour = sample_size), size = 1) +
    scale_colour_manual(name = "", values = met.brewer("Demuth")) +
    stat_function(fun = pnorm, size = .5, linetype = "dashed", colour = '#171110')  +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, title = "Sample Size")) +
    theme(legend.position = "bottom",
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.margin = margin(c(1,5,5,5)),
          plot.title = element_text(size = 16, face="bold"))
  
  ## Exponential
  gg_exp <- ggplot(temp_exp, aes(x = Qy, y = Fy)) +
    theme_bw(base_size = 25) +
    theme(plot.margin = unit(c(5, 7, 2, 2), "mm")) +
    xlab("Normalized Coefficient") + ylab("CDF") +
    ggtitle("Convergence in Distribution - Exponential Errors") + 
    geom_line(aes(colour = sample_size), size = 1) +
    scale_colour_manual(name = "", values = met.brewer("Moreau")) +
    stat_function(fun = pnorm, size = .5, linetype = "dashed", colour = '#171110')  +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, title = "Sample Size")) +
    theme(legend.position = "bottom",
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.margin = margin(c(1,5,5,5)),
          plot.title = element_text(size = 16, face="bold"))
  
  # With facets
  ## Normal
  gg_normal_facet <- ggplot(temp_normal, aes(x = Qy, y = Fy)) +
    theme_bw(base_size = 25) +
    theme(plot.margin = unit(c(5, 7, 2, 2), "mm")) +
    xlab("Normalized Coefficient") + ylab("CDF") +
    ggtitle("Convergence in Distribution - Normal Errors") + 
    geom_line(aes(colour = sample_size), size = 1) +
    scale_colour_manual(name = "", values = met.brewer("Demuth")) +
    stat_function(fun = pnorm, size = .5, linetype = "dashed", colour = '#171110')  +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, title = "Sample Size")) +
    theme(legend.position = "bottom",
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.margin = margin(c(1,5,5,5)),
          plot.title = element_text(size = 16, face="bold")) + 
    facet_wrap(vars(sample_size), labeller = label_parsed, ncol = 2, nrow = 3)
  
  ## Exponential
  gg_exp_facet <- ggplot(temp_exp, aes(x = Qy, y = Fy)) +
    theme_bw(base_size = 25) +
    theme(plot.margin = unit(c(5, 7, 2, 2), "mm")) +
    xlab("Normalized Coefficient") + ylab("CDF") +
    ggtitle("Convergence in Distribution - Exponential Errors") + 
    geom_line(aes(colour = sample_size), size = 1) +
    scale_colour_manual(name = "", values = met.brewer("Demuth")) +
    stat_function(fun = pnorm, size = .5, linetype = "dashed", colour = '#171110')  +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, title = "Sample Size")) +
    theme(legend.position = "bottom",
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.margin = margin(c(1,5,5,5)),
          plot.title = element_text(size = 16, face="bold")) + 
    facet_wrap(vars(sample_size), labeller = label_parsed, ncol = 2, nrow = 3)
  
  # Returning
  output_list <- list("gg_normal" = gg_normal,
                      "gg_exp" = gg_exp,
                      "gg_normal_facet" = gg_normal_facet,
                      "gg_exp_facet" = gg_exp_facet)
  return(output_list)
}

test_size <- function(results, sample_sizes_vec){
  # # In ARMA, some sample sizes are run with error (message says something
  # about non-stationary AR with KSS test), so `results` is of a smaller length
  # This fixes it: I take only the sample sizes thar throw no error
  if (length(results) != length(sample_sizes_vec)){
    estimated_samples <- c()
    for (i in 1:length(results)){
      estimated_samples[i] <- results[[i]]$Normal$sample_size[1]
    }
    sample_sizes_vec <- estimated_samples
  }
  
  # Create a matrix to store the results
  ## Normal
  rej_rate_normal <- data.frame(
    "sample_size" = sample_sizes_vec,
    "rej" = rep(NA, length(sample_sizes_vec)))
  
  ## Exponential
  rej_rate_exp <- data.frame(
    "sample_size" = sample_sizes_vec,
    "rej" = rep(NA, length(sample_sizes_vec)))
  
  # Loop over sample size
  for (t in 1:length(sample_sizes_vec)) {
    # Compute the test size
    rej_rate_normal[t, 2] <- mean(results[[t]]$Normal$reject, na.rm = TRUE)
    rej_rate_exp[t, 2] <- mean(results[[t]]$Exponential$reject, na.rm = TRUE)
    
  }
  
  # Create a plot with the results
  ## Normal
  gg_normal <- ggplot(data = rej_rate_normal, aes(x = sample_size)) +
    xlab("Sample Size") + ylab("Rejection Rate") +
    ggtitle("Test Size - Normal Errors") + 
    # scale_y_continuous(limits = c(0, 0.5)) +
    geom_line(aes(y = rej), size = 1.5, colour = "#558F07") +
    geom_hline(yintercept = 0.05, size = .5, 
               linetype = "dashed", color = "red") + 
    theme_bw() +
    theme(axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(size = 16, face="bold"))
  
  
  ## Exponential
  gg_exp <- ggplot(data = rej_rate_exp, aes(x = sample_size)) +
    theme_bw(base_size = 25) +
    theme(plot.margin = unit(c(5, 7, 2, 2), "mm")) +
    xlab("Sample Size") + ylab("Rejection Rate") +
    ggtitle("Test Size - Exponential Errors") + 
    # scale_y_continuous(limits = c(0, 0.5)) +
    geom_line(aes(y = rej), size = 1.5, colour = "#558F07") +
    geom_hline(yintercept = 0.05, size = .5, 
               linetype = "dashed", color = "red") + 
    theme_bw() + 
    theme(axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(size = 16, face="bold"))
  
  # Returning
  output_list <- list("gg_normal" = gg_normal,
                      "gg_exp" = gg_exp)
  return(output_list)
}



################################################################################
############################MA(1) Intercept#####################################

results_ma1_intercept <- monte_carlo(
  c = 0, 
  theta = 0.5, 
  phi = 0.3,
  n_sim = 100,
  sample_sizes_vec = c(30, 35, 50, 100, 250, 500),
  delta_vec = c(0.25, 0.20, 0.15, 0.10),
  order_arma = c(0, 0, 1),
  coefficient = "intercept"
)


results_ma1_intercept <- monte_carlo(
  c = 0, 
  theta = 0.5, 
  phi = 0.3,
  n_sim = 300,
  sample_sizes_vec = c(
    seq(10, 20, by = 1),
    seq(20, 50, by = 5),
    seq(50, 60, by = 1),
    seq(60, 90, by = 5),
    seq(90, 150, by = 10),
    seq(150, 300, by = 25),
    seq(300, 500, by = 50)
  ),
  delta_vec = c(0.25, 0.20, 0.15, 0.10),  # esse parâmetro não é usado na função
  order_arma = c(0, 0, 1),
  coefficient = "intercept"
)
################################################################################
#######################Convergence in probability###############################
conv_prob_ma1_intercept <- convergence_in_prob(results_ma1_intercept,
                                               delta_vec = c(0.25, 0.20, 0.15, 0.10),
                                               sample_sizes_vec = c(10, 20, 50, 60, 90, 150, 300, 400, 500))

conv_prob_ma1_intercept$gg_normal
conv_prob_ma1_intercept$gg_exp

ggsave("gg_cov_prob_normal_plot_ma1_intercept.png", plot = conv_prob_ma1_intercept$gg_normal, width = 10, height = 8, dpi = 300)
ggsave("gg_cov_prob_exp_plot_ma1_intercept.png", plot = conv_prob_ma1_intercept$gg_exp, width = 10, height = 8, dpi = 300)


################################################################################
#######################Convergence in Distribution##############################
conv_dist_ma_1_intercept <- convergence_in_dist(results_ma1_intercept,
                                                sample_sizes_vec = c(30, 35, 50, 100, 250, 500))

conv_dist_ma_1_intercept$gg_exp
conv_dist_ma_1_intercept$gg_normal
conv_dist_ma_1_intercept$gg_normal_facet
conv_dist_ma_1_intercept$gg_exp_facet

ggsave("gg_cov_dist_expt_plot_ma1_intercept.png", plot = conv_dist_ma_1_intercept$gg_exp, width = 10, height = 8, dpi = 300)
ggsave("gg_cov_dist_normal_plot_ma1_intercept.png", plot = conv_dist_ma_1_intercept$gg_normal, width = 10, height = 8, dpi = 300)


################################################################################
#######################Test Size###############################################

test_size_ma1_intercept <- test_size(results_ma1_intercept,
                                     sample_sizes_vec = c(30, 35, 50, 100, 250, 500))

test_size_ma1_intercept$gg_normal
test_size_ma1_intercept$gg_exp
ggsave("gg_teste_size_normal_plot_ma1_intercept.png", plot = test_size_ma1_intercept$gg_normal, width = 10, height = 8, dpi = 300)
ggsave("gg_teste_size_exp_plot_ma1_intercept.png", plot = test_size_ma1_intercept$gg_exp, width = 10, height = 8, dpi = 300)











################################################################################
############################AR(1) Intercept#####################################


results_ar1_intercept <- monte_carlo(
  c = 0, 
  theta = 0.5, 
  phi = 0.3,
  n_sim = 100,
  sample_sizes_vec = c(
    seq(10, 20, by = 1),
    seq(20, 50, by = 5),
    seq(50, 60, by = 1),
    seq(60, 90, by = 5),
    seq(90, 150, by = 10),
    seq(150, 300, by = 25),
    seq(300, 500, by = 50)
  ),
  delta_vec = c(0.25, 0.20, 0.15, 0.10),  # esse parâmetro não é usado na função
  order_arma = c(1, 0, 0),  # AR(1)
  coefficient = "intercept"
)


################################################################################
#######################Convergence in probability###############################
conv_prob_ar1_intercept <- convergence_in_prob(results_ar1_intercept,
                                               delta_vec = c(0.25, 0.20, 0.15, 0.10),
                                               sample_sizes_vec = c(10, 20, 50, 60, 90, 150, 300, 400, 500))

conv_prob_ar1_intercept$gg_normal
conv_prob_ar1_intercept$gg_exp

ggsave("gg_cov_prob_normal_plot_ar1_intercept.png", plot = conv_prob_ar1_intercept$gg_normal, width = 10, height = 8, dpi = 300)
ggsave("gg_cov_prob_exp_plot_ar1_intercept.png", plot = conv_prob_ar1_intercept$gg_exp, width = 10, height = 8, dpi = 300)


################################################################################
#######################Convergence in Distribution##############################
conv_dist_ar1_intercept <- convergence_in_dist(results_ar1_intercept,
                                               sample_sizes_vec = c(30, 35, 50, 100, 250, 500))

conv_dist_ar1_intercept$gg_exp
conv_dist_ar1_intercept$gg_normal
conv_dist_ar1_intercept$gg_normal_facet
conv_dist_ar1_intercept$gg_exp_facet

ggsave("gg_cov_dist_exp_plot_ar1_intercept.png", plot = conv_dist_ar1_intercept$gg_exp, width = 10, height = 8, dpi = 300)
ggsave("gg_cov_dist_normal_plot_ar1_intercept.png", plot = conv_dist_ar1_intercept$gg_normal, width = 10, height = 8, dpi = 300)


################################################################################
#######################Test Size###############################################

test_size_ar1_intercept <- test_size(results_ar1_intercept,
                                     sample_sizes_vec = c(30, 35, 50, 100, 250, 500))

test_size_ar1_intercept$gg_normal
test_size_ar1_intercept$gg_exp

ggsave("gg_test_size_normal_plot_ar1_intercept.png", plot = test_size_ar1_intercept$gg_normal, width = 10, height = 8, dpi = 300)
ggsave("gg_test_size_exp_plot_ar1_intercept.png", plot = test_size_ar1_intercept$gg_exp, width = 10, height = 8, dpi = 300)








################################################################################



############################AR(1) Coefficient#####################################


results_ar1_coefficient <- monte_carlo(
  c = 0, 
  theta = 0.5, 
  phi = 0.3,
  n_sim = 100,
  sample_sizes_vec = c(
    seq(10, 20, by = 1),
    seq(20, 50, by = 5),
    seq(50, 60, by = 1),
    seq(60, 90, by = 5),
    seq(90, 150, by = 10),
    seq(150, 300, by = 25),
    seq(300, 500, by = 50)
  ),
  delta_vec = c(0.25, 0.20, 0.15, 0.10),  # esse parâmetro não é usado na função
  order_arma = c(1, 0, 0),  # AR(1)
  coefficient = "ar1"
)


################################################################################
####################### Convergence in Probability #############################
conv_prob_ar1_coefficient <- convergence_in_prob(results_ar1_coefficient,
                                                 delta_vec = c(0.25, 0.20, 0.15, 0.10),
                                                 sample_sizes_vec = c(10, 20, 50, 60, 90, 150, 300, 400, 500))

conv_prob_ar1_coefficient$gg_normal
conv_prob_ar1_coefficient$gg_exp

ggsave("gg_cov_prob_normal_plot_ar1_coefficient.png", plot = conv_prob_ar1_coefficient$gg_normal, width = 10, height = 8, dpi = 300)
ggsave("gg_cov_prob_exp_plot_ar1_coefficient.png", plot = conv_prob_ar1_coefficient$gg_exp, width = 10, height = 8, dpi = 300)

################################################################################
####################### Convergence in Distribution ############################
conv_dist_ar1_coefficient <- convergence_in_dist(results_ar1_coefficient,
                                                 sample_sizes_vec = c(30, 35, 50, 100, 250, 500))

conv_dist_ar1_coefficient$gg_exp
conv_dist_ar1_coefficient$gg_normal
conv_dist_ar1_coefficient$gg_normal_facet
conv_dist_ar1_coefficient$gg_exp_facet

ggsave("gg_cov_dist_exp_plot_ar1_coefficient.png", plot = conv_dist_ar1_coefficient$gg_exp, width = 10, height = 8, dpi = 300)
ggsave("gg_cov_dist_normal_plot_ar1_coefficient.png", plot = conv_dist_ar1_coefficient$gg_normal, width = 10, height = 8, dpi = 300)

################################################################################
####################### Test Size ##############################################
test_size_ar1_coefficient <- test_size(results_ar1_coefficient,
                                       sample_sizes_vec = c(30, 35, 50, 100, 250, 500))

test_size_ar1_coefficient$gg_normal
test_size_ar1_coefficient$gg_exp

ggsave("gg_test_size_normal_plot_ar1_coefficient.png", plot = test_size_ar1_coefficient$gg_normal, width = 10, height = 8, dpi = 300)
ggsave("gg_test_size_exp_plot_ar1_coefficient.png", plot = test_size_ar1_coefficient$gg_exp, width = 10, height = 8, dpi = 300)




################################################################################





############################ARMA(1, 1) Intercept#####################################


results_arma11_intercept <- monte_carlo(
  c = 0, 
  theta = 0.5, 
  phi = 0.3,
  n_sim = 100,
  sample_sizes_vec = c(
    seq(10, 20, by = 1),
    seq(20, 50, by = 5),
    seq(50, 60, by = 1),
    seq(60, 90, by = 5),
    seq(90, 150, by = 10),
    seq(150, 300, by = 25),
    seq(300, 500, by = 50)
  ),
  delta_vec = c(0.25, 0.20, 0.15, 0.10),  # esse parâmetro não é usado na função
  order_arma = c(1, 0, 1),  # ARMA(1)
  coefficient = "intercept"
)


################################################################################
####################### Convergence in Probability #############################
conv_prob_arma11_intercept <- convergence_in_prob(results_arma11_intercept,
                                                  delta_vec = c(0.25, 0.20, 0.15, 0.10),
                                                  sample_sizes_vec = c(10, 20, 50, 60, 90, 150, 300, 400, 500))

conv_prob_arma11_intercept$gg_normal
conv_prob_arma11_intercept$gg_exp

ggsave("gg_cov_prob_normal_plot_arma11_intercept.png", plot = conv_prob_arma11_intercept$gg_normal, width = 10, height = 8, dpi = 300)
ggsave("gg_cov_prob_exp_plot_arma11_intercept.png", plot = conv_prob_arma11_intercept$gg_exp, width = 10, height = 8, dpi = 300)

################################################################################
####################### Convergence in Distribution ############################
conv_dist_arma11_intercept <- convergence_in_dist(results_arma11_intercept,
                                                  sample_sizes_vec = c(30, 35, 50, 100, 250, 500))

conv_dist_arma11_intercept$gg_exp
conv_dist_arma11_intercept$gg_normal
conv_dist_arma11_intercept$gg_normal_facet
conv_dist_arma11_intercept$gg_exp_facet

ggsave("gg_cov_dist_exp_plot_arma11_intercept.png", plot = conv_dist_arma11_intercept$gg_exp, width = 10, height = 8, dpi = 300)
ggsave("gg_cov_dist_normal_plot_arma11_intercept.png", plot = conv_dist_arma11_intercept$gg_normal, width = 10, height = 8, dpi = 300)

################################################################################
####################### Test Size ##############################################
test_size_arma11_intercept <- test_size(results_arma11_intercept,
                                        sample_sizes_vec = c(30, 35, 50, 100, 250, 500))

test_size_arma11_intercept$gg_normal
test_size_arma11_intercept$gg_exp

ggsave("gg_test_size_normal_plot_arma11_intercept.png", plot = test_size_arma11_intercept$gg_normal, width = 10, height = 8, dpi = 300)
ggsave("gg_test_size_exp_plot_arma11_intercept.png", plot = test_size_arma11_intercept$gg_exp, width = 10, height = 8, dpi = 300)








############################ARMA(1, 1) AR Coefficient#####################################


results_arma11_coefficient <- monte_carlo(
  c = 0, 
  theta = 0.5, 
  phi = 0.3,
  n_sim = 100,
  sample_sizes_vec = c(
    seq(10, 20, by = 1),
    seq(20, 50, by = 5),
    seq(50, 60, by = 1),
    seq(60, 90, by = 5),
    seq(90, 150, by = 10),
    seq(150, 300, by = 25),
    seq(300, 500, by = 50)
  ),
  delta_vec = c(0.25, 0.20, 0.15, 0.10),  # esse parâmetro não é usado na função
  order_arma = c(1, 0, 1),  # ARMA(1)
  coefficient = "ar1"
)


################################################################################
####################### Convergence in Probability #############################
conv_prob_arma11_coefficient <- convergence_in_prob(results_arma11_coefficient,
                                                    delta_vec = c(0.25, 0.20, 0.15, 0.10),
                                                    sample_sizes_vec = c(10, 20, 50, 60, 90, 150, 300, 400, 500))

conv_prob_arma11_coefficient$gg_normal
conv_prob_arma11_coefficient$gg_exp

ggsave("gg_cov_prob_normal_plot_arma11_coefficient.png", plot = conv_prob_arma11_coefficient$gg_normal, width = 10, height = 8, dpi = 300)
ggsave("gg_cov_prob_exp_plot_arma11_coefficient.png", plot = conv_prob_arma11_coefficient$gg_exp, width = 10, height = 8, dpi = 300)

################################################################################
####################### Convergence in Distribution ############################
conv_dist_arma11_coefficient <- convergence_in_dist(results_arma11_coefficient,
                                                    sample_sizes_vec = c(30, 35, 50, 100, 250, 500))

conv_dist_arma11_coefficient$gg_exp
conv_dist_arma11_coefficient$gg_normal
conv_dist_arma11_coefficient$gg_normal_facet
conv_dist_arma11_coefficient$gg_exp_facet

ggsave("gg_cov_dist_exp_plot_arma11_coefficient.png", plot = conv_dist_arma11_coefficient$gg_exp, width = 10, height = 8, dpi = 300)
ggsave("gg_cov_dist_normal_plot_arma11_coefficient.png", plot = conv_dist_arma11_coefficient$gg_normal, width = 10, height = 8, dpi = 300)

################################################################################
####################### Test Size ##############################################
test_size_arma11_coefficient <- test_size(results_arma11_coefficient,
                                          sample_sizes_vec = c(30, 35, 50, 100, 250, 500))

test_size_arma11_coefficient$gg_normal
test_size_arma11_coefficient$gg_exp

ggsave("gg_test_size_normal_plot_arma11_coefficient.png", plot = test_size_arma11_coefficient$gg_normal, width = 10, height = 8, dpi = 300)
ggsave("gg_test_size_exp_plot_arma11_coefficient.png", plot = test_size_arma11_coefficient$gg_exp, width = 10, height = 8, dpi = 300)






############################ARMA(1, 1) MA Coefficient#####################################


results_armama11_coefficient <- monte_carlo(
  c = 0, 
  theta = 0.5, 
  phi = 0.3,
  n_sim = 100,
  sample_sizes_vec = c(
    seq(10, 20, by = 1),
    seq(20, 50, by = 5),
    seq(50, 60, by = 1),
    seq(60, 90, by = 5),
    seq(90, 150, by = 10),
    seq(150, 300, by = 25),
    seq(300, 500, by = 50)
  ),
  delta_vec = c(0.25, 0.20, 0.15, 0.10),  # esse parâmetro não é usado na função
  order_arma = c(1, 0, 1),  # ARMA(1)
  coefficient = "ma1"
)


################################################################################
####################### Convergence in Probability #############################
conv_prob_armama11_coefficient <- convergence_in_prob(results_armama11_coefficient,
                                                      delta_vec = c(0.25, 0.20, 0.15, 0.10),
                                                      sample_sizes_vec = c(10, 20, 50, 60, 90, 150, 300, 400, 500))

conv_prob_armama11_coefficient$gg_normal
conv_prob_armama11_coefficient$gg_exp

ggsave("gg_cov_prob_normal_plot_armama11_coefficient.png", plot = conv_prob_armama11_coefficient$gg_normal, width = 10, height = 8, dpi = 300)
ggsave("gg_cov_prob_exp_plot_armama11_coefficient.png", plot = conv_prob_armama11_coefficient$gg_exp, width = 10, height = 8, dpi = 300)

################################################################################
####################### Convergence in Distribution ############################
conv_dist_armama11_coefficient <- convergence_in_dist(results_armama11_coefficient,
                                                      sample_sizes_vec = c(30, 35, 50, 100, 250, 500))

conv_dist_armama11_coefficient$gg_exp
conv_dist_armama11_coefficient$gg_normal
conv_dist_armama11_coefficient$gg_normal_facet
conv_dist_armama11_coefficient$gg_exp_facet

ggsave("gg_cov_dist_exp_plot_armama11_coefficient.png", plot = conv_dist_armama11_coefficient$gg_exp, width = 10, height = 8, dpi = 300)
ggsave("gg_cov_dist_normal_plot_armama11_coefficient.png", plot = conv_dist_armama11_coefficient$gg_normal, width = 10, height = 8, dpi = 300)

################################################################################
####################### Test Size ##############################################
test_size_armama11_coefficient <- test_size(results_armama11_coefficient,
                                            sample_sizes_vec = c(30, 35, 50, 100, 250, 500))

test_size_armama11_coefficient$gg_normal
test_size_armama11_coefficient$gg_exp

ggsave("gg_test_size_normal_plot_armama11_coefficient.png", plot = test_size_armama11_coefficient$gg_normal, width = 10, height = 8, dpi = 300)
ggsave("gg_test_size_exp_plot_armama11_coefficient.png", plot = test_size_armama11_coefficient$gg_exp, width = 10, height = 8, dpi = 300)

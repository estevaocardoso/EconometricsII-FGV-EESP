################################################################################
#Estevão Cardoso - Pset 3#######################################################

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
library(boot)
library(scales)
library(geomtextpath)
library(broom)
library(patchwork)
library(vars)
#algumas bibliotecas talvez nção estejam instaladas, então só verificar


################################################################################
######################Question 1################################################
rm(list = ls())
set.seed(735115)

brazil_data <- read.csv("data\\data_brazil.csv")
print(colnames(brazil_data))

brazil_data


#########################question 1 overall#####################################
#ADL(2,1)

brazil_data_reduced <- brazil_data %>%
  filter(date >= 1942 & date <= 2019) %>%
  mutate(
    gdp_lag1 = lag(real_gdp_growth_pct, 1),
    gdp_lag2 = lag(real_gdp_growth_pct, 2),
    exr_lag1 = lag(exchange_rate_real_dolar_annual_average, 1)
  ) %>%
  drop_na(real_gdp_growth_pct, gdp_lag1, gdp_lag2, exr_lag1)

#estimating the model
adl21_model <- lm(real_gdp_growth_pct ~ gdp_lag1 + gdp_lag2 + exr_lag1, data = brazil_data_reduced)

#predict fitted values + confidence intervals 
pred_sample <- predict(adl21_model, newdata = brazil_data_reduced, interval = "confidence") %>%
  as.data.frame()

pred_sample$date <- brazil_data_reduced$date

# getting the data for 2020 prediction
gdp_2019 <- brazil_data$real_gdp_growth_pct[brazil_data$date == 2019]
gdp_2018 <- brazil_data$real_gdp_growth_pct[brazil_data$date == 2018]
exr_2019 <- brazil_data$exchange_rate_real_dolar_annual_average[brazil_data$date == 2019]

predict_2020 <- data.frame(
  gdp_lag1 = gdp_2019,
  gdp_lag2 = gdp_2018,
  exr_lag1 = exr_2019
)

pred_2020 <- predict(adl21_model, newdata = predict_2020, interval = "confidence") %>%
  as.data.frame()

pred_2020$date <- 2020

# combining sample predictions with 2020 forecast
pred_all <- bind_rows(pred_sample, pred_2020)

# combining the actual GDP for plotting
plot_data_adl21 <- brazil_data %>%
  filter(date >= 1944 & date <= 2019) %>%
  select(date, real_gdp_growth_pct) %>%
  bind_rows(data.frame(date = 2020, real_gdp_growth_pct = NA)) %>%
  left_join(pred_all, by = "date") %>%
  rename(
    fitted_adl21 = fit,
    lwr_adl21 = lwr,
    upr_adl21 = upr
  )

#the plot itself

adl21plot <- ggplot(plot_data_adl21, aes(x = date)) +
  #real GDP
  geom_line(aes(y = real_gdp_growth_pct, color = "Real GDP"), size = 1) +
  #fitted values + forecast
  geom_line(aes(y = fitted_adl21, color = "Fitted by ADL(2,1)"), size = 1) +
  #confidence band continuous for sample + forecast
  geom_ribbon(aes(ymin = lwr_adl21, ymax = upr_adl21), fill = "#FF4500", alpha = 0.2) +
  #forecast point for 2020 highlighted
  geom_point(data = filter(plot_data_adl21, date == 2020),
             aes(y = fitted_adl21), color = "black", size = 2) +
  scale_color_manual(
    name = NULL,
    values = c("Real GDP" = "#698B69", "Fitted by ADL(2,1)" = "#8B2500")
  ) +
  labs(
    #title = "ADL(2,1): Real GDP Growth, Fitted Values, and Forecast with Confidence Band",
    x = "Year",
    y = "Real GDP Growth (%)"
  ) +
  theme_minimal()+
  theme(legend.position = "bottom")

ggsave("adl21.png", plot = adl21plot, width = 8, height = 6, dpi = 300)






################################################################################


brazil_data_reduced <- brazil_data %>% filter(date >= 1942 & date <= 2019)

# defininf the dataframe and variables
brazil_model_data <- brazil_data_reduced %>%
  mutate(
    gdp_lag1 = lag(real_gdp_growth_pct, 1),
    gdp_lag2 = lag(real_gdp_growth_pct, 2),
    exr_lag1 = lag(exchange_rate_real_dolar_annual_average, 1),
    inf_lag1 = lag(ipc_fipe_pct, 1),
    inf_lag2 = lag(ipc_fipe_pct, 2)
  ) %>% drop_na()

# this is just a function to fit model and collect results
fit_model <- function(formula, name, forecast_data) {
  model <- lm(formula, data = brazil_model_data)
  pred <- predict(model, newdata = forecast_data, interval = "confidence")
  
  coefs <- tidy(model) %>%
    mutate(term = paste0(term, "\n(", round(std.error, 3), ")")) %>%
    dplyr::select(term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    mutate(Model = name)
  
  info <- data.frame(
    Model = name,
    AIC = AIC(model),
    BIC = BIC(model),
    Forecast_2020 = pred[1],
    Forecast_Lower = pred[2],
    Forecast_Upper = pred[3]
  )
  
  return(list(model = model, coefs = coefs, info = info, forecast = pred))
}

# constructing a dataframe to store a bunch of lag variables
forecast_2020 <- data.frame(
  gdp_lag1 = brazil_data$real_gdp_growth_pct[brazil_data$date == 2019],
  gdp_lag2 = brazil_data$real_gdp_growth_pct[brazil_data$date == 2018],
  exr_lag1 = brazil_data$exchange_rate_real_dolar_annual_average[brazil_data$date == 2019],
  inf_lag1 = brazil_data$ipc_fipe_pct[brazil_data$date == 2019],
  inf_lag2 = brazil_data$ipc_fipe_pct[brazil_data$date == 2018]
)

# fitting all the models together
mdl1 <- fit_model(real_gdp_growth_pct ~ gdp_lag1 + gdp_lag2 + exr_lag1, "ADL(2,1)", forecast_2020)
mdl2 <- fit_model(real_gdp_growth_pct ~ gdp_lag1 + gdp_lag2 + inf_lag1 + inf_lag2, "ADL(2,2)", forecast_2020)
mdl3 <- fit_model(real_gdp_growth_pct ~ gdp_lag1 + gdp_lag2 + exr_lag1 + inf_lag1, "General", forecast_2020)
mdl4 <- fit_model(real_gdp_growth_pct ~ gdp_lag1 + gdp_lag2, "ARMA(2,0)", forecast_2020)

#now combining
all_coefs_long <- bind_rows(
  mdl1$coefs %>% mutate(Model = "Model 1"),
  mdl2$coefs %>% mutate(Model = "Model 2"),
  mdl3$coefs %>% mutate(Model = "Model 3"),
  mdl4$coefs %>% mutate(Model = "Model 4")
) %>%
  pivot_longer(cols = -Model, names_to = "Term", values_to = "Value") %>%
  pivot_wider(names_from = Model, values_from = Value)

# joining (AIC, BIC, forecast, etc.) ding the same s***
all_infos_long <- bind_rows(
  mdl1$info %>% mutate(Model = "Model 1"),
  mdl2$info %>% mutate(Model = "Model 2"),
  mdl3$info %>% mutate(Model = "Model 3"),
  mdl4$info %>% mutate(Model = "Model 4")
) %>%
  pivot_longer(cols = -Model, names_to = "Term", values_to = "Value") %>%
  pivot_wider(names_from = Model, values_from = Value)

#putting everythin in one table

#percebe que a tabela que sai daqui é meio ruim, então vc vai ter que pedir pro
#chatgpt só deixar ela no formato daquelas que estão no overleaf ou no formato que 
#tu quiser.

final_results <- bind_rows(all_infos_long, all_coefs_long)
print(xtable(final_results, type = "latex",
             caption = "Results for GDP forecast using DL and AR models",
             digits = 2,
             include.rownames = FALSE))


#this function is to plot individually each graph for each model in one single graph
plot_models_individualy <- function(model_obj, name) {
  fitted_vals <- model_obj$model$fitted.values
  conf_int <- predict(model_obj$model, interval = "confidence")
  plot_data <- brazil_model_data %>%
    dplyr::select(date, real_gdp_growth_pct) %>%
    slice_tail(n = length(fitted_vals)) %>%
    mutate(
      fitted = fitted_vals,
      lwr = conf_int[, "lwr"],
      upr = conf_int[, "upr"]
    )
  

  plot_data <- bind_rows(plot_data, data.frame(
    date = 2020,
    real_gdp_growth_pct = NA,
    fitted = model_obj$forecast[1],
    lwr = model_obj$forecast[2],
    upr = model_obj$forecast[3]
  ))
  
  p <- ggplot(plot_data, aes(x = date)) +
    geom_line(aes(y = real_gdp_growth_pct, color = "Real GDP"), size = 1) +
    geom_line(aes(y = fitted, color = "Fitted"), size = 1) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#FF4500", alpha = 0.2) +
    geom_point(data = filter(plot_data, date == 2020),
               aes(y = fitted), color = "black", size = 2) +
    scale_color_manual(
      values = c("Real GDP" = "#698B69", "Fitted" = "#8B2500")
    ) +
    theme_minimal() +
    labs(title = name, x = "Year", y = "Real GDP Growth (%)") +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  
  ggsave(paste0(tolower(gsub("[()]", "", gsub(" ", "", name))), ".png"), plot = p, width = 8, height = 5)
  
  return(p)
}

p1 <- plot_models_individualy(mdl1, "ADL(2,1)")
p2 <- plot_models_individualy(mdl2, "ADL(2,2)")
p3 <- plot_models_individualy(mdl3, "General")
p4 <- plot_models_individualy(mdl4, "ARMA(2,0)")

#doing a 2x2 grid layout to combine all the graphs in only because i have OCD
combined_plot <- (p1 | p2) / (p3 | p4)  # 
ggsave("all_models_plot_question1.png", plot = combined_plot, width = 16, height = 10)



#this other function is to plot every model forecast in one single figure (not grid)

comparison_plot_data <- brazil_model_data %>%
  dplyr::select(date, real_gdp_growth_pct) %>%
  slice_tail(n = length(mdl1$model$fitted.values)) %>%
  mutate(
    fitted_adl21 = mdl1$model$fitted.values,
    fitted_adl22 = mdl2$model$fitted.values,
    fitted_general = mdl3$model$fitted.values,
    fitted_ar = mdl4$model$fitted.values
  )


comparison_plot_data <- bind_rows(comparison_plot_data, data.frame(
  date = 2020,
  real_gdp_growth_pct = NA,
  fitted_adl21 = mdl1$forecast[1],
  fitted_adl22 = mdl2$forecast[1],
  fitted_general = mdl3$forecast[1],
  fitted_ar = mdl4$forecast[1]
))

plot_long <- comparison_plot_data %>%
  pivot_longer(cols = starts_with("fitted"),
               names_to = "Model", values_to = "Fitted") %>%
  mutate(
    Model = case_when(
      Model == "fitted_adl21" ~ "ADL(2,1)",
      Model == "fitted_adl22" ~ "ADL(2,2)",
      Model == "fitted_general" ~ "General",
      Model == "fitted_ar" ~ "ARMA(2,0)",
      TRUE ~ Model
    )
  )

comparasion_plot <- ggplot() +
  
  geom_line(data = comparison_plot_data,
            aes(x = date, y = real_gdp_growth_pct, color = "Real GDP"), size = 0.8) +
  geom_line(data = plot_long,
            aes(x = date, y = Fitted, color = Model), size = 1) +
  geom_point(data = filter(plot_long, date == 2020),
             aes(x = date, y = Fitted, color = Model), size = 2.5, shape = 21, fill = "white") +
  
  scale_color_manual(
    values = c(
      "Real GDP" = "black",
      "ADL(2,1)" = "#1b9e77",
      "ADL(2,2)" = "#d95f02",
      "General" = "#7570b3",
      "ARMA(2,0)" = "#e7298a"
    )
  ) +
  labs(
    #title = "Comparison of GDP Growth Fitted Values and 2020 Forecasts by Model",
    x = "Year",
    y = "Real GDP Growth (%)",
    color = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("comparasion_all_models_plot_question1.png", plot = comparasion_plot, width = 10, height = 7, dpi = 300)





################################################################################
#############################Question 6 - 1#####################################

#setting the dataframe properly
brazil_data_var <- brazil_data %>%
  filter(date >= 1942 & date <= 2019) %>%
  dplyr::select(
    date,
    real_gdp_growth_pct,
    exchange_rate_real_dolar_annual_average,
    ipc_fipe_pct
  ) %>%
  drop_na()

print(colnames(brazil_data_var))

#defining a ts
ts_var_data <- ts(
  brazil_data_var[, -1],   
  start = 1942,
  end = 2019,
  frequency = 1
)

#function to prepare dataframe with fitted values, ci and forecast for plotting
prepare_plot_data_var <- function(var_model, var_name = "real_gdp_growth_pct") {
  fitted_var <- fitted(var_model)[, var_name]
  
  mean_se <- mean(abs(summary(var_model)$varresult[[var_name]]$residuals))
  fitted_years <- (start(ts_var_data)[1] + var_model$p):(end(ts_var_data)[1])
  
  plot_data <- data.frame(
    date = fitted_years,
    real_gdp_growth_pct = brazil_data_var %>%
      filter(date %in% fitted_years) %>%
      pull(real_gdp_growth_pct),
    fitted = fitted_var,
    lwr = fitted_var - 1.96 * mean_se,
    upr = fitted_var + 1.96 * mean_se
  )
  
  forecast_var <- predict(var_model, n.ahead = 1, ci = 0.95)
  fcst_gdp <- forecast_var$fcst[[var_name]]
  
  forecast_2020 <- data.frame(
    date = 2020,
    real_gdp_growth_pct = NA,
    fitted = fcst_gdp[1],
    lwr = fcst_gdp[2],
    upr = fcst_gdp[3]
  )
  
  plot_data <- bind_rows(plot_data, forecast_2020)
  
  return(list(plot_data = plot_data, forecast = fcst_gdp))
}

#function to plot individual var graphs
plot_var_model <- function(plot_data, model_name) {
  p <- ggplot(plot_data, aes(x = date)) +
    geom_line(aes(y = real_gdp_growth_pct, color = "Real GDP"), size = 1) +
    geom_line(aes(y = fitted, color = paste0("Fitted by VAR(", model_name, ")")), size = 1) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#E066FF", alpha = 0.2) +
    geom_point(
      data = filter(plot_data, date == 2020),
      aes(y = fitted), color = "black", size = 2
    ) +
    scale_color_manual(
      name = NULL,
      values = setNames(
        c("#CD4F39", "#7A378B"),
        c("Real GDP", paste0("Fitted by VAR(", model_name, ")"))
      )
    ) +
    labs(x = "Year", y = "Real GDP Growth (%)") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}


#estimate var(2)
var_model2 <- VAR(ts_var_data, p = 2, type = "const")

plot_data_var2_res <- prepare_plot_data_var(var_model2)

#estimate var(3)
var_model3 <- VAR(ts_var_data, p = 3, type = "const")

plot_data_var3_res <- prepare_plot_data_var(var_model3)

#estimate var(1)
var_model1 <- VAR(ts_var_data, p = 1, type = "const")
plot_data_var1_res <- prepare_plot_data_var(var_model1)






#plotting the graphs
p1 <- plot_var_model(plot_data_var1_res$plot_data, 1)
p2 <- plot_var_model(plot_data_var2_res$plot_data, 2)
p3 <- plot_var_model(plot_data_var3_res$plot_data, 3)

combined_plot <- p1 | p2 | p3
ggsave("all_var_models_plot.png", plot = combined_plot, width = 18, height = 5)



#######################VAR(1,2,3) reduced form######################################

#Here I split the estimation of each reduced for for each VAR(p). So in each stage
#you advance, you will get a table will 3 reduced forms for one VAR(P). For ex:
#for VAR(1) the result will be the reduced form of GDP, Iflation and EXR along this
#BIC and AIC criteria and forecast and lower and upper bounds for gdp. So
#three reduced form for one VAR(1). So the other models will follow this method.


#as tabelas são geradas pra cada VAR(p) e bem desorganizadas, então o melhor é vc
#pegar elas com o codigo latex e jogar no chatgpt que ele organiza direitinho.
#só toma cuidado depois pra ver se ele não mudou nenhum coeficiente. ou pega um
#modelo das que estão no overeleaf e diz pra ele deixar a que vc gerar aqui igual.


#######Var (1)
extract_all_var_coefs <- function(var_model) {
  summary_var <- summary(var_model)
  coefs_list <- lapply(names(summary_var$varresult), function(var_name) {
    coefs <- summary_var$varresult[[var_name]]$coefficients
    coefs_df <- as.data.frame(coefs)
    coefs_df <- coefs_df %>%
      rownames_to_column(var = "Term") %>%
      mutate(
        Variable = var_name,
        StdError = round(`Std. Error`, 3),
        Estimate = round(Estimate, 3),
        Term = paste0(Term, " (", StdError, ")")
      ) %>%
      dplyr::select(Variable, Term, Estimate)
    return(coefs_df)
  })
  return(bind_rows(coefs_list))
}

coefs_var1 <- extract_all_var_coefs(var_model1)

aic_bic_df <- purrr::map_dfr(names(var_model1$varresult), function(var_name) {
  model_eq <- var_model1$varresult[[var_name]] 
  tibble(
    Term = c("AIC", "BIC"),
    Variable = var_name,
    Value = c(AIC(model_eq), BIC(model_eq))
  )
})


fc_var1 <- predict(var_model1, n.ahead = 1, ci = 0.95)
gdp_fc <- fc_var1$fcst$real_gdp_growth_pct 
forecast_df <- tibble(
  Variable = "real_gdp_growth_pct",
  Term = c("GDP Forecast 2020", "Forecast Lower Bound", "Forecast Upper Bound"),
  Estimate = round(c(gdp_fc[1], gdp_fc[2], gdp_fc[3]), 3)
)
final_var1_df <- bind_rows(aic_bic_df, coefs_var1, forecast_df) %>%
  arrange(Variable, match(Term, c("AIC", "BIC", "GDP Forecast 2020", 
                                  "Forecast Lower Bound", "Forecast Upper Bound")))


print(xtable(final_var1_df,
             type = "latex",
             caption = "Reduced Form VAR(1)",
             digits = 3,
             include.rownames = FALSE))



########Var (2)


coefs_var2 <- extract_all_var_coefs(var_model2)

aic_bic_df <- purrr::map_dfr(names(var_model2$varresult), function(var_name) {
  model_eq <- var_model2$varresult[[var_name]] 
  tibble(
    Term = c("AIC", "BIC"),
    Variable = var_name,
    Value = c(AIC(model_eq), BIC(model_eq))
  )
})


fc_var2 <- predict(var_model2, n.ahead = 1, ci = 0.95)
gdp_fc <- fc_var2$fcst$real_gdp_growth_pct 
forecast_df <- tibble(
  Variable = "real_gdp_growth_pct",
  Term = c("GDP Forecast 2020", "Forecast Lower Bound", "Forecast Upper Bound"),
  Estimate = round(c(gdp_fc[1], gdp_fc[2], gdp_fc[3]), 3)
)

final_var2_df <- bind_rows(aic_bic_df, coefs_var2, forecast_df) %>%
  arrange(Variable, match(Term, c("AIC", "BIC", "GDP Forecast 2020", 
                                  "Forecast Lower Bound", "Forecast Upper Bound")))

print(xtable(final_var2_df,
             type = "latex",
             caption = "Reduced Form VAR(2)",
             digits = 3,
             include.rownames = FALSE))


########Var (3)


coefs_var3 <- extract_all_var_coefs(var_model3)

aic_bic_df <- purrr::map_dfr(names(var_model3$varresult), function(var_name) {
  model_eq <- var_model3$varresult[[var_name]] 
  tibble(
    Term = c("AIC", "BIC"),
    Variable = var_name,
    Value = c(AIC(model_eq), BIC(model_eq))
  )
})

fc_var3 <- predict(var_model3, n.ahead = 1, ci = 0.95)
gdp_fc <- fc_var3$fcst$real_gdp_growth_pct 
forecast_df <- tibble(
  Variable = "real_gdp_growth_pct",
  Term = c("GDP Forecast 2020", "Forecast Lower Bound", "Forecast Upper Bound"),
  Estimate = round(c(gdp_fc[1], gdp_fc[2], gdp_fc[3]), 3)
)

final_var3_df <- bind_rows(aic_bic_df, coefs_var3, forecast_df) %>%
  arrange(Variable, match(Term, c("AIC", "BIC", "GDP Forecast 2020", 
                                  "Forecast Lower Bound", "Forecast Upper Bound")))

print(xtable(final_var3_df,
             type = "latex",
             caption = "Reduced Form VAR(3)",
             digits = 3,
             include.rownames = FALSE))






################################################################################
#############################Question 6 - 2#####################################

#defining another df to dont mess anything
irf_brazil <- brazil_data

irf_brazil <- irf_brazil %>%
  rename(
    GDP = real_gdp_growth_pct,
    EXR = exchange_rate_real_dolar_annual_average,
    INFL = ipc_fipe_pct)

Y <- irf_brazil[, c("INFL", "EXR", "GDP")]
colnames(Y) <- c("INFL", "EXR", "GDP")

#we dont need but anyways, estimating the VAR(2)
Y <- na.omit(Y)

var_model <- VAR(Y, p = 2, type = "const")


irf_struct <- irf(
  var_model,
  n.ahead = 5,
  boot = TRUE,
  ci = 0.90,
  runs = 1000,
  ortho = TRUE
)


#function to plot the IPR graph 


plot_irf <- function(irf_obj, filename, width, height, dpi) {
  variables <- colnames(Y)
  shocks <- variables
  horizons <- 0:(nrow(irf_obj$irf[[1]]) - 1)
  
  irf_data <- data.frame()
  for (resp in variables) {
    for (shock in shocks) {
      irf_values <- irf_obj$irf[[shock]][, resp]
      lower <- irf_obj$Lower[[shock]][, resp]
      upper <- irf_obj$Upper[[shock]][, resp]
      
      temp <- data.frame(
        Response = resp,
        Shock = shock,
        Horizon = horizons,
        IRF = irf_values,
        Lower = lower,
        Upper = upper
      )
      irf_data <- rbind(irf_data, temp)
    }
  }
  

  irf_data$FacetLabel <- paste0("Response: ", irf_data$Response, " | Shock: ", irf_data$Shock)
 

  ordered_labels <- irf_data %>%
    distinct(FacetLabel, Response, Shock) %>%
    arrange(
      factor(Response, levels = c("INFL", "EXR", "GDP")),
      factor(Shock, levels = c("INFL", "EXR", "GDP"))
    ) %>%
    pull(FacetLabel)

  irf_data$FacetLabel <- factor(irf_data$FacetLabel, levels = ordered_labels)
 
  plot <- ggplot(irf_data, aes(x = Horizon, y = IRF)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey60", alpha = 0.4) +
    geom_line(size = 1, color = "#7A378B") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    facet_wrap(~ FacetLabel, scales = "free_y", ncol = 3) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey80", size = 0.3),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(size = 11, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    ) +
    labs(y = "Impulse Response", x = "Year")
  ggsave(filename = filename, plot = plot, width = width, height = height, dpi = dpi)
}

plot_irf(irf_struct, filename = "irf_plot.png", width = 12, height = 9, dpi = 300)


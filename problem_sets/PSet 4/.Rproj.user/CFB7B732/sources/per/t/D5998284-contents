################################################################################
#Estevão Cardoso - Pset 4#######################################################

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

brazil_data <- read.csv("data\\brazil_data.csv")
usa_data <- read.csv("data\\usa_data.csv")

colnames(brazil_data)
str(brazil_data)
colnames(usa_data)


#########################item1#####################################

#we first adjust the colnames and then we define the period
colnames(brazil_data) <- c(
  "date",
  "exchange_rate",
  "ipca_index",
  "drop"
)

brazil_data <- brazil_data[, 1:3]

colnames(usa_data) <- c(
  "date",
  "us_cpi_index"
)

brazil_data <- brazil_data %>%
  mutate(
    date = as.yearmon(as.character(date), "%Y.%m"),
    date = as.Date(date)
  )

usa_data <- usa_data %>%
  mutate(
    date = as.Date(date)
  )

start_date <- as.Date("1995-01-01")
end_date <- as.Date("2019-12-31")


brazil_data_flt <- brazil_data %>%
  filter(date >= start_date & date <= end_date)

usa_data_flt <- usa_data %>%
  filter(date >= start_date & date <= end_date)

#monthly series for brazil timeserie format
brazil_ts <- ts(
  brazil_data_flt %>%
    dplyr::select(exchange_rate, ipca_index) %>%
    as.matrix(),
  start = c(1995, 1),
  end = c(2019, 12),
  frequency = 12
)
#monthly series for EUA timeserie format
usa_ts <- ts(
  usa_data_flt$us_cpi_index,
  start = c(1995, 1),
  end = c(2019, 12),
  frequency = 12
)


#########################the answer for item2 is here###########################

brazil_data_flt <- brazil_data_flt %>%
  mutate(
    exchange_rate = as.numeric(exchange_rate),
    ipca_index = as.numeric(ipca_index)
  )

usa_data_flt <- usa_data_flt %>%
  mutate(
    us_cpi_index = as.numeric(us_cpi_index)
  )

#calculating the log-variation since Jan/1995

#brazil
exchange_rate_base <- log(brazil_data_flt$exchange_rate[1])
ipca_index_base <- log(brazil_data_flt$ipca_index[1])

brazil_data_flt <- brazil_data_flt %>%
  mutate(
    log_change_exchange = 100 * (log(exchange_rate) - exchange_rate_base),
    log_change_ipca = 100 * (log(ipca_index) - ipca_index_base)
  )


#EUA
us_cpi_index_base <- log(usa_data_flt$us_cpi_index[1])

usa_data_flt <- usa_data_flt %>%
  mutate(
    log_change_us_cpi = 100 * (log(us_cpi_index) - us_cpi_index_base)
  )


#Brazil timeserie format
exchange_rate_base <- log(brazil_ts[1, 1])
ipca_index_base <- log(brazil_ts[1, 2])

log_change_exchange_ts <- 100 * (log(brazil_ts[, 1]) - exchange_rate_base)
log_change_ipca_ts <- 100 * (log(brazil_ts[, 2]) - ipca_index_base)

#USA timeserie format
us_cpi_index_base <- log(usa_ts[1])

log_change_us_cpi_ts <- 100 * (log(usa_ts) - us_cpi_index_base)
###############################item4############################################

#we need to combine both dataframes to get the cointegrated variable 
zt_data <- brazil_data_flt %>%
  dplyr::select(date, log_change_exchange, log_change_ipca) %>%
  left_join(
    usa_data_flt %>% dplyr::select(date, log_change_us_cpi),
    by = "date"
  )

zt_data <- zt_data %>%
  mutate(
    z_t =  log_change_ipca - log_change_exchange - log_change_us_cpi
  )
zt_ts <- ts(zt_data$z_t, start = c(1995, 1), frequency = 12)


#calculating the cointegratioon vector using timeserie format

zt_ts <-  log_change_ipca_ts - log_change_exchange_ts - log_change_us_cpi_ts

#creating a new dataframe to store the serie because there was a problem with the dates
zt_data <- data.frame(
  date = seq(as.Date("1995-01-01"), as.Date("2019-12-01"), by = "month"),
  log_change_exchange = log_change_exchange_ts,
  log_change_ipca = log_change_ipca_ts,
  log_change_us_cpi = log_change_us_cpi_ts,
  z_t = zt_ts
)

#############################item5##############################################

#we then do a graph usint the timeserie version of the data

question1item5 <- ggplot(zt_data, aes(x = date)) +
  geom_line(aes(y = log_change_exchange, color = "Y1"), size = 1) +
  geom_line(aes(y = log_change_ipca, color = "Y2"), size = 1) +
  geom_line(aes(y = log_change_us_cpi, color = "Y3"), size = 1) +
  labs(
    #title = expression("Transformed Series " ~ Y[k,t] ~ " for " ~ k %in% "{1,2,3}"),
    x = "Date",
    y = expression(Y[k,t]),
    color = "Variable"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c(
      "Y1" = "#c19a6b",   
      "Y2" = "#006400",  
      "Y3" = "#004e64"    
    ),
    labels = c(
      "Y1" = "Exchange Rate (BRL/USD)",
      "Y2" = "Brazil IPCA Index",
      "Y3" = "US CPI Index"
    )
  )

ggsave(
  filename = "quesion1item5.png",  
  plot = question1item5,
  width = 10,   
  height = 6,   
  dpi = 300    
)


#################################item6##########################################
#defining the series to be used in the ADF
Y1_ts <- log_change_exchange_ts
Y2_ts <- log_change_ipca_ts
Y3_ts <- log_change_us_cpi_ts

test_Y1_exr <- ur.df(Y1_ts, type = "drift", selectlags = "BIC")
test_Y2_ipca <- ur.df(Y2_ts, type = "trend", selectlags = "BIC")
test_Y3_cpi <- ur.df(Y3_ts, type = "trend", selectlags = "BIC")
summary(test_Y3_cpi)

#creating a function to extract the coeficients form the tests

extract_test_stats_from_summary <- function(test_obj) {
  
  summary_text <- capture.output(summary(test_obj))
  line_idx <- grep("Value of test-statistic is:", summary_text)
  
  if(length(line_idx) == 0) {
    warning("Linha com 'Value of test-statistic is:' não encontrada no summary")
    return(NULL)
  }
  
  line <- summary_text[line_idx]
  nums_str <- sub(".*Value of test-statistic is:\\s*", "", line)
  nums <- as.numeric(strsplit(nums_str, "\\s+")[[1]])
  
  return(nums)
}


get_full_adf_stats <- function(test_obj, var_name) {
  crit_vals <- test_obj@cval
  test_names <- rownames(crit_vals)
  test_stats <- extract_test_stats_from_summary(test_obj)
  if(length(test_stats) < length(test_names)) {
    test_stats <- c(test_stats, rep(NA, length(test_names) - length(test_stats)))
  }
  
  df <- data.frame(
    Variable = rep(var_name, length(test_names)),
    Test = test_names,
    Statistic = round(test_stats, 4),
    Critical_1pct = round(crit_vals[, "1pct"], 4),
    Critical_5pct = round(crit_vals[, "5pct"], 4),
    Critical_10pct = round(crit_vals[, "10pct"], 4),
    stringsAsFactors = FALSE
  )
  
  return(df)
}

full_adf_table <- rbind(
  get_full_adf_stats(test_Y1_exr, "Exchange Rate"),
  get_full_adf_stats(test_Y2_ipca, "Brazil IPCA"),
  get_full_adf_stats(test_Y3_cpi, "US CPI")
)

xtable(full_adf_table)



#################################item7##########################################

zt_plot <- ggplot(zt_data, aes(x = date)) +
  geom_line(aes(y = z_t, color = "z_t"), size = 1) +
  labs(
    x = "Date",
    y = expression(z[t]),
    color = "Serie"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c("z_t" = "#800000"),  # Vinho escuro, elegante
    labels = c("z_t" = expression(z[t]))
  )

# Exporta imagem com mesmo padrão de qualidade
ggsave(
  filename = "zt_series_plot.png",
  plot = zt_plot,
  width = 10,
  height = 6,
  dpi = 300
)

###########################Item8################################################

adf_zt <- ur.df(zt_data$z_t, type = "drift", selectlags = "BIC")
summary(adf_zt)

zttable <- get_full_adf_stats(adf_zt, "Cointegrated $Z_t$")

xtable(zttable)


##########################Question2#############################################
################################################################################

#######################item 2###################################################

#we will use the data contained in the zt_data dataframe
zt_data <- data.frame(
  date = seq(as.Date("1995-01-01"), as.Date("2019-12-01"), by = "month"),
  log_change_exchange = log_change_exchange_ts,
  log_change_ipca = log_change_ipca_ts,
  log_change_us_cpi = log_change_us_cpi_ts,
  z_t = zt_ts
)

#running a ols
unknown_zt <- lm(log_change_exchange ~ log_change_ipca + log_change_us_cpi, data = zt_data)
summary(unknown_zt)

xtable(unknown_zt)
#getting the residuals
resid_zt <- residuals(unknown_zt)

#testing wiht ADF
adf.test(residuos)

#plotting a graph for the residuals 
resid_data <- data.frame(
  date = zt_data$date,
  resid_zt = resid_zt
)

resid_plot <- ggplot(resid_data, aes(x = date)) +
  geom_line(aes(y = resid_zt, color = "resid_zt"), size = 1) +
  labs(
    x = "Date",
    y = expression(hat(u)[t]),
    color = "Series"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c("resid_zt" = "#455000"),
    labels = c("resid_zt" = expression(hat(u)[t]))
  )

ggsave(
  filename = "residuals_st_plot.png",
  plot = resid_plot,
  width = 10,
  height = 6,
  dpi = 300
)

#########################item3##################################################

PO_test_serie <- cbind(
  exchange_rate = log_change_exchange_ts,
  ipca_index = log_change_ipca_ts,
  us_cpi_index = log_change_us_cpi_ts
)

po_test <- ca.po(PO_test_serie, demean = "const", type = "Pu")

summary(po_test)

PO_test_table <- get_full_adf_stats(po_test, "PO Test")

xtable(PO_test_table)

#########################item5##################################################

#repeating the analysis for the two other variable orderings
#we test if our conclusions are robust to reordering the variables

#1st alternative ordering: log_change_ipca ~ log_change_exchange + log_change_us_cpi

#running OLS
alt1_model <- lm(log_change_ipca ~ log_change_exchange + log_change_us_cpi, data = zt_data)
summary(alt1_model)
xtable(alt1_model)

#getting residuals
alt1_resid <- residuals(alt1_model)

#ADF test on residuals
adf.test(alt1_resid)

#potting residuals
alt1_resid_data <- data.frame(
  date = zt_data$date,
  alt1_resid = alt1_resid
)

alt1_plot <- ggplot(alt1_resid_data, aes(x = date)) +
  geom_line(aes(y = alt1_resid, color = "alt1_resid"), size = 1) +
  labs(
    x = "Date",
    y = expression(hat(u)[t]),
    color = "Series"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c("alt1_resid" = "#1A73E8"),  
    labels = c("alt1_resid" = expression(hat(u)[t]))
  )

ggsave(
  filename = "residuals_alt1_plot.png",
  plot = alt1_plot,
  width = 10,
  height = 6,
  dpi = 300
)

#phillips-Ouliaris Test for first alternative ordering
PO_test_alt1 <- cbind(
  ipca_index = log_change_ipca_ts,
  exchange_rate = log_change_exchange_ts,
  us_cpi_index = log_change_us_cpi_ts
)

po_alt1 <- ca.po(PO_test_alt1, demean = "const", type = "Pu")
summary(po_alt1)


#2nd alternative ordering: log_change_us_cpi ~ log_change_exchange + log_change_ipca

#running OLS
alt2_model <- lm(log_change_us_cpi ~ log_change_exchange + log_change_ipca, data = zt_data)
summary(alt2_model)
xtable(alt2_model)

#getting residuals
alt2_resid <- residuals(alt2_model)

#ADF test on residuals
adf.test(alt2_resid)

#plotting residuals
alt2_resid_data <- data.frame(
  date = zt_data$date,
  alt2_resid = alt2_resid
)

alt2_plot <- ggplot(alt2_resid_data, aes(x = date)) +
  geom_line(aes(y = alt2_resid, color = "alt2_resid"), size = 1) +
  labs(
    x = "Date",
    y = expression(hat(u)[t]),
    color = "Series"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c("alt2_resid" = "#D93025"),  # vermelho elegante
    labels = c("alt2_resid" = expression(hat(u)[t]))
  )

ggsave(
  filename = "residuals_alt2_plot.png",
  plot = alt2_plot,
  width = 10,
  height = 6,
  dpi = 300
)

#phillips-Ouliaris Test for second alternative ordering
PO_test_alt2 <- cbind(
  us_cpi_index = log_change_us_cpi_ts,
  exchange_rate = log_change_exchange_ts,
  ipca_index = log_change_ipca_ts
)

po_alt2 <- ca.po(PO_test_alt2, demean = "const", type = "Pu")
summary(po_alt2)

#gttin the results of the PO tests together


PO_test_alt_table <- rbind(
  get_full_adf_stats(po_alt1, "Brazil IPCA"),
  get_full_adf_stats(po_alt2, "USA CPI")
)

xtable(PO_test_alt_table)



########################Question 2 test series em level#########################
################################################################################

#######################item 2###################################################

Y1_ts_level <- brazil_ts[, 1]  
Y2_ts_level <- brazil_ts[, 2] 
Y3_ts_level <- usa_ts          

df_level <- data.frame(
  date = seq(as.Date("1995-01-01"), as.Date("2019-12-01"), by = "month"),
  Y1_ts = Y1_ts_level,
  Y2_ts = Y2_ts_level,
  Y3_ts = Y3_ts_level
)

#running a ols
unknown_zt <- lm(Y1_ts ~ Y2_ts + Y3_ts, data = df)
summary(unknown_zt)

xtable(unknown_zt)
#getting the residuals
resid_zt <- residuals(unknown_zt)

#testing wiht ADF
adf.test(resid_zt)
summary(resid_zt)


#plotting a graph for the residuals 
resid_data <- data.frame(
  date = df_level$date,
  resid_zt = resid_zt
)

resid_plot <- ggplot(resid_data, aes(x = date)) +
  geom_line(aes(y = resid_zt, color = "resid_zt"), size = 1) +
  labs(
    x = "Date",
    y = expression(hat(u)[t]),
    color = "Series"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c("resid_zt" = "#455000"),
    labels = c("resid_zt" = expression(hat(u)[t]))
  )

ggsave(
  filename = "residuals_st_plot.png",
  plot = resid_plot,
  width = 10,
  height = 6,
  dpi = 300
)

#########################item3##################################################

PO_test_serie <- cbind(
  exchange_rate = Y1_ts,
  ipca_index = Y2_ts,
  us_cpi_index = Y3_ts
)

po_test <- ca.po(PO_test_serie, demean = "const", type = "Pu")

summary(po_test)

PO_test_table <- get_full_adf_stats(po_test, "PO Test")

xtable(PO_test_table)

#########################item5##################################################

#repeating the analysis for the two other variable orderings
#we test if our conclusions are robust to reordering the variables

#1st alternative ordering: log_change_ipca ~ log_change_exchange + log_change_us_cpi

#running OLS
alt1_model <- lm(Y2_ts ~ Y1_ts + Y3_ts, data = df_level)
summary(alt1_model)
xtable(alt1_model)

#getting residuals
alt1_resid <- residuals(alt1_model)

#ADF test on residuals
adf.test(alt1_resid)

#potting residuals
alt1_resid_data <- data.frame(
  date = df_level$date,
  alt1_resid = alt1_resid
)

alt1_plot <- ggplot(alt1_resid_data, aes(x = date)) +
  geom_line(aes(y = alt1_resid, color = "alt1_resid"), size = 1) +
  labs(
    x = "Date",
    y = expression(hat(u)[t]),
    color = "Series"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c("alt1_resid" = "#1A73E8"),  
    labels = c("alt1_resid" = expression(hat(u)[t]))
  )

ggsave(
  filename = "residuals_alt1_plot.png",
  plot = alt1_plot,
  width = 10,
  height = 6,
  dpi = 300
)

#phillips-Ouliaris Test for first alternative ordering
PO_test_alt1 <- cbind(
  ipca_index = Y2_ts,
  exchange_rate = Y1_ts,
  us_cpi_index = Y3_ts
)

po_alt1 <- ca.po(PO_test_alt1, demean = "const", type = "Pu")
summary(po_alt1)


#2nd alternative ordering: log_change_us_cpi ~ log_change_exchange + log_change_ipca

#running OLS
alt2_model <- lm(Y3_ts ~ Y1_ts + Y2_ts, data = df_level)
summary(alt2_model)
xtable(alt2_model)

#getting residuals
alt2_resid <- residuals(alt2_model)

#ADF test on residuals
adf.test(alt2_resid)

#plotting residuals
alt2_resid_data <- data.frame(
  date = df_level$date,
  alt2_resid = alt2_resid
)

alt2_plot <- ggplot(alt2_resid_data, aes(x = date)) +
  geom_line(aes(y = alt2_resid, color = "alt2_resid"), size = 1) +
  labs(
    x = "Date",
    y = expression(hat(u)[t]),
    color = "Series"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c("alt2_resid" = "#D93025"),  # vermelho elegante
    labels = c("alt2_resid" = expression(hat(u)[t]))
  )

ggsave(
  filename = "residuals_alt2_plot.png",
  plot = alt2_plot,
  width = 10,
  height = 6,
  dpi = 300
)

#phillips-Ouliaris Test for second alternative ordering
PO_test_alt2 <- cbind(
  us_cpi_index = Y3_ts,
  exchange_rate = Y1_ts,
  ipca_index = Y2_ts
)

po_alt2 <- ca.po(PO_test_alt2, demean = "const", type = "Pu")
summary(po_alt2)

#gettin the results of the PO tests together


PO_test_alt_table <- rbind(
  get_full_adf_stats(po_alt1, "Brazil IPCA"),
  get_full_adf_stats(po_alt2, "USA CPI")
)

xtable(PO_test_alt_table)


################################################################################
##########################question3#############################################

###############################item1############################################

johansen_data <- cbind(
  exchange_rate = Y1_ts_level,
  ipca_index = Y2_ts_level,
  us_cpi_index = Y3_ts_level
)


VARselect(johansen_data, lag.max = 12, type = "none")

johansen_test <- ca.jo(
  johansen_data,
  type = "eigen",    
  ecdet = "none",   
  K = 12              
)


summary(johansen_test)







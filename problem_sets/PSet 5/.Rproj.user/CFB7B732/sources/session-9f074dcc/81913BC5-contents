################################################################################
#Estevão Cardoso - Pset 5#######################################################

#Importing the libraries we need ###############################################
library(gmm)
library(MetBrewer)
library(forecast)
library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
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
library(readr)
library(lubridate)
#algumas bibliotecas talvez nção estejam instaladas, então só verificar


################################################################################
######################Question 1################################################

#item 2

rm(list = ls())

set.seed(220122)

n <- 100000
lambda0 <- 5
data <- rexp(n = n, rate = lambda0)




#moment function
gmm_moment <- function(theta, x){
  m1 <- 1 / theta - x
  m2 <- 1 / (theta**2) - (x - 1 / theta)**2
  f <- cbind(m1, m2)
  return(f)
}

#gradient of the moment function
gmm_gradient <- function(theta, x) {
  G1 <- -1 / (theta ** 2)
  G2 <- -2*mean(x) / (theta ** 2)
  G <- rbind(G1, G2)
  return(G)
}

#running the gmm using the brent method
gmm_result <- gmm(g = gmm_moment,
                  x = data,
                  t0 = c(1),                
                  gradv = gmm_gradient,
                  #wmatrix = "ident",       
                  method = "Brent",         
                  lower = 1,           
                  upper = 7)               

summary(gmm_result)

stargazer(gmm_result)


################################################################################
######################Question 1################################################

#item 4
x_arma <- arima.sim(n = n,
                    model = list(ar = c(0.2), ma = c(0.1, 0.1)))

#setting a dataframe with Yt (Y0), Y_{t-1} (Y1), Y_{t-3} (Y3) and the intercept

df <- data.frame(
  "Y0" = x_arma,
  "Y1" = c(NA, x_arma[1:(length(x_arma) - 1)]),
  "Y3" = c(NA, NA, NA, x_arma[1:(length(x_arma) - 3)]),
  "intercept" = rep(1, length(x_arma))
)

df <- df[4:nrow(df), ]


#running the gmm
gmm_arma <- gmm(
  g = df$Y0 ~ df$Y1,
  x = df$Y3
)

summary(gmm_arma)

stargazer(gmm_arma)


################################################################################
##############################Question3#########################################

#item 1

#improting the datasets

bvsp  <- read_csv("data/ABEV3.csv")
bbdc3 <- read_csv("data/BBDC3.csv")
abev3 <- read_csv("data/ABEV3.csv")
itub3 <- read_csv("data/ITUB3.csv")


#creating a function to calculate the daily log

calc_returns <- function(df, ticker) {
  df %>%
    mutate(
      Date = ymd(Date)
    ) %>%
    filter(Date >= as.Date("2021-01-01") & Date <= as.Date("2021-12-31")) %>%
    arrange(Date) %>%
    mutate(
      ret = 100 * (`Adj Close` / lag(`Adj Close`) - 1)  # retorno percentual simples * 100
    ) %>%
    dplyr::select(Date, ret) %>%
    rename(!!ticker := ret)
}

#using the function calc_returns to calculate the returns
bvsp_r  <- calc_returns(bvsp,  "r_bvsp")
bbdc3_r <- calc_returns(bbdc3, "r_bbdc3")
abev3_r <- calc_returns(abev3, "r_abev3")
itub3_r <- calc_returns(itub3, "r_itub3")

#calculating the returns by date

daily_returns_2021 <- bvsp_r %>%
  inner_join(bbdc3_r, by = "Date") %>%
  inner_join(abev3_r, by = "Date") %>%
  inner_join(itub3_r, by = "Date")



#joining the results and ploting a graph 

returns_all <- reduce(
  list(abev3_r, bbdc3_r, itub3_r, bvsp_r),
  full_join,
  by = "Date"
) %>%
  arrange(Date)

returns_long <- returns_all %>%
  pivot_longer(
    cols = starts_with("r_"),
    names_to = "Asset",
    values_to = "Return"
  )

pallet <- c("#D84315", "#FF6F00", "#FDD835", "#6D4C41")

returns_plot <- ggplot(returns_long, aes(x = Date, y = Return, color = Asset)) +
  geom_line(na.rm = TRUE, size = 0.5) +
  scale_color_manual(values = pallet) +
  theme_minimal() +
  labs(
    x = "Date",
    y = "Daily Return",
    color = "Asset"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_blank()
  )

ggsave("returns_plot.png", plot = returns_plot, width = 8, height = 5, dpi = 300)


#item 2

#importing the selic dataframe
selic <- read_delim("data/selic.csv", delim = ";", 
                    locale = locale(decimal_mark = ","),
                    trim_ws = TRUE)

#adjusting it properly
selic <- selic %>%
  mutate(Data = dmy(Data)) %>%      
  filter(year(Data) == 2021) 

selic <- selic %>%
  mutate(`Fator diário` = as.numeric(`Fator diário`))

selic <- selic %>%
  mutate(Return_Selic = 100 * (`Fator diário` - 1))

selic_daily <- selic %>%
  dplyr::select(Date = Data, Return_Selic)

#doing the merge with the main dataframe
daily_returns_selic_2021 <- daily_returns_2021 %>%
  inner_join(selic_daily, by = "Date")

daily_returns_selic_2021 <- daily_returns_selic_2021[2:nrow(daily_returns_selic_2021), ]
#item 3


excess_returns <- daily_returns_selic_2021 %>%
  transmute(
    date = Date,
    Market_Excess = r_bvsp - Return_Selic,
    Excess_ABEV3 = r_abev3 - Return_Selic,
    Excess_BBDC3 = r_bbdc3 - Return_Selic,
    Excess_ITUB3 = r_itub3 - Return_Selic
  )

#item 4

df_capm <- excess_returns %>% 
  mutate(`Rm-Rf` = Market_Excess,
         `R_ABEV-Rf` =Excess_ABEV3,
         `R_BBDC-Rf` = Excess_BBDC3,
         `R_ITUB-Rf` = Excess_ITUB3) %>% 
  dplyr::select(date, `Rm-Rf`, starts_with("R_"))



moment_capm <- function(theta, data) {
  alpha <- theta[1]
  beta  <- theta[2]
  
  Ri_Rf <- data[, 1]  
  Rm_Rf <- data[, 2]  
  
  m1 <- Ri_Rf - alpha - beta * Rm_Rf            
  m2 <- m1 * 1                                  
  m3 <- m1 * Rm_Rf                              
  
  cbind(m2, m3)                              
}

#model gmm for abev
data_abev <- df_capm %>%
  dplyr::select(`R_ABEV-Rf`, `Rm-Rf`) %>%
  as.matrix()
start_theta <- c(0, 0)
gmm_abev <- gmm(
  g = moment_capm,
  x = data_abev,
  t0 = start_theta
)

summary(gmm_abev)

#model gmm for bbdc
data_bbdc <- df_capm %>%
  dplyr::select(`R_BBDC-Rf`, `Rm-Rf`) %>%
  as.matrix()

gmm_bbdc <- gmm(
  g = moment_capm,
  x = data_bbdc,
  t0 = start_theta
)

summary(gmm_bbdc)

#model gmm for itub
data_itub <- df_capm %>%
  dplyr::select(`R_ITUB-Rf`, `Rm-Rf`) %>%
  as.matrix()

gmm_itub <- gmm(
  g = moment_capm,
  x = data_itub,
  t0 = start_theta
)
summary(gmm_itub)


stargazer(gmm_itub, gmm_bbdc, gmm_abev)

#gerenating a table for the results

createTexreg.gmm <- function(gmm_obj, model.name = "GMM Model") {
  coefs <- coef(gmm_obj)
  se <- sqrt(diag(vcov(gmm_obj)))
  pvalues <- 2 * (1 - pnorm(abs(coefs / se)))  
  
  texreg::createTexreg(
    coef.names = names(coefs),
    coef = coefs,
    se = se,
    pvalues = pvalues,
    model.name = model.name
  )
}

gmm_abev_texreg <- createTexreg.gmm(gmm_abev, model.name = "ABEV3")
gmm_bbdc_texreg <- createTexreg.gmm(gmm_bbdc, model.name = "BBDC3")
gmm_itub_texreg <- createTexreg.gmm(gmm_itub, model.name = "ITUB3")

screenreg(
  list(gmm_abev_texreg, gmm_bbdc_texreg, gmm_itub_texreg),
  digits = 4,
  stars = c(0.01, 0.05, 0.1)
)


#item 5

moment_capm_joint <- function(theta, data) {
  
  alpha1 <- theta[1]
  beta1  <- theta[2]
  alpha2 <- theta[3]
  beta2  <- theta[4]
  alpha3 <- theta[5]
  beta3  <- theta[6]
  
  
  Ri1_Rf <- data[, 1]  
  Ri2_Rf <- data[, 2]  
  Ri3_Rf <- data[, 3]  
  Rm_Rf  <- data[, 4]  
  
  
  m1 <- (Ri1_Rf - alpha1 - beta1 * Rm_Rf)
  m2 <- m1 * 1
  m3 <- m1 * Rm_Rf
  
  m4 <- (Ri2_Rf - alpha2 - beta2 * Rm_Rf)
  m5 <- m4 * 1
  m6 <- m4 * Rm_Rf
  
  m7 <- (Ri3_Rf - alpha3 - beta3 * Rm_Rf)
  m8 <- m7 * 1
  m9 <- m7 * Rm_Rf
  
  cbind(m2, m3, m5, m6, m8, m9)
}

data_joint <- df_capm %>%
  dplyr::select(`R_ABEV-Rf`, `R_BBDC-Rf`, `R_ITUB-Rf`, `Rm-Rf`) %>%
  as.matrix()

start_theta_joint <- c(0, 0, 0, 0, 0, 0)  

gmm_joint <- gmm(
  g = moment_capm_joint,
  x = data_joint,
  t0 = start_theta_joint
)


summary(gmm_joint)

#testing the hypotheis

n_restrictions <- 3

R_capm <- cbind(diag(n_restrictions), 
                matrix(0, nrow = n_restrictions, ncol = n_restrictions))

# Create the vector of null values
theta_null <- rep(0, n_restrictions)

# Test the null hypothesis of interest
print(linearHypothesis(
  model = gmm_joint,
  hypothesis.matrix = R_capm,
  rhs = theta_null,
  test = "Chisq"
))

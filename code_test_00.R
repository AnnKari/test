install.packages(c("fredr", "tseries", "urca", "vars", "forecast", "ggplot2", "dplyr", "lubridate", "stargazer", "usethis"
#FRED API access
# Time series analysis
# Unit root tests
# VAR models (Tutorial 2)
# ARIMA modeling & forecasting
# Plotting
# Data manipulation
# Date handling
# LaTeX tables
))

sapply(c("fredr", "tseries", "urca", "vars", "forecast", "ggplot2", "dplyr", "lubridate", "stargazer", "usethis"),require , character.only = TRUE)

library(fredr)
library(tseries)
library(urca)
library(forecast)
library(ggplot2)
library(lubridate)
library(usethis)

usethis::edit_r_environ()

fred_key= Sys.getenv("FRED_API_KEY")

library(fredr)

# Set your API key (reads from .Renviron)
fredr_set_key(Sys.getenv("FRED_API_KEY"))
# Download Real GDP (quarterly, seasonally adjusted)
gdp_raw = fredr(
series_id = "GDPC1", observation_start = as.Date("1947-01-01"), observation_end = as.Date("2024-12-31"), frequency = "q"
)
# Inspect the data
head(gdp_raw)

# Compute annualized quarterly GDP growth rate
gdp = gdp_raw[, c("date", "value")] # Keep only date and value 
colnames(gdp) = c("date", "gdp") # Rename columns
# Log-difference * 400 = annualized growth rate (%) 
gdp$growth <- c(NA, 400 * diff(log(gdp$gdp)))
# Remove the first observation (which is NA)
gdp <- gdp[-1, ]
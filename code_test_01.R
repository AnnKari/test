# Macroeconometrics

#This is a test

# Tutorial 1  : AR Model on US Real GDP Growth
# Author      : Juan Pablo Ugarte Checura
# Date        : April 2026

# This script walks through a complete AR modeling workflow:
#   1. Data acquisition from FRED
#   2. Data transformation (levels -> growth rates)
#   3. Visualization (time series plot, ACF, PACF)
#   4. Stationarity testing (KPSS)
#   5. Lag order selection (AIC, BIC)
#   6. AR model estimation
#   7. Diagnostic checks
#   8. Forecasting
#   9. Saving output

# Before running, make sure you have:
#   - Installed required packages (see below)
#   - Obtained FRED API key (https://fred.stlouisfed.org/docs/api/api_key.html)
#   - Open R terminal and modify your .Renviron file: usethis::edit_r_environ()
#   - Set your API key in .Renviron: FRED_API_KEY=your_key_here


# 0. Setup --------------------------------------------------------------------

# Install packages (run once, then comment out)
# install.packages(c(
#   "fredr",        # FRED API access
#   "tseries",      # Time series analysis
#   "urca",         # Unit root tests
#   "vars",         # VAR models (Tutorial 2)
#   "forecast",     # ARIMA modeling & forecasting
#   "ggplot2",      # Plotting
#   "dplyr",        # Data manipulation
#   "lubridate",    # Date handling
#   "stargazer"     # LaTeX tables
# ))

library(fredr)
library(tseries)
library(urca)
library(forecast)
library(ggplot2)
library(lubridate)
library(here)


# Set FRED API key (reads from .Renviron)
fredr_set_key(Sys.getenv("FRED_API_KEY"))

# Create output directories

dir.create(here("output/figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output/tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("data/processed"), recursive = TRUE, showWarnings = FALSE)


# 1. Data acquisition from FRED ------------------------------------------------

# Download Real GDP (quarterly, seasonally adjusted)
gdp_raw <- fredr(
  series_id         = "GDPC1",
  observation_start = as.Date("1947-01-01"),
  observation_end   = as.Date("2024-12-31"),
  frequency         = "q"
)

# Inspect the data
head(gdp_raw)


# 2. Data transformations ------------------------------------------------------

# Keep only date and value columns
gdp <- gdp_raw[, c("date", "value")]      # Keep only date and value
colnames(gdp) <- c("date", "gdp")         # Rename columns

# Log-difference * 400 = annualized growth rate (%)
gdp$growth <- c(NA, 400 * diff(log(gdp$gdp)))

# Remove the first observation (which is NA)
gdp <- gdp[-1, ]


# 3. Step 1: Visualize the series ----------------------------------------------

# Time series plot
p_ts <- ggplot(gdp, aes(x = date, y = growth)) +
  geom_line(color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(title = "US Real GDP Growth (Annualized)",
       x = NULL, y = "Percent") +
  theme_minimal()

print(p_ts)

# Convert to a ts object for time series functions
gdp_ts <- ts(gdp$growth,
             start = c(year(min(gdp$date)), quarter(min(gdp$date))),
             frequency = 4)

# Select a subset of the sample (e.g., 1980s up to pre-COVID)
gdp_ts <- window(gdp_ts, start = c(1986, 1), end = c(2019, 4))

# ACF and PACF plots (remember the issue of ACF's conf. int. seen in class!)
par(mfrow = c(1, 2))
acf(gdp_ts, main = "ACF of GDP Growth", lag.max = 20)
pacf(gdp_ts, main = "PACF of GDP Growth", lag.max = 20)
par(mfrow = c(1, 1))


# 4. Step 2: Test for stationarity ---------------------------------------------

# KPSS test (null: stationary)
kpss_result <- kpss.test(gdp_ts, null = "Level")
print(kpss_result)


# 5. Step 3: Select the lag order ----------------------------------------------

# Fit AR models for p = 1, ..., 8 and compare information criteria
max_p <- 8
ic_table <- data.frame(
  p   = 1:max_p,
  AIC = NA,
  BIC = NA,
  HQ  = NA
)

for (p in 1:max_p) {
  fit <- ar(gdp_ts, order.max = p,
            aic = FALSE, method = "ols")
  n <- length(gdp_ts) - p
  k <- p + 1
  ss <- sum(fit$resid^2, na.rm = TRUE)
  ic_table$AIC[p] <- n * log(ss / n) + 2 * k
  ic_table$BIC[p] <- n * log(ss / n) + k * log(n)
  ic_table$HQ[p]  <- n * log(ss / n) + 2 * k * log(log(n))
}

print(ic_table)


# 6. Step 4: Estimate the AR model ---------------------------------------------

# AIC, BIC, and HQ select p = 2
p_star <- 2

# Estimate AR(2) model
ar_model <- Arima(gdp_ts, order = c(p_star, 0, 0))
summary(ar_model)


# 7. Step 5: Diagnostic checking -----------------------------------------------

# Collect residuals
resid <- residuals(ar_model)

# Check autocorrelation (Ljung-Box test, H0: no autocorrelation up to lag 10)
Box.test(resid, lag = 10, type = "Ljung-Box", fitdf = p_star)

# Plot residuals
checkresiduals(ar_model)


# 8. Step 6: Forecasting -------------------------------------------------------

# Produce forecasts for 8 quarters ahead
fc <- forecast(ar_model, h = 8)
print(fc)

# Plot forecasts
p_fc <- autoplot(fc) +
  labs(title = "AR(2) Forecast of GDP Growth",
       x = NULL, y = "Percent (annualized)") +
  theme_minimal()

print(p_fc)


# 9. Saving output -------------------------------------------------------------

# Save figures
ggsave("output/figures/gdp_growth_ts.pdf", p_ts, width = 8, height = 4)

# Save ACF/PACF as a single PDF
pdf("output/figures/acf_pacf.pdf", width = 8, height = 3.5)
par(mfrow = c(1, 2))
acf(gdp_ts, main = "ACF of GDP Growth", lag.max = 20)
pacf(gdp_ts, main = "PACF of GDP Growth", lag.max = 20)
par(mfrow = c(1, 1))
dev.off()

# Save residual diagnostics
pdf("output/figures/residuals_diag.pdf", width = 8, height = 4)
checkresiduals(ar_model)
dev.off()

# Save forecast plot (same look as on-screen p_fc)
ggsave("output/figures/ar2_forecast.pdf", p_fc, width = 8, height = 4)

# Save model summary to text file
sink("output/tables/ar2_model_summary.txt")
summary(ar_model)
sink()

# Save data for reuse
write.csv(gdp, "data/processed/gdp_growth.csv", row.names = FALSE)


# EXERCISE: Try it yourself!

# 1. Download a different FRED series and repeat this workflow.
#
# 2. What happens if you use AIC instead of BIC to select the lag order?
#    Do the diagnostics improve?
#
# 3. Try estimating the model with a full and a truncated sample.
#    How do the AR coefficients and forecasts change?
#
# 4. What happens if you estimate the model on GDP in log-levels
#    Do the stationarity tests give a different result?

# end



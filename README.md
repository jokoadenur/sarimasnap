![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)
![CRAN Version](https://img.shields.io/badge/CRAN-0.1.0-brightgreen)
![Open Issues](https://img.shields.io/badge/open%20issues-0-brightgreen)
![License](https://img.shields.io/badge/License-MIT-blue)

<img width="1024" height="1024" alt="ChatGPT Image Feb 16, 2026, 08_45_53 PM" src="https://github.com/user-attachments/assets/78407e97-a484-440b-bcdf-5fbddd41651b" />


# sarimasnap

`sarimasnap` is an R package for **automatic Seasonal ARIMA (SARIMA) modeling**
using structured grid search combined with **statistical diagnostic filtering**.

The package is designed for:
- Official statistics
- Macroeconomic indicators
- Financial time series
- Production and administrative data
- Any quarterly or monthly time-series dataset

---

## Dataset Requirements

A time series will be analyzed by `sarimasnap()` if:

1. Data type: numeric vector or `ts` object  
2. Minimum observations: ≥ 20 recommended  
3. Ordered data (chronological order required)  
4. Missing values should be removed before modeling  

---

## Recommended Structure

### Option 1: ts object (Recommended)

```r
ts_data <- ts(data_vector, start = c(2010, 1), frequency = 4)
sarimasnap(ts_data) # h = 10 by default
```
### Option 2: ts object (adjusted forecast period)

```r
ts_data <- ts(data_vector, start = c(2010, 1), frequency = 4)
sarimasnap(ts_data, h = 15)
```
## What sarimasnap Does?

The function performs: (1) Seasonal detection; (2) 80/20 train-test split; (3) Grid search over SARIMA(p,d,q)(P,D,Q)[s]; (4) Ljung-Box residual filtering; (5) RMSE ranking (Top 10 models); (6) Stationarity testing (DF & ADF); (7) Final estimation on full dataset; (8) Residual diagnostics; (9) Forecast generation

## Installation

From Github `devtools::install_github("jokoadenur/sarimasnap")` then activate the library `library(sarimasnap)`

## Example 1

```r
library(sarimasnap)
data(AirPassengers)
result <- sarimasnap(AirPassengers, h = 12)
```
## Outputs
<img width="748" height="778" alt="image" src="https://github.com/user-attachments/assets/cdc1b68b-d6c2-465e-aae2-dc3a319d6654" />
<img width="1259" height="782" alt="image" src="https://github.com/user-attachments/assets/0e771f33-ed2b-4d02-bd6e-ec8aad6d2578" />

## Example 2
```r
# Install packages directly from GitHub (run once)
remotes::install_github("jokoadenur/sarimasnap")
remotes::install_github("jokoadenur/bpsuseR")

# Load required libraries
library(sarimasnap)
library(bpsuseR)
library(tidyr)
library(dplyr)

# Load quarterly GDP dataset (wide format: Qtr1, Qtr2, Qtr3, Qtr4)
dataku <- pdb_adhk_triwulanan

# Convert from wide format to long format
# This stacks Qtr1–Qtr4 into one column called "value"
dataku_long <- dataku %>% 
  pivot_longer(
    cols = starts_with("Qtr"),
    names_to = "quarter",
    values_to = "value"
  ) %>%
  arrange(tahun, quarter)   # Make sure data is ordered by year and quarter

# Remove missing values (e.g., 2025 Q3 and Q4 not yet available)
dataku_long <- dataku_long %>%
  filter(!is.na(value))

# Convert the cleaned data into a quarterly time series object
# start = first year and first quarter
# frequency = 4 because the data is quarterly
ts_data <- ts(
  dataku_long$value,
  start = c(min(dataku_long$tahun), 1),
  frequency = 4
)

# Plot the time series to visually inspect trend and seasonality
plot(ts_data)

# Run automatic SARIMA modeling
# The function will:
# - Search for the best SARIMA model
# - Apply diagnostic filtering
# - Rank candidate models
# - Select the best model
# - Produce forecasts and evaluation metrics
sarimasnap(ts_data)
```
## Outputs
<img width="659" height="749" alt="image" src="https://github.com/user-attachments/assets/c1f45c55-82b8-4eb5-bf38-c870c14be305" />
<img width="1257" height="784" alt="image" src="https://github.com/user-attachments/assets/39ccb002-3dfa-409f-b730-917060071d89" />

## Author

Joko Ade Nursiyono (East Java Data Analyst, BPS-Statistics Indonesia)

## Citation
If you use sarimasnap in research or official reports, please cite:
```r
Nursiyono, J. A. (2026). sarimasnap: Automatic Seasonal ARIMA Modeling with Diagnostic Filtering in R.
```


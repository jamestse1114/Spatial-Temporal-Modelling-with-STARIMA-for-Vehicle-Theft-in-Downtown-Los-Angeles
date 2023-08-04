# Spatial-Temporal-Modelling-with-STARIMA-for-Vehicle-Theft-in-Downtown-Los-Angeles

This repository contains code and methodologies used for spatio-temporal crime data analysis. The main focus is on weekly crime events data, performing an in-depth exploratory analysis, STL decomposition, spatial and temporal autocorrelation examination, and applying space-time autoregressive integrated moving average (STARIMA) modeling.

## Data

The analysis is based on the crime data from a large US city. The data includes time, location, and type of crime events that happened during a 10-year period. The original data was aggregated weekly for each neighborhood in the city.

## Methodology

1. **Data Cleaning & Preprocessing:** Performed necessary data preprocessing steps to clean the data and format it for analysis. Converted the data into a format compatible with spatial analysis tools in R.

2. **Exploratory Data Analysis (EDA):** Performed an in-depth EDA to understand the nature of the data. This includes examining the distribution of crime events and checking for patterns and outliers.

3. **Seasonal-Trend Decomposition (STL):** Used STL to decompose the time series data into three components: trend, seasonality, and residuals. The results were visualized to interpret the trend and seasonality in the data.

4. **Autocorrelation Examination:** Performed ACF and PACF on the STL residuals to examine the temporal autocorrelation in the data. 

5. **Spatial Autocorrelation Examination:** Calculated and interpreted Global and Local Moran's I to examine the spatial autocorrelation.

6. **STARIMA Modeling:** Applied STARIMA modeling to the data to capture the spatio-temporal autocorrelation. 

## Code

The project is entirely coded in R. You can find the code in the `.R` files in this repository. To run the code, you need to have R installed and it's recommended to use an IDE such as RStudio.

## Results

This analysis provided insights into the spatio-temporal patterns of crime in the city, showing the presence of both spatial and temporal autocorrelation. The results of the STARIMA model showed that this model is a good fit for the data, capturing the dependencies present in the data.

## Usage

To use this code, clone this GitHub repository and open the `.R` files in your R environment. Make sure to install all the necessary R packages using `install.packages()`.

## Contributing

If you'd like to contribute, please fork the repository and make changes as you'd like. Pull requests are warmly welcome.

## License

This project is open source, under the terms of the [MIT license](https://choosealicense.com/licenses/mit/).

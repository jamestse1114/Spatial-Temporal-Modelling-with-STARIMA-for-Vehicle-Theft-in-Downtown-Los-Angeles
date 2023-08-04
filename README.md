# Spatial-Temporal-Modelling-with-STARIMA-for-Vehicle-Theft-in-Downtown-Los-Angeles

This repository is dedicated to the research and analysis of vehicle theft in downtown Los Angeles using spatiotemporal data. The emphasis lies on weekly crime events data, including comprehensive exploratory analysis, STL decomposition, investigation of spatial and temporal autocorrelation, and application of the Space-Time AutoRegressive Integrated Moving Average (STARIMA) model.

## Dataset Description

The analysis utilizes crime data derived from a populous US city. This data incorporates the time, location, and nature of crime events documented over a decade. For ease of analysis, the raw data has been aggregated weekly for each neighbourhood within the city.

## Methodology

1. **Data Sanitation & Preprocessing:** We executed the necessary preprocessing steps to cleanse the data and prepare it for the subsequent analysis. The data was converted into a format that is compatible with spatial analysis tools provided in R.

2. **Exploratory Data Analysis (EDA):** A thorough EDA was carried out to understand the data's characteristics. This includes investigating the distribution of crime events and identifying patterns and potential outliers.

3. **Seasonal-Trend Decomposition (STL):** Utilizing STL, we divided the time series data into three components: trend, seasonality, and residuals. These results were then visualized to interpret the data's trend and seasonal fluctuations.

4. **Autocorrelation Examination:** The ACF and PACF were applied to the STL residuals to investigate the temporal autocorrelation present in the data.

5. **Spatial Autocorrelation Examination:** We computed and interpreted both Global and Local Moran's I to investigate spatial autocorrelation.
6. 
7. **STARIMA Modeling:** STARIMA modelling was applied to the data with the intention of capturing the spatiotemporal autocorrelation.

## Code

This project is entirely written in the R programming language. You can find the source code in the `.R` files within this repository. For executing this code, an installation of R is required and an Integrated Development Environment (IDE) such as RStudio is recommended.

## Findings

The analysis revealed insights into the spatio-temporal patterns of vehicle theft in the city, indicating the existence of both spatial and temporal autocorrelation. The STARIMA model's results demonstrated that this model suits the data well, accurately capturing the dependencies present in the data.

## Usage

To use this code, clone this GitHub repository and open the `.R` files in your R environment. Ensure to install all the necessary R packages with the `install.packages()` function.

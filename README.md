# Time Series Analysis and Forecasting of Sedimentary Deposits

## Description
This project is to forecast the thickness of the sedimentary deposits. The data set records the thickness of the sedimentary deposits each year over a period of 624 years at one particular location. By data transformation, model selection and model adequacy check, I built an ARIMA model to predict the thickness of the sedimentary deposits in R. Since different temperatures lead to different levels of melt, which, of course, affect the variation in the thickness of the sedimentary deposits. Thus, looking at these layers of silt and sediment can give an indication of the temperature changes from one year to the next.

## Data 
* Deposits: the thickness of the sedimentary deposits each year at one particular location


## Process, Method and Analysis
### Data Manipulation
* Missing Value check
* Visualize data by generating scatterplot over time 
* Transform data: do log transformation to Deposits
* Do Augmented Dickey-Fuller(adf) test to check whether the log data is stationary. The result is non-stationary.
* Difference the log data and do adf test again. Now the transformed data is stationary.
* Separate data into train and test set. The first 564 points as train data, the last 60 points as test data.


### Model Building
* Create plot of the sample autocorrelation(acf) and the sample partial autocorrelation(pacf) of the data
* Use ARIMA(p,1,q) model to fit the data with different p and q value and get AIC value for each model. Choose nine candidate 
  models with small AIC values.
* Review significance of estimated coefficients of all nine models. Only estimated coefficients of ARIMA(0, 1, 2), ARIMA(1,1,1)   are all significant.
* Compare the MSE of train data, AIC and MSE of test data. The final model is arima(1,1,1)
* Analyze the final model by checking the residual plots, qq plot, acf of residual and Ljung–Box test of residual. 
* Generate forecasting and its interval for next 12 points.
 

## Technologies
* R version: 3.5.2
* Package version
  * tseries version: 0.10-46
  * DMwR version: 3.5.3

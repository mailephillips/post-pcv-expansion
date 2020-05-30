# post-pcv-expansion

This project includes the code used for the manuscript "Evaluating post-vaccine expansion patterns of pneumococcal serotypes" (link). In this study, we evaluated patterns of change in serotype prevalence after pneumococcal conjugate vaccine introduction in Israel. For more details about the methods, please refer to the full manuscript.

## Files included in this project

### Data
*NOTE:* Data provided in this project are NOT the original data and were not the data used in the manuscript. These data have been simulated for the purposes of privacy and data sharing. 
- `carriage.csv`: contains aggregated carriage counts by serotype and epidemiological year
- `IPD.csv`: contains aggregated invasive pneumococaccal disease counts by serotype, epidemiological year, and age
- `all sts.csv`: contains a list of all serotypes used from either dataset, including vaccine-type serotypes

### R code
- `carriage model covmatx.R`: this code fits the carriage model and saves the output
- `cases gained.R`: this code estimates the number of cases gained due to increased carriage in children

## Running the code
To run this project, first run the code in `carriage model covmatx.R`. This takes several hours to run to allow for convergence of the covariance matrix. Next, run the code in `cases gained.R`. This code also takes several hours to run all the way through.

## Authors

* **Maile T Phillips**
* **Joshua L Warren**
* **Daniel M Weinberger**

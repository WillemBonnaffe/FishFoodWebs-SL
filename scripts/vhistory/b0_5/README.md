# Effects of temperature and DBO on lake and river trophic networks 

Author: Willem Bonnaffé (w.bonnaffe@gmail.com)

## Update log
* 21-06-2022 - uploaded repository on Github
* 05-07-2022 - created b0_5 
             - Corrected under-dispersion of missing DBO observations
             - introduced effect of richness in models
             - introduced interaction between dbo and temperature
             - removed average trophic level variable
             - introduced connectance as response variable
             - Cleaned repo structure 
* 07-07-2022 - Performed model selection 

## Next
* xx-xx-xxxx - Test robustness of results to spatial autocorrelation by modelling spatial covariation 

## Abstract

## Aims
* Estimate the effect of temperature and DBO on maximum and mean trophic level
* Compare effects between lake and river habitats
* Infer missing DBO observations
* Account for heterogeneity in hydrographic bassins

## Method

We use a hierarchical Bayesian modelling approach. The linear predictive model is given below:

$$ Yobs_{ij} \sim \mathcal{N}(\hat{Y}obs_{i},\Sigma_j) $$

$$ \hat{Y}obs_{i} = \beta_0 + year_i + temp_i + temp_i^2 + type_i + type_i \times temp_i + dbo_i + dbo_i^2 + type_i \times dbo_i + temp_i \times dbo_i + alt_i + rich_i $$

where 
$\beta_0$
is the intercept,
$year$
is the year since 201(?),
$temp$
is the temperature,
$type$
is the habitat,
$dbo$
the biochemical oxygen demand,
$alt$
the altitude,
$rich$
the species richness,
$\epsilon_i \sim \mathcal{N}(0,\Sigma_j)$
is a random effect of hydrographic bassin.

We use a slightly different model to account for missing DBO observations:

$$ Ymis_{ij} \sim \mathcal{N}(\hat{Y}mis_{i},\Sigma_j) $$

$$ \hat{Y}mis_{i} = \beta_0 + year_i + temp_i + temp_i^2 + type_i + type_i \times temp_i + x_{mis,i} + x_{mis,i}^2 + type_i \times x_{mis,i} + temp_i \times x_{mis,i} + alt_i + rich_i$$

$$ x_{mis,i} \sim \mathcal{N}(\mu_{mis},\sigma_{mis}) $$

where 
$x_{mis,i}$
are the missing DBO observations, which are assumed to be normally distributed. 
We use informative priors to constrain the distribution of missing values to be close to the distribution of observed DBO.

Combining these we obtain the full posterior distribution:

$$ p(\beta, \Sigma, \mu_{mis}, \sigma_{mis}| Yobs, Ymis) \propto ~ 
\prod_{i,j} ~
p(Yobs_{ij} | \beta, \Sigma_{j}) ~
p(Ymis_{ij} |\beta, \Sigma_{j}, x_{mis}) ~
p(x_{mis} | \mu_{mis}, \sigma_{mis}) ~
p(\beta) ~
p(\Sigma) ~
p(\mu_{mis}) ~ 
p(\sigma_{mis}) $$

where 
$\beta$ 
are the mean parameters, 
$\Sigma$ 
are the standard deviations of each hydrographic bassin,
$\mu_{mis}$ 
and 
$\sigma_{mis}$ 
are the mean and standard deviation of the missing DBO observations, 
$Yobs$ 
and 
$Ymis$ 
are the response variables (either max or mean trophic level) for the missing and observed datasets, 
$x_{mis}$ 
are the missing DBO observations.
$Yobs_{ij}$ 
indicates the 
$i^{th}$
observation of the 
$j^{th}$
hydrographic bassin.

## Results

For connectance, we find limited support for an effect of temperature as the interaction temperature by type, the linear and quadratic effect of temperature, as well as interaction with dbo were not found to be significant.
We find a linear effect of dbo, but not quadratic.
We find no evidence either for an effect of year on connectance.
For maximum trophic level, only the quadratic effect of temperature is not significant. 

![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_conne_model_0/fig_1.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_maxTL_model_0/fig_1.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_conne_model_0/fig_2.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_maxTL_model_0/fig_2.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_conne_model_0/fig_7.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_maxTL_model_0/fig_7.png)

Parameters mean estimates and confidence interval can be found here:
https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_conne_model_0/summary.csv
https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_maxTL_model_0/summary.csv

## Missing DBO distributions

The missing DBO distributions are slightly under-dispersed.

![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_conne_model_0/fig_3.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_maxTL_model_0/fig_3.png)

## Residuals

The under-dispersion of missing DBO residuals has been cured.

![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_conne_model_0/fig_4.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_conne_model_0/fig_5.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_maxTL_model_0/fig_4.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_maxTL_model_0/fig_5.png)

There is slight evidence for spatial autocorrelations in the residuals.

![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_conne_model_0/fig_18.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_maxTL_model_0/fig_18.png)

## Model selection

We performed model selection by removing non-significant effects from the models.
For the maximum trophic level, we removed the interaction temperature-dbo, the effect dbo^2, the interaction type-temperature, and the effect of the year.
For the connectance, we removed only the quadratic effect of dbo^2.

![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_conne_model_4/fig_1.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_maxTL_model_1/fig_1.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_conne_model_4/fig_2.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_maxTL_model_1/fig_2.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_conne_model_4/fig_7.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_maxTL_model_1/fig_7.png)

Parameters mean estimates and confidence interval can be found here:
https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_conne_model_4/summary.csv
https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_maxTL_model_1/summary.csv

## Missing DBO distributions

![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_conne_model_4/fig_3.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_maxTL_model_1/fig_3.png)

## Residuals

![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_conne_model_4/fig_4.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_conne_model_4/fig_5.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_maxTL_model_1/fig_4.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_maxTL_model_1/fig_5.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_conne_model_4/fig_18.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_5/out_maxTL_model_1/fig_18.png)


# Summary of analysis

Author: Willem Bonnaffé (w.bonnaffe@gmail.com)

## Update log
* 21-06-2022 - uploaded repository on Github

## Abstract

## Aims
* Estimate the effect of temperature and DBO on maximum and mean trophic level
* Compare effects between lake and river habitats

## Method

We use a hierarchical Bayesian modelling approach. The full posterior distribution is defined as follows:

$$ p(\beta, \Sigma, \mu_{mis}, \sigma_{mis}| Yobs, Ymis) \propto ~ 
\prod^{I,J} ~
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

The linear predictive model is given below:

$$ Y_{ij} = \beta_0 + year_i + temp_i + temp_i^2 + type_i + type_i \times temp_i + dbo_i + dbo_i^2 + type_i \times dbo_i + alt_i + \epsilon_j $$

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
$\epsilon_i \sim \mathcal{N}(0,\Sigma_j)$
is a random effect of hydrographic bassin.




## Results

![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_4/out_maxTL/fig_1.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_4/out_maxTL/fig_2.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_4/out_avgTL/fig_1.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_4/out_avgTL/fig_2.png)

Parameters mean estimates and confidence interval can be found here:
https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_4/out_maxTL/summary.csv
https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_4/out_avgTL/summary.csv

## Missing DBO distributions

![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_4/out_maxTL/fig_3.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_4/out_avgTL/fig_3.png)

## Residuals

![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_4/out_maxTL/fig_4.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_4/out_maxTL/fig_5.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_4/out_avgTL/fig_4.png)
![This is an image](https://github.com/WillemBonnaffe/RESOTRO/blob/main/riverlake/BM/b0_4/out_avgTL/fig_5.png)



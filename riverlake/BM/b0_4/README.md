# Summary of analysis

Author: Willem Bonnaffé (w.bonnaffe@gmail.com)

## Update log
* 21-06-2022 - uploaded repository on Github

## Abstract

## Aims
* Estimate the effect of temperature and DBO on maximum and mean trophic level
* Compare effects between lake and river habitats

## Method

We use a hierarchical Bayesian modelling approach.

$$ p(\beta, \Sigma | Y^{(obs)}) \propto p(Y^{(obs)}_{ij} | \beta, \Sigma_{j}) $$ 
p(Y^{(mis)}_{ij}|\beta, \Sigma_{j}, x^{(mis)}) p(x^{(mis)} | \mu^{(mis)}, \sd^{(mis)}) p(\beta) p(\Sigma) p(\mu^{(mis)}) p(\sd^{(mis)}) $$

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



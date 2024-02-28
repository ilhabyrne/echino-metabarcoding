# Model 1 - unformatted matrix
Call:
  rda(formula = pa ~ year + site + temperature + salinity + chl,      data = env1, scale = TRUE) 

Partitioning of correlations:
  Inertia Proportion
Total           7.000     1.0000
Constrained     1.769     0.2527
Unconstrained   5.231     0.7473

### Eigenvalues 
##### Unconstrained eigenvalyues
Importance of components:
  RDA1    RDA2    RDA3    RDA4     RDA5    PC1
Eigenvalue            0.9159 0.56032 0.16428 0.10395 0.024609 2.0399
Proportion Explained  0.1308 0.08005 0.02347 0.01485 0.003516 0.2914
Cumulative Proportion 0.1308 0.21089 0.23436 0.24921 0.252729 0.5441
PC2    PC3    PC4     PC5     PC6      PC7
Eigenvalue            1.0273 0.8727 0.5313 0.40869 0.30846 0.042564
Proportion Explained  0.1468 0.1247 0.0759 0.05838 0.04407 0.006081
Cumulative Proportion 0.6909 0.8156 0.8915 0.94985 0.99392 1.000000

##### Accumulated constrained eigenvalues
Importance of components:
  RDA1   RDA2    RDA3    RDA4    RDA5
Eigenvalue            0.9159 0.5603 0.16428 0.10395 0.02461
Proportion Explained  0.5177 0.3167 0.09286 0.05876 0.01391
Cumulative Proportion 0.5177 0.8345 0.92733 0.98609 1.00000

## Scaling 2 for species and site scores
* Species are scaled proportional to eigenvalues
* Sites are unscaled: weighted dispersion equal on all dimensions
* General scaling constant of scores:  3.395963 

##### Species scores
RDA1     RDA2     RDA3     RDA4     RDA5      PC1
Acanthaster -0.1370 -0.40918  0.08107 -0.11407  0.04222 -1.11136
Echinaster   0.6700  0.55667  0.06552 -0.19970  0.03737 -0.17601
Holothuria   0.4701 -0.41566  0.07065  0.10494 -0.03207  0.01794
Linckia     -0.2092 -0.22699 -0.34132 -0.18019  0.07745 -1.07470
Ophiactis    0.1097  0.18005 -0.33385  0.03472 -0.12091  0.53282
Ophionereis -0.7457  0.43230  0.03752  0.11414  0.03904 -0.42877
Ophiura     -0.4566 -0.04717  0.15946 -0.24629 -0.11914  0.68673

##### Site scores (weighted sums of species scores)
RDA1     RDA2    RDA3    RDA4     RDA5     PC1
Undi_2019 -0.1725  0.01209  1.3774  0.6231   0.1929  0.2332
Oste_2019 -0.4441 -0.46965 -1.0933 -1.4381   3.9357 -0.4020
Liza_2019 -0.4441 -0.46965 -1.0933 -1.4381   3.9357 -0.4087
Gree_2019  0.6233  0.03570  1.2768  3.1847  -3.4851  0.8718
Sudb_2019  0.1662 -1.35181 -0.5819 -0.2376   2.3859 -0.3974
Gibs_2019 -0.4441 -0.46965 -1.0933 -1.4381   3.9357 -0.4034
Hall_2019  0.1662 -1.35181 -0.5819 -0.2376   2.3859 -0.3879
Eddy_2019 -0.4441 -0.46965 -1.0933 -1.4381   3.9357 -0.4305
Otte_2019 -1.3402  0.68936  3.4003 -4.4470 -15.0771  1.7709
Bram_2019 -2.6542  1.62461 -0.4733  1.5424   8.2418 -1.1057
Bram_2017  1.4554  2.87705  1.5518 -1.8042   1.0591  0.1707
Gibs_2017  0.1662 -1.35181 -0.5819 -0.2376   2.3859 -0.6515
Tong_2017 -0.1725  0.01209  1.3774  0.6231   0.1929 -0.1791
Undi_2018  1.6087  0.60738  0.2046 -4.0259   5.3804 -0.7955
Gree_2018  0.6233  0.03570  1.2768  3.1847  -3.4851  0.7816
Sudb_2018  0.1662 -1.35181 -0.5819 -0.2376   2.3859 -0.6184
Eddy_2018  0.3380  1.79012 -4.7511  2.8909 -15.2719  1.3740
Yama_2018  0.6233  0.03570  1.2768  3.1847  -3.4851  0.6492
Brit_2018  0.6233  0.03570  1.2768  3.1847  -3.4851  0.6538
Fore_2018 -0.4441 -0.46965 -1.0933 -1.4381   3.9357 -0.7250


##### Site constraints (linear combinations of constraining variables)
RDA1      RDA2     RDA3     RDA4     RDA5     PC1
Undi_2019 -0.15440 -0.980236  1.05502  0.06187  0.81348  0.2332
Oste_2019 -0.65240 -0.502711 -0.31399 -0.91857 -0.43017 -0.4020
Liza_2019 -0.60663 -0.360201 -0.56432 -0.71938 -0.27734 -0.4087
Gree_2019 -0.07589 -0.221346  0.28750  0.90182  1.53504  0.8718
Sudb_2019 -0.39854 -0.699674  0.70879 -0.08000  0.53660 -0.3974
Gibs_2019  0.39137 -0.515082 -0.46748  0.78676  1.49289 -0.4034
Hall_2019  0.09024 -0.552025 -0.51215  0.06099  0.74988 -0.3879
Eddy_2019 -0.54017  0.093404 -1.38359 -0.55919 -0.14828 -0.4305
Otte_2019 -1.17746 -0.121636  0.41122 -0.63512 -0.30723  1.7709
Bram_2019 -1.92305  1.114787  0.09676  0.29435  0.10067 -1.1057
Bram_2017  0.84566  2.016615  0.25065  0.74192  0.13744  0.1707
Gibs_2017 -0.06859 -1.010686 -0.04125  2.00246 -1.79012 -0.6515
Tong_2017  0.61284  0.007152  1.41165  0.35499 -1.03111 -0.1791
Undi_2018  1.53256 -0.040616 -0.01808 -1.45078 -0.00480 -0.7955
Gree_2018  0.91661 -0.692155 -1.29936  0.03734 -0.53106  0.7816
Sudb_2018  0.78957 -0.093708  0.53703 -0.36235  0.07638 -0.6184
Eddy_2018  0.28283  0.464307 -0.86092  0.08954 -0.31179  1.3740
Yama_2018  0.01460  0.172435  1.32059 -0.90025 -0.55843  0.6492
Brit_2018 -0.03337  0.691017 -0.56719  0.40853 -0.20228  0.6538
Fore_2018  0.15425  1.230359 -0.05087 -0.11491  0.15022 -0.7250

##### Biplot scores for constraining variables
RDA1    RDA2    RDA3     RDA4     RDA5 PC1
year        -0.58356 -0.3407 -0.2088 -0.35410  0.61189   0
site         0.21859 -0.5260  0.6430 -0.50879 -0.05614   0
temperature  0.72818 -0.3911 -0.3758 -0.40791 -0.09602   0
salinity    -0.01864  0.4800 -0.1206 -0.20030 -0.84534   0
chl          0.67015  0.5760  0.2697  0.04363 -0.38012   0

## Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = pa ~ year + site + temperature + salinity + chl, data = env1, scale = TRUE)
Df Variance     F Pr(>F)
Model     5   1.7691 0.947   0.55
Residual 14   5.2309             

## ANOVAs
anova.cca(m2, by="axis", step=10000) # statistical significance of each axis
Permutation test for rda under reduced model
Forward tests for axes
Permutation: free
Number of permutations: 999

Model: rda(formula = pa ~ year + site + temperature + salinity + chl, data = env1, scale = TRUE)
Df Variance      F Pr(>F)
RDA1      1   0.9159 2.4514  0.489
RDA2      1   0.5603 1.4996  0.811
RDA3      1   0.1643 0.4397  0.998
RDA4      1   0.1040 0.2782  0.997
RDA5      1   0.0246 0.0659  0.997
Residual 14   5.2309              

anova.cca(m2, by="term", step=10000) # statistical significance of each term
Permutation test for rda under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

Model: rda(formula = pa ~ year + site + temperature + salinity + chl, data = env1, scale = TRUE)
Df Variance      F Pr(>F)
year         1   0.4064 1.0876  0.398
site         1   0.2991 0.8006  0.621
temperature  1   0.4675 1.2512  0.291
salinity     1   0.2210 0.5916  0.754
chl          1   0.3751 1.0038  0.399
Residual    14   5.2309              

vif.cca(m2) # anything above 10/20 should be avoided
year        site temperature    salinity         chl 
5.241758    1.195992    1.169185    1.565909    5.135178 

RsquareAdj(m2)$adj.r.squared
-0.01415366

# Model 2 - rare species removed

Call:
rda(formula = fin ~ year + site + temperature + salinity + chl,      data = env3, scale = TRUE) 

Partitioning of correlations:
              Inertia Proportion
Total           4.000     1.0000
Constrained     1.507     0.3768
Unconstrained   2.493     0.6232

## Eigenvalues
##### Unconstrained
Importance of components:
                        RDA1   RDA2    RDA3    PC1    PC2     PC3
Eigenvalue            0.9076 0.4987 0.10102 1.7298 0.4217 0.34117
Proportion Explained  0.2269 0.1247 0.02525 0.4325 0.1054 0.08529
Cumulative Proportion 0.2269 0.3516 0.37682 0.8093 0.9147 1.00000

##### Accumulated constrained eigenvalues
Importance of components:
                        RDA1   RDA2    RDA3
Eigenvalue            0.9076 0.4987 0.10102
Proportion Explained  0.6021 0.3309 0.06702
Cumulative Proportion 0.6021 0.9330 1.00000

## Scaling 2 for species and site scores
* Species are scaled proportional to eigenvalues
* Sites are unscaled: weighted dispersion equal on all dimensions
* General scaling constant of scores:  2.783158 


##### Species scores
               RDA1     RDA2      RDA3     PC1     PC2    PC3
Acanthaster -0.6037  0.06027 -0.237707  1.2106 -0.1637 0.1401
Echinaster   0.9925 -0.13055 -0.287283  0.2763  0.8248 0.3086
Holothuria   0.2085  0.97027 -0.009122 -0.5853 -0.2877 0.7254
Linckia     -0.6037  0.06027 -0.237707  1.2106 -0.1637 0.1401


##### Site scores (weighted sums of species scores)
             RDA1    RDA2    RDA3      PC1      PC2      PC3
Oste_2019 -0.7854 -0.9721 -0.9978  0.15168  0.62591 -0.86426
Liza_2019 -0.7854 -0.9721 -0.9978  0.17721  0.46417 -0.64785
Gree_2019  0.8813  0.7281  3.9878 -1.07467  0.38566  0.11310
Sudb_2019 -0.5538  0.9889 -1.0888  0.18001  0.13063  0.55756
Gibs_2019 -0.7854 -0.9721 -0.9978  0.73159  0.10628 -1.19928
Hall_2019 -0.5538  0.9889 -1.0888  0.15993 -0.42723  1.03823
Eddy_2019 -0.7854 -0.9721 -0.9978  0.14138 -0.04839  0.04854
Bram_2019 -0.7854 -0.9721 -0.9978  0.27251  0.17813  0.76075
Bram_2017  2.3031 -1.6286 -0.2207  0.05928  1.24379 -0.38261
Gibs_2017 -0.5538  0.9889 -1.0888  0.68037  0.17136  0.18939
Undi_2018  1.0995  0.5931 -5.3883  0.67165  0.93853  1.19908
Gree_2018  0.8813  0.7281  3.9878 -1.20071 -0.16628 -0.28327
Sudb_2018 -0.5538  0.9889 -1.0888  0.53332 -1.31805  0.03153
Yama_2018  0.8813  0.7281  3.9878 -1.13783 -0.18676 -0.72629
Brit_2018  0.8813  0.7281  3.9878 -1.13492 -0.66952  0.75941
Fore_2018 -0.7854 -0.9721 -0.9978  0.78919 -1.42822 -0.59401


##### Site constraints (linear combinations of constraining variables)
             RDA1     RDA2     RDA3      PC1      PC2      PC3
Oste_2019 -0.8712 -0.02776 -0.70864  0.15168  0.62591 -0.86426
Liza_2019 -0.7997 -0.23621 -0.57416  0.17721  0.46417 -0.64785
Gree_2019 -0.1478  0.30108  0.82168 -1.07467  0.38566  0.11310
Sudb_2019 -0.5823  0.73641  0.07484  0.18001  0.13063  0.55756
Gibs_2019 -0.1014  0.26460  0.52384  0.73159  0.10628 -1.19928
Hall_2019 -0.3725  0.15340 -0.04041  0.15993 -0.42723  1.03823
Eddy_2019 -0.6617 -0.95961 -0.52692  0.14138 -0.04839  0.04854
Bram_2019 -0.8041 -1.30386  0.67915  0.27251  0.17813  0.76075
Bram_2017  1.7333 -0.80733  0.65475  0.05928  1.24379 -0.38261
Gibs_2017 -0.1712  1.25028  1.46409  0.68037  0.17136  0.18939
Undi_2018  0.8926  0.46194 -1.41483  0.67165  0.93853  1.19908
Gree_2018  0.1236  0.28823 -0.36649 -1.20071 -0.16628 -0.28327
Sudb_2018  0.5555  0.66478 -0.34657  0.53332 -1.31805  0.03153
Yama_2018  0.2543  0.62156 -0.54073 -1.13783 -0.18676 -0.72629
Brit_2018  0.2754 -0.62703  0.31223 -1.13492 -0.66952  0.75941
Fore_2018  0.6772 -0.78048 -0.01182  0.78919 -1.42822 -0.59401


##### Biplot scores for constraining variables
                RDA1     RDA2     RDA3 PC1 PC2 PC3
year        -0.76184 -0.19552 -0.24128   0   0   0
site        -0.02123  0.68426 -0.62017   0   0   0
temperature  0.20102  0.43079 -0.68025   0   0   0
salinity     0.22350 -0.32157 -0.21393   0   0   0
chl          0.94215 -0.01841 -0.06963   0   0   0

## Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = fin ~ year + site + temperature + salinity + chl, data = env3, scale = TRUE)
         Df Variance      F Pr(>F)
Model     5   1.5073 1.2093  0.303
Residual 10   2.4927              

## ANOVAs
anova.cca(m3, by="axis", step=10000) # statistical significance of each axis
Permutation test for rda under reduced model
Forward tests for axes
Permutation: free
Number of permutations: 999

Model: rda(formula = fin ~ year + site + temperature + salinity + chl, data = env3, scale = TRUE)
         Df Variance      F Pr(>F)
RDA1      1  0.90756 4.3690  0.375
RDA2      1  0.49870 2.4008  0.709
RDA3      1  0.10102 0.4863  0.998
Residual 12  2.49272              

anova.cca(m3, by="term", step=10000) # statistical significance of each term
Permutation test for rda under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

Model: rda(formula = fin ~ year + site + temperature + salinity + chl, data = env3, scale = TRUE)
            Df Variance      F Pr(>F)  
year         1  0.55169 2.2132  0.095 .
site         1  0.28619 1.1481  0.379  
temperature  1  0.03646 0.1463  0.932  
salinity     1  0.17792 0.7138  0.561  
chl          1  0.45501 1.8253  0.202  
Residual    10  2.49272                

Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

vif.cca(m3) # anything above 10/20 should be avoided
       year        site temperature    salinity         chl 
   4.352530    1.429738    1.379851    1.470407    4.117496 

RsquareAdj(m3)$adj.r.squared
0.0652291

Call:
  rda(formula = pa ~ Year + Season + Temperature + Salinity, data = env1,      
      scale = TRUE) 

Partitioning of correlations:
  Inertia Proportion
Total          23.000     1.0000
Constrained     7.409     0.3221
Unconstrained  15.591     0.6779

Eigenvalues, and their contribution to the correlations 

Importance of components:
  RDA1   RDA2    RDA3    RDA4    PC1    PC2    PC3     PC4
Eigenvalue            3.5029 1.6423 1.32536 0.93817 5.1507 3.2170 2.6910 1.64875
Proportion Explained  0.1523 0.0714 0.05762 0.04079 0.2239 0.1399 0.1170 0.07168
Cumulative Proportion 0.1523 0.2237 0.28133 0.32212 0.5461 0.6859 0.8029 0.87462
PC5     PC6     PC7     PC8
Eigenvalue            1.06765 0.91214 0.61844 0.28557
Proportion Explained  0.04642 0.03966 0.02689 0.01242
Cumulative Proportion 0.92104 0.96070 0.98758 1.00000

Accumulated constrained eigenvalues
Importance of components:
  RDA1   RDA2   RDA3   RDA4
Eigenvalue            3.5029 1.6423 1.3254 0.9382
Proportion Explained  0.4728 0.2217 0.1789 0.1266
Cumulative Proportion 0.4728 0.6945 0.8734 1.0000

Scaling 2 for species and site scores
* Species are scaled proportional to eigenvalues
* Sites are unscaled: weighted dispersion equal on all dimensions
* General scaling constant of scores:  4.075935 


Species scores

RDA1      RDA2      RDA3     RDA4      PC1      PC2
Acanthaster     0.59785 -0.196551 -0.221648  0.07778  0.18469 -0.17428
Amphiodia       0.06959  0.397419 -0.332356 -0.41722 -0.17464  0.38224
Asteropsis     -0.27737  0.117307 -0.015303 -0.11276  0.26427  0.73103
Astropyga      -0.27737  0.117307 -0.015303 -0.11276  0.26427  0.73103
Dactylosaster  -0.41511  0.052678  0.137460  0.06566 -0.17913 -0.22824
Echinaster      0.49894 -0.080337  0.012220 -0.24460 -0.09243  0.04736
Eucidaris       0.09900 -0.003788 -0.048120  0.36729 -0.10839  0.06884
Fromia          0.35209  0.354503 -0.233562  0.14727 -0.23601  0.10290
Holothuria      0.51030  0.153318  0.348254 -0.04628  0.31900  0.24196
Koehleraster   -0.22261 -0.476004 -0.519399  0.05243 -0.05408 -0.08329
Labidodemas    -0.18308  0.260777 -0.105204 -0.02980  0.71947  0.27789
Linckia         0.55791  0.036165 -0.030241  0.01992  0.19833 -0.09431
Macrophiothrix  0.02948  0.235787 -0.127144  0.07240  0.70990 -0.35476
Mespilia        0.02948  0.235787 -0.127144  0.07240  0.70990 -0.35476
Ophiocoma      -0.18308  0.260777 -0.105204 -0.02980  0.71947  0.27789
Ophiomastix     0.02375  0.293882  0.081571  0.19971  0.46886 -0.33381
Ophionereis    -0.17595  0.257329 -0.044945  0.16107  0.45392 -0.29062
Ophiothrix      0.42531  0.108315  0.096315  0.31396  0.08493 -0.17302
Ophiura        -0.08137  0.150752  0.287020 -0.02694 -0.25431 -0.15373
Stichopus       0.65875  0.131335 -0.123692 -0.13826 -0.22064  0.10264
Synaptula       0.03782 -0.104679  0.384925 -0.19569 -0.05588 -0.07481
Thelenota      -0.28480  0.213045  0.007619  0.10197  0.39200 -0.43058
Tripneustes    -0.18308  0.260777 -0.105204 -0.02980  0.71947  0.27789


Site scores (weighted sums of species scores)

RDA1     RDA2     RDA3    RDA4     PC1     PC2
Dec_2015 -1.9815 -0.48748  0.53723  0.5117 -0.8254 -1.0517
Jan_2016 -1.9416  2.33556 -0.88799 -3.2269  1.2177  3.3683
Mar_2016 -0.3258 -0.68350  3.15517 -0.6330 -0.2575 -0.3447
Jun_2016  0.5852 -0.92359  0.50366  0.7822  0.6385  0.1810
Aug_2016 -0.9408 -2.89984 -2.41509  0.2673 -0.2492 -0.3838
Dec_2016 -0.2864 -0.20321  0.50813  1.7388 -0.4390 -0.1785
Jan_2017  0.1119 -0.03285  1.47511  1.4197 -0.3458 -0.4479
Mar_2017  1.2529 -1.02353  0.06991 -1.3585 -0.1617  0.2017
Dec_2017 -0.7824  0.03734  0.15322 -0.7177 -1.2494 -0.3273
Jan_2018 -0.3687  3.98150 -1.20888  3.0197  3.2710 -1.6346
Mar_2018  1.3101 -0.32613 -0.65278 -2.6402 -0.2745 -0.4646
Nov_2019  1.6810 -0.85605  0.12996  0.9776 -0.2372  0.6080
Jan_2020  1.6861  1.08178 -1.36765 -0.1407 -1.0875  0.4741


Site constraints (linear combinations of constraining variables)

RDA1    RDA2     RDA3    RDA4     PC1     PC2
Dec_2015 -1.9127  0.2427  0.63337  0.3025 -0.8254 -1.0517
Jan_2016 -1.2780  0.5405 -0.07051 -0.5196  1.2177  3.3683
Mar_2016  0.1743 -0.4823  1.77361 -0.9017 -0.2575 -0.3447
Jun_2016  0.6497 -2.3418  1.00687  0.2703  0.6385  0.1810
Aug_2016 -1.0257 -2.1933 -2.39323  0.2416 -0.2492 -0.3838
Dec_2016 -1.2336  0.5190  0.30544  0.6713 -0.4390 -0.1785
Jan_2017  0.0123  0.7470  1.09475  0.9124 -0.3458 -0.4479
Mar_2017  0.7933 -0.2068  1.42313 -0.1204 -0.1617  0.2017
Dec_2017 -0.7794  0.8336 -0.77731 -0.2069 -1.2494 -0.3273
Jan_2018  0.1358  1.0864 -0.58584  0.3336  3.2710 -1.6346
Mar_2018  0.9904  0.1642 -0.72845 -3.2818 -0.2745 -0.4646
Nov_2019  1.8512 -0.5426 -0.60565  1.6202 -0.2372  0.6080
Jan_2020  1.6223  1.6334 -1.07618  0.6786 -1.0875  0.4741


Biplot scores for constraining variables

RDA1    RDA2    RDA3     RDA4 PC1 PC2
Year         0.8130  0.4276 -0.3828  0.09820   0   0
Season       0.2807 -0.9424 -0.1797  0.02923   0   0
Temperature  0.3651  0.7540  0.5164 -0.17741   0   0
Salinity    -0.3582  0.1742 -0.1020  0.91158   0   0

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = pa ~ Year + Season + Temperature + Salinity, data = env1, scale = TRUE)
Df Variance      F Pr(>F)
Model     4   7.4087 0.9504  0.514
Residual  8  15.5913              
> anova.cca(m2, by="axis", step=10000) # statistical significance of each axis
Permutation test for rda under reduced model
Forward tests for axes
Permutation: free
Number of permutations: 999

Model: rda(formula = pa ~ Year + Season + Temperature + Salinity, data = env1, scale = TRUE)
Df Variance      F Pr(>F)
RDA1      1   3.5029 1.7974  0.403
RDA2      1   1.6423 0.8427  0.957
RDA3      1   1.3254 0.6801  0.934
RDA4      1   0.9382 0.4814  0.831
Residual  8  15.5913              

anova.cca(m2, by="term", step=10000) # statistical significance of each term
Permutation test for rda under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

Model: rda(formula = pa ~ Year + Season + Temperature + Salinity, data = env1, scale = TRUE)
Df Variance      F Pr(>F)
Year         1   2.8190 1.4464  0.106
Season       1   1.8756 0.9624  0.483
Temperature  1   1.7310 0.8882  0.562
Salinity     1   0.9832 0.5045  0.885
Residual     8  15.5913              

vif.cca(m2) # anything above 10/20 should be avoided
Year      Season Temperature    Salinity 
1.343434    3.144518    3.660236    1.521375 

RsquareAdj(m2)$adj.r.squared
-0.01682204

###############################################################################

Call:
  rda(formula = fin ~ Year + Season + Temperature + Salinity, data = env3,      
      scale = TRUE) 

Partitioning of correlations:
  Inertia Proportion
Total          15.000      1.000
Constrained     5.266      0.351
Unconstrained   9.734      0.649

Eigenvalues, and their contribution to the correlations 

Importance of components:
  RDA1    RDA2    RDA3    RDA4    PC1    PC2     PC3     PC4
Eigenvalue            2.9709 1.12926 0.72282 0.44266 3.5500 2.4719 1.30475 1.19193
Proportion Explained  0.1981 0.07528 0.04819 0.02951 0.2367 0.1648 0.08698 0.07946
Cumulative Proportion 0.1981 0.27334 0.32153 0.35104 0.5877 0.7525 0.83949 0.91895
PC5     PC6     PC7
Eigenvalue            0.59175 0.46402 0.16001
Proportion Explained  0.03945 0.03093 0.01067
Cumulative Proportion 0.95840 0.98933 1.00000

Accumulated constrained eigenvalues
Importance of components:
  RDA1   RDA2   RDA3    RDA4
Eigenvalue            2.9709 1.1293 0.7228 0.44266
Proportion Explained  0.5642 0.2145 0.1373 0.08407
Cumulative Proportion 0.5642 0.7787 0.9159 1.00000

Scaling 2 for species and site scores
* Species are scaled proportional to eigenvalues
* Sites are unscaled: weighted dispersion equal on all dimensions
* General scaling constant of scores:  3.584025 


Species scores

RDA1      RDA2      RDA3     RDA4      PC1      PC2
Acanthaster -0.74151 -0.009133  0.098748 -0.09941 -0.25526 -0.04911
Amphiodia    0.01785 -0.747452 -0.072282 -0.06998  0.21299  0.20112
Echinaster  -0.54720 -0.152790 -0.175959  0.10603  0.11032 -0.46776
Eucidaris   -0.06234  0.092974  0.471023  0.13354  0.13973 -0.52380
Holothuria  -0.53524  0.141285 -0.133214 -0.12735 -0.46133 -0.22123
Labidodemas  0.27074 -0.242207  0.106843 -0.06222 -0.77542  0.12136
Linckia     -0.57006 -0.146474  0.191960  0.11477 -0.23178 -0.52373
Ophiocoma    0.27074 -0.242207  0.106843 -0.06222 -0.77542  0.12136
Ophiomastix  0.04015  0.114932  0.107349 -0.41837 -0.52111  0.29611
Ophionereis  0.26695 -0.066753  0.230657 -0.13503 -0.49005 -0.33175
Ophiothrix  -0.40094  0.154298  0.353787 -0.09076 -0.01663  0.66881
Ophiura      0.12489  0.205262 -0.221214 -0.25444  0.33301  0.63650
Stichopus   -0.68385 -0.270266 -0.006912 -0.15940  0.21323 -0.40381
Thelenota    0.37735 -0.047135  0.153196 -0.05063 -0.39887 -0.01322
Tripneustes  0.27074 -0.242207  0.106843 -0.06222 -0.77542  0.12136


Site scores (weighted sums of species scores)

RDA1    RDA2    RDA3    RDA4      PC1     PC2
Dec_2015  1.62160  0.5261 -0.5290  1.9698  0.74060 -0.4532
Jan_2016  1.53449 -2.4735 -0.6245 -0.1900 -1.14218  0.2197
Mar_2016  0.50255  1.8059 -1.3979 -0.4992  0.23924  1.0449
Jun_2016 -0.60717  1.0240  0.3043  1.0635 -0.43136  0.7485
Dec_2016  0.47152  0.8022  1.7968  2.0277  0.28467 -2.0459
Jan_2017 -0.02293  2.0929 -0.6426 -3.8439  0.12936  1.0934
Mar_2017 -1.32263 -0.2028 -1.4569  1.2929  0.01139 -1.4737
Dec_2017  1.05729 -0.1393 -1.0822 -0.0111  1.56604  1.5594
Jan_2018  0.90416 -0.9035  3.6734 -4.2951 -2.73493  0.3871
Mar_2018 -1.30838 -1.7726 -1.6940  0.9180  0.21553 -0.6706
Nov_2019 -1.69164  0.3541  1.6081  1.7329  0.41396 -0.5731
Jan_2020 -1.13886 -1.1134  0.0445 -0.1654  0.70768  0.1635


Site constraints (linear combinations of constraining variables)

RDA1    RDA2    RDA3    RDA4      PC1     PC2
Dec_2015  1.78907  0.3759  0.1236  0.3885  0.74060 -0.4532
Jan_2016  1.25600 -0.5995 -0.1082  0.3306 -1.14218  0.2197
Mar_2016 -0.26460  1.0886 -1.7220 -0.2763  0.23924  1.0449
Jun_2016 -0.78030  1.0742  0.1601  2.4424 -0.43136  0.7485
Dec_2016  1.23707  0.2778  0.5109 -0.0335  0.28467 -2.0459
Jan_2017  0.10304  1.1862 -0.1056 -1.4502  0.12936  1.0934
Mar_2017 -0.74606  1.1491 -0.9965 -0.6879  0.01139 -1.4737
Dec_2017  0.88712 -1.0823  0.5425  0.2483  1.56604  1.5594
Jan_2018  0.09768 -0.6115  0.6424 -0.6417 -2.73493  0.3871
Mar_2018 -0.88411 -2.2238 -1.8700  0.6027  0.21553 -0.6706
Nov_2019 -1.54879  0.1871  1.8442  0.7012  0.41396 -0.5731
Jan_2020 -1.14612 -0.8217  0.9785 -1.6241  0.70768  0.1635


Biplot scores for constraining variables

RDA1    RDA2     RDA3    RDA4 PC1 PC2
Year        -0.7197 -0.4382  0.36767 -0.3934   0   0
Season      -0.6265  0.3085 -0.03585  0.7149   0   0
Temperature -0.2282 -0.1331 -0.38405 -0.8847   0   0
Salinity     0.4040  0.3806  0.79121 -0.2567   0   0

Model: rda(formula = fin ~ Year + Season + Temperature + Salinity, data = env3, scale = TRUE)
Df Variance      F Pr(>F)
Model     4   5.2656 0.9466  0.557
Residual  7   9.7344   

Permutation test for rda under reduced model
Forward tests for axes
Permutation: free
Number of permutations: 999

Model: rda(formula = fin ~ Year + Season + Temperature + Salinity, data = env3, scale = TRUE)
Df Variance      F Pr(>F)
RDA1      1   2.9709 2.1364  0.350
RDA2      1   1.1293 0.8120  0.970
RDA3      1   0.7228 0.5198  0.977
RDA4      1   0.4427 0.3183  0.943
Residual  7   9.7344              

anova.cca(m3, by="term", step=10000) # statistical significance of each term
Permutation test for rda under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

Model: rda(formula = fin ~ Year + Season + Temperature + Salinity, data = env3, scale = TRUE)
Df Variance      F Pr(>F)
Year         1   1.9221 1.3822  0.187
Season       1   1.4573 1.0479  0.389
Temperature  1   0.9434 0.6784  0.781
Salinity     1   0.9428 0.6780  0.734
Residual     7   9.7344              

vif.cca(m3) # anything above 10/20 should be avoided
Year      Season Temperature    Salinity 
1.459767    2.334003    2.634808    1.621219 

RsquareAdj(m3)$adj.r.squared
-0.01979395
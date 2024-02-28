# Model 1 - unformatted matrix
Call:
  rda(formula = pa ~ year + site + temperature + salinity + chl, data = env1, scale = TRUE) 

Partitioning of correlations:
| Category      | Inertia | Proportion |
|---------------|---------|------------|
| Total         | 7.000   | 1.0000     |
| Constrained   | 1.769   | 0.2527     |
| Unconstrained | 5.231   | 0.7473     |

### Eigenvalues 
#### Unconstrained eigenvalues
**Importance of components:**

| Metric                | RDA1    | RDA2    | RDA3    | RDA4    | RDA5    | PC1    | PC2    | PC3    | PC4    | PC5    | PC6    | PC7    |
|-----------------------|---------|---------|---------|---------|---------|--------|--------|--------|--------|--------|--------|--------|
| Eigenvalue            | 0.9159  | 0.56032 | 0.16428 | 0.10395 | 0.024609| 2.0399 | 1.0273 | 0.8727 | 0.5313 | 0.40869| 0.30846| 0.042564|
| Proportion Explained  | 0.1308  | 0.08005 | 0.02347 | 0.01485 | 0.003516| 0.2914 | 0.1468 | 0.1247 | 0.0759 | 0.05838| 0.04407| 0.006081|
| Cumulative Proportion | 0.1308  | 0.21089 | 0.23436 | 0.24921 | 0.252729| 0.5441 | 0.6909 | 0.8156 | 0.8915 | 0.94985| 0.99392| 1.000000|

#### Accumulated constrained eigenvalues
**Importance of components for RDAs 1 & 2:**

| Metric                | Value 1 | Value 2 |
|-----------------------|---------|---------|
| Eigenvalue            | 0.9159  | 0.5603  |
| Proportion Explained  | 0.5177  | 0.3167  |
| Cumulative Proportion | 0.5177  | 0.8345  |


## Scaling 2 for species and site scores
* Species are scaled proportional to eigenvalues
* Sites are unscaled: weighted dispersion equal on all dimensions
* General scaling constant of scores:  3.395963 

#### Species scores
| Species     | Value 1 | Value 2  |
|-------------|---------|----------|
| Acanthaster | -0.1370 | -0.40918 |
| Echinaster  |  0.6700 |  0.55667 |
| Holothuria  |  0.4701 | -0.41566 |
| Linckia     | -0.2092 | -0.22699 |
| Ophiactis   |  0.1097 |  0.18005 |
| Ophionereis | -0.7457 |  0.43230 |
| Ophiura     | -0.4566 | -0.04717 |


#### Biplot scores for constraining variables
| Variable    | Value 1 | Value 2 |
|-------------|---------|---------|
| year        | -0.58356| -0.3407 |
| site        |  0.21859| -0.5260 |
| temperature |  0.72818| -0.3911 |
| salinity    | -0.01864|  0.4800 |
| chl         |  0.67015|  0.5760 |

## Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = pa ~ year + site + temperature + salinity + chl, data = env1, scale = TRUE)
|          | Df | Variance |    F   | Pr(>F) |
|----------|----|----------|--------|--------|
| Model    | 5  | 1.7691   | 0.947  | 0.55   |
| Residual | 14 | 5.2309   |        |        |
     

## ANOVAs
anova.cca(m2, by="axis", step=10000) # statistical significance of each axis
Permutation test for rda under reduced model
Forward tests for axes
Permutation: free
Number of permutations: 999

Model: rda(formula = pa ~ year + site + temperature + salinity + chl, data = env1, scale = TRUE)
|          | Df | Variance |    F    | Pr(>F) |
|----------|----|----------|---------|--------|
| RDA1     | 1  | 0.9159   | 2.4514  | 0.489  |
| RDA2     | 1  | 0.5603   | 1.4996  | 0.811  |
| RDA3     | 1  | 0.1643   | 0.4397  | 0.998  |
| RDA4     | 1  | 0.1040   | 0.2782  | 0.997  |
| RDA5     | 1  | 0.0246   | 0.0659  | 0.997  |
| Residual | 14 | 5.2309   |         |        |

anova.cca(m2, by="term", step=10000) # statistical significance of each term
Permutation test for rda under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

Model: rda(formula = pa ~ year + site + temperature + salinity + chl, data = env1, scale = TRUE)
| Variable    | Df | Variance |    F    | Pr(>F) |
|-------------|----|----------|---------|--------|
| year        | 1  | 0.4064   | 1.0876  | 0.398  |
| site        | 1  | 0.2991   | 0.8006  | 0.621  |
| temperature | 1  | 0.4675   | 1.2512  | 0.291  |
| salinity    | 1  | 0.2210   | 0.5916  | 0.754  |
| chl         | 1  | 0.3751   | 1.0038  | 0.399  |
| Residual    | 14 | 5.2309   |         |        |


vif.cca(m2) # anything above 10/20 should be avoided
| Variable    | Value    |
|-------------|----------|
| year        | 5.241758 |
| site        | 1.195992 |
| temperature | 1.169185 |
| salinity    | 1.565909 |
| chl         | 5.135178 |


RsquareAdj(m2)$adj.r.squared
-0.01415366

# Model 2 - rare species removed

Call:
rda(formula = fin ~ year + site + temperature + salinity + chl,      data = env3, scale = TRUE) 

Partitioning of correlations:
| Category      | Inertia | Proportion |
|---------------|---------|------------|
| Total         | 4.000   | 1.0000     |
| Constrained   | 1.507   | 0.3768     |
| Unconstrained | 2.493   | 0.6232     |

## Eigenvalues
##### Unconstrained
Importance of components:
| Metric                | RDA1   | RDA2   | RDA3   | PC1    | PC2    | PC3    |
|-----------------------|--------|--------|--------|--------|--------|--------|
| Eigenvalue            | 0.9076 | 0.4987 | 0.10102| 1.7298 | 0.4217 | 0.34117|
| Proportion Explained  | 0.2269 | 0.1247 | 0.02525| 0.4325 | 0.1054 | 0.08529|
| Cumulative Proportion | 0.2269 | 0.3516 | 0.37682| 0.8093 | 0.9147 | 1.00000|

##### Accumulated constrained eigenvalues
Importance of components:
| Metric                | RDA1   | RDA2   | RDA3    |
|-----------------------|--------|--------|---------|
| Eigenvalue            | 0.9076 | 0.4987 | 0.10102 |
| Proportion Explained  | 0.6021 | 0.3309 | 0.06702 |
| Cumulative Proportion | 0.6021 | 0.9330 | 1.00000 |

## Scaling 2 for species and site scores
* Species are scaled proportional to eigenvalues
* Sites are unscaled: weighted dispersion equal on all dimensions
* General scaling constant of scores:  2.783158 

##### Species scores
|             | RDA1   | RDA2    | RDA3     | PC1    | PC2    | PC3   |
|-------------|--------|---------|----------|--------|--------|-------|
| Acanthaster | -0.6037| 0.06027 | -0.237707| 1.2106 | -0.1637| 0.1401|
| Echinaster  |  0.9925| -0.13055| -0.287283| 0.2763 |  0.8248| 0.3086|
| Holothuria  |  0.2085|  0.97027| -0.009122| -0.5853| -0.2877| 0.7254|
| Linckia     | -0.6037|  0.06027| -0.237707| 1.2106 | -0.1637| 0.1401|

##### Biplot scores for constraining variables
| Variable    | RDA1     | RDA2     | RDA3     | PC1 | PC2 | PC3 |
|-------------|----------|----------|----------|-----|-----|-----|
| year        | -0.76184 | -0.19552 | -0.24128 | 0   | 0   | 0   |
| site        | -0.02123 |  0.68426 | -0.62017 | 0   | 0   | 0   |
| temperature |  0.20102 |  0.43079 | -0.68025 | 0   | 0   | 0   |
| salinity    |  0.22350 | -0.32157 | -0.21393 | 0   | 0   | 0   |
| chl         |  0.94215 | -0.01841 | -0.06963 | 0   | 0   | 0   |

## Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = fin ~ year + site + temperature + salinity + chl, data = env3, scale = TRUE)
|          | Df | Variance |    F   | Pr(>F) |
|----------|----|----------|--------|--------|
| Model    | 5  | 1.5073   | 1.2093 | 0.303  |
| Residual | 10 | 2.4927   |        |        | 

## ANOVAs
anova.cca(m3, by="axis", step=10000) # statistical significance of each axis
Permutation test for rda under reduced model
Forward tests for axes
Permutation: free
Number of permutations: 999

Model: rda(formula = fin ~ year + site + temperature + salinity + chl, data = env3, scale = TRUE)
|          | Df | Variance |    F   | Pr(>F) |
|----------|----|----------|--------|--------|
| RDA1     | 1  | 0.90756  | 4.3690 | 0.375  |
| RDA2     | 1  | 0.49870  | 2.4008 | 0.709  |
| RDA3     | 1  | 0.10102  | 0.4863 | 0.998  |
| Residual | 12 | 2.49272  |        |        |

anova.cca(m3, by="term", step=10000) # statistical significance of each term
Permutation test for rda under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

Model: rda(formula = fin ~ year + site + temperature + salinity + chl, data = env3, scale = TRUE)
| Variable    | Df | Variance |    F   | Pr(>F) |
|-------------|----|----------|--------|--------|
| year        | 1  | 0.55169  | 2.2132 | 0.095  |
| site        | 1  | 0.28619  | 1.1481 | 0.379  |
| temperature | 1  | 0.03646  | 0.1463 | 0.932  |
| salinity    | 1  | 0.17792  | 0.7138 | 0.561  |
| chl         | 1  | 0.45501  | 1.8253 | 0.202  |
| Residual    | 10 | 2.49272  |        |        |

vif.cca(m3) # anything above 10/20 should be avoided
| Variable    | Value    |
|-------------|----------|
| year        | 4.352530 |
| site        | 1.429738 |
| temperature | 1.379851 |
| salinity    | 1.470407 |
| chl         | 4.117496 |

RsquareAdj(m3)$adj.r.squared
0.0652291

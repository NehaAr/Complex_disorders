# Complex_disorders
This repository contains different computational approaches to be used in studying complex disorders

![Comparison between different GWAS methods](/graphs/comparison_graph.png)

Comparison after 10 replicates

| Method      | Power Analysis                | FDR Handling                                | Characteristics                                                                                                       |
| ----------- | ---------------------------- | ------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| **GLM**     | Low to moderate (0.2 → 0.7)  | Poor at controlling FDR; power rises slowly | Simple regression, no population structure or kinship correction. Has high false positives                                        |
| **MLM**     | Moderate (0.2 → 0.45)        | Better than GLM, but still limited          | Mixed Linear Model; accounts for kinship, reduces false positives                                                             |
| **CMLM**    | Undefined  | ----                 | Compressed MLM; reduces computational load while modeling relatedness                                                         |
| **MLMM**    | High (0.6 → 0.95)            | Good FDR control                            | Multiple Locus Mixed Model; considers multiple SNPs simultaneously                                                            |
| **FarmCPU** | High (0.6 → 0.9)             | Good                                        | Fixed and random effect model; iterative, reduces confounding between testing marker and covariates                           |
| **BLINK**   | High (0.6 → 0.9)             | Good                                        | Bayesian-information and LD Iteratively Nested Keyway computational efficiency                                      |
| **SUPER**   | Moderate → high (0.2 → 0.95) | Very good                                   | Settlement of Multiple Effects Using Regression in SUPER GWAS when used with FAST-LLM, improves power for large-scale GWAS |


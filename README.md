# MELODY
Cell deconvolution for whole blood RNAseq samples



Installation
You can install the released version of MELODY from github with:

``` r
install.packages("devtools")
devtools::install_github("mackaay/MELODY")

library(e1071)
library(furrr)
```

Example
``` r
library(e1071)
library(furrr)
library(MELODY)
data(blood_signature)
data(example)
results <- MELODY(sig_matrix = blood_signature, mixture_file = example)
```

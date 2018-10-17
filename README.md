# Bayesian Model Optimized Reference Correction  (BaMORC)

# BaMORC <img src="man/figures/logo.png" height="20%" width="20%" align="right" />
The BaMORC package was desinged to faciliate Protein NMR researchers with a easy tool revoluting the traditional protein NMR research pipeline by provide several new methods to allow detect and correct 13C referencing error at early data analysis step. 

## Installation

The latest stable version can be installed from CRAN:

``` r
install.packages('BaMORC')
```

The latest development version can be installed from github:

``` r
# install.packages("devtools")
devtools::install_github('xxxx/BaMORC')
```

### Prerequisites
To reference correct assigned protein NMR spectra, following packages are required.

``` r
install.packages(c("data.table", "dplyr", "DEoptim", "httr", "docopt"))
```

To use unassigned protein NMR reference correction method, SSC is required and user need to run following code to get the R script

#### For Mac and Linux:

* Find the R CLI script location

Open terminal and type the following code:
```
> R -e 'system.file("exec", "bamorc.R", package = "BaMORC")'

```

You will see the R script location print out in the terminal as shown in following image.
<img src="man/figures/script_loc.png" height="68%" width="68%" align="center" />

And to test the R CLI script using the following pattern.
```
> <path to the R ClI scirpt>/bamorc.R -h
```

In the example code, it should be like:
```
/Users/bill/Library/R/3.5/library/BaMORC/exec/bamorc.R -h
```



* Install the SSC

```

```



## BaMORC Examples
```
library(BaMORC)

sequence <- paste(BaMORC::RefDB_data$carbonDat[[1]]$AA,collapse = "")
secondaryStructure <- paste(BaMORC::RefDB_data$carbonDat[[1]]$SS,collapse = "")
chemicaeclShifts_input <- BaMORC::RefDB_data$carbonDat[[1]][,c(4,5)]
from=-5
to=5

BaMORC(sequence, secondaryStructure, chemicaeclShifts_input, from=-5, to=5)
```
And the expected result is:

```
[1] 0.07073937 # Expected output

```



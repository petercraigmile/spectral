## Spectral

Email questions about the code to pfc &lt;AT&gt; stat.osu.edu

A R package for carrying out the spectral analysis of univariate time
series.  The package is contained in the `spectral` folder.  Different
versions of the R package are in the `releases` folder.

Now includes R functions for carrying out spectral-based tests of
periodicity.  Extra code used in Wei and Craigmile (2010) are
contained in the folder `spectral_tests_of_periodicity`.  Code testing
some of the function in the R package are in `tests`.

Supported by the National Science Foundation under award number DMS-0604963 and DMS-0906864.

### Installation

Please make sure you have a C and Fortran compiler installed for your version of R. Then, after installing the devtools R package, you can install this R package using

```
devtools::install_github("petercraigmile/spectral/spectral")
```


After installation, type `library(spectral)` to use the R library.


### Issues:

The following functions used for spectral-based tests of periodicity are undocumented.

```
 Deni.Wald.crit Deni.Wald.setup Deni.Wald.stat Thomson.test
 global.Ftest hcrits hearing.stat hearing.stat.crit
 hearing.stat.pvalue local.Ftest local.Ftest.Satterthwaite
````

The documentation for `Thomson.stat` is incorrect.




### References:

D. R. Brillinger. Time Series: Data Analysis and Theory. Holt, New York, NY, 1981.

D. Percival and A. Walden. Spectral Analysis for Univariate Time Series. Cambridge
University Press, Cambridge, 2020.

M. B. Priestley. Spectral Analysis and Time Series. (Vol. 1): Univariate Series. Academic
Press, London, UK, 1981.

A. T. Walden, E. J. McCoy and D. B. Percival (1995), The Effective Bandwidth of a Multitaper Spectral Estimator, Biometrika, 82, 201-14

L. Wei and P. F. Craigmile (2010). Global and local spectral-based
tests for periodicities. Biometrika, 97, 223-230.  

P. F. Craigmile, and W. M. King (2004). Periodogram based tests for
distortion products otoacoustic emission. The Journal of the
Acoustical Society of America, 116, 442-451

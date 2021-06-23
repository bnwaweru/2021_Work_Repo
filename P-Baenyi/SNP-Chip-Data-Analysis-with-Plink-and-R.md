SNP-Chip Data Analysis with Plink and R
================
Bernice Waweru
Fri 19, Mar 2021

-   [Session Info](#session-info)

### Session Info

Details of packages used with the work flow

``` r
devtools::session_info()
```

    ## - Session info ---------------------------------------------------------------
    ##  setting  value                       
    ##  version  R version 4.0.3 (2020-10-10)
    ##  os       Windows 10 x64              
    ##  system   x86_64, mingw32             
    ##  ui       RTerm                       
    ##  language (EN)                        
    ##  collate  English_United States.1252  
    ##  ctype    English_United States.1252  
    ##  tz       Africa/Nairobi              
    ##  date     2021-03-19                  
    ## 
    ## - Packages -------------------------------------------------------------------
    ##  package     * version date       lib source        
    ##  assertthat    0.2.1   2019-03-21 [2] CRAN (R 4.0.3)
    ##  callr         3.5.1   2020-10-13 [2] CRAN (R 4.0.3)
    ##  cli           2.2.0   2020-11-20 [2] CRAN (R 4.0.3)
    ##  crayon        1.3.4   2017-09-16 [2] CRAN (R 4.0.3)
    ##  desc          1.2.0   2018-05-01 [2] CRAN (R 4.0.3)
    ##  devtools      2.3.2   2020-09-18 [2] CRAN (R 4.0.3)
    ##  digest        0.6.27  2020-10-24 [1] CRAN (R 4.0.3)
    ##  ellipsis      0.3.1   2020-05-15 [2] CRAN (R 4.0.3)
    ##  evaluate      0.14    2019-05-28 [2] CRAN (R 4.0.3)
    ##  fansi         0.4.2   2021-01-15 [2] CRAN (R 4.0.3)
    ##  fs            1.5.0   2020-07-31 [2] CRAN (R 4.0.3)
    ##  glue          1.4.2   2020-08-27 [2] CRAN (R 4.0.3)
    ##  htmltools     0.5.1   2021-01-12 [2] CRAN (R 4.0.3)
    ##  knitr         1.30    2020-09-22 [2] CRAN (R 4.0.3)
    ##  lifecycle     0.2.0   2020-03-06 [2] CRAN (R 4.0.3)
    ##  magrittr      2.0.1   2020-11-17 [2] CRAN (R 4.0.3)
    ##  memoise       1.1.0   2017-04-21 [2] CRAN (R 4.0.3)
    ##  pkgbuild      1.2.0   2020-12-15 [2] CRAN (R 4.0.3)
    ##  pkgload       1.1.0   2020-05-29 [2] CRAN (R 4.0.3)
    ##  prettyunits   1.1.1   2020-01-24 [2] CRAN (R 4.0.3)
    ##  processx      3.4.5   2020-11-30 [2] CRAN (R 4.0.3)
    ##  ps            1.5.0   2020-12-05 [2] CRAN (R 4.0.3)
    ##  purrr         0.3.4   2020-04-17 [2] CRAN (R 4.0.3)
    ##  R6            2.5.0   2020-10-28 [2] CRAN (R 4.0.3)
    ##  remotes       2.2.0   2020-07-21 [2] CRAN (R 4.0.3)
    ##  rlang         0.4.10  2020-12-30 [2] CRAN (R 4.0.3)
    ##  rmarkdown     2.6     2020-12-14 [2] CRAN (R 4.0.3)
    ##  rprojroot     2.0.2   2020-11-15 [2] CRAN (R 4.0.3)
    ##  sessioninfo   1.1.1   2018-11-05 [2] CRAN (R 4.0.3)
    ##  stringi       1.5.3   2020-09-09 [2] CRAN (R 4.0.3)
    ##  stringr       1.4.0   2019-02-10 [2] CRAN (R 4.0.3)
    ##  testthat      3.0.1   2020-12-17 [2] CRAN (R 4.0.3)
    ##  usethis       2.0.0   2020-12-10 [2] CRAN (R 4.0.3)
    ##  withr         2.4.0   2021-01-16 [2] CRAN (R 4.0.3)
    ##  xfun          0.20    2021-01-06 [2] CRAN (R 4.0.3)
    ##  yaml          2.2.1   2020-02-01 [2] CRAN (R 4.0.3)
    ## 
    ## [1] C:/Users/BWaweru/OneDrive - CGIAR/Documents/R/win-library/4.0
    ## [2] C:/R/R-4.0.3/library

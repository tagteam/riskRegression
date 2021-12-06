```{=html}
<a href="http://cran.rstudio.com/web/packages/riskRegression/index.html"><img src="http://www.r-pkg.org/badges/version/riskRegression"></a>
<a href="http://cranlogs.r-pkg.org/downloads/total/last-month/riskRegression"><img src="http://cranlogs.r-pkg.org/badges/riskRegression"></a>
<a href="https://ci.appveyor.com/project/bozenne/riskRegression"><img src="https://ci.appveyor.com/api/projects/status/github/tagteam/riskRegression?svg=true" alt="Build status"></a>
<a href="https://github.com/tagteam/riskRegression/actions"><img src="https://github.com/tagteam/riskRegression/workflows/R-CMD-check/badge.svg" alt="R build status"></a>
```
# Install

``` {.r org-language="R" exports="both" eval="never"}
library(devtools)
install_github("tagteam/riskRegression")
```

# References

The following references provide the methodological framework for the
features of riskRegression.

1.  T.A. Gerds and M. Schumacher. Consistent estimation of the expected
    Brier score in general survival models with right-censored event
    times. Biometrical Journal, 48(6):1029--1040, 2006.

2.  T.A. Gerds and M. Schumacher. Efron-type measures of prediction
    error for survival analysis. Biometrics, 63(4):1283--1287, 2007.

3.  T.A. Gerds, T. Cai, and M. Schumacher. The performance of risk
    prediction models. Biometrical Journal, 50(4):457--479, 2008.

4.  U B Mogensen, H. Ishwaran, and T A Gerds. Evaluating random forests
    for survival analysis using prediction error curves. Journal of
    Statistical Software, 50(11), 2012.

5.  P. Blanche, J-F Dartigues, and H. Jacqmin-Gadda. Estimating and
    comparing time-dependent areas under receiver operating
    characteristic curves for censored event times with competing risks.
    Statistics in Medicine, 32(30): 5381--5397, 2013.

6.  Paul Blanche, Ce\'cile Proust-Lima, Lucie Loube\`re, Claudine Berr,
    Jean- Franc,ois Dartigues, and He\'le\`ne Jacqmin-Gadda. Quantifying
    and comparing dynamic predictive accuracy of joint models for
    longitudinal marker and time-to-event in presence of censoring and
    competing risks. Biometrics, 71 (1):102--113, 2015.

Functions `predict.CauseSpecificCox`{.verbatim}, `predictCox`{.verbatim}
and `iidCox`{.verbatim}:

-   Brice Ozenne, Anne Lyngholm Sorensen, Thomas Scheike, Christian
    Torp-Pedersen and Thomas Alexander Gerds. riskRegression: Predicting
    the Risk of an Event using Cox Regression Models. The R
    Journal (2017) 9:2, pages 440-460.

```{=latex}
@article{gerds2006consistent,
  title =    {Consistent Estimation of the Expected {B}rier Score
                  in General Survival Models with Right-Censored Event
                  Times},
  author =   {Gerds, T.A. and Schumacher, M.},
  journal =  {Biometrical Journal},
  volume =   48,
  number =   6,
  pages =    {1029--1040},
  year =     2006,
  publisher =    {Wiley Online Library}
}

@article{gerds2007efron,
  title =    {Efron-Type Measures of Prediction Error for Survival
                  Analysis},
  author =   {Gerds, T.A. and Schumacher, M.},
  journal =  {Biometrics},
  volume =   63,
  number =   4,
  pages =    {1283--1287},
  year =     2007,
  publisher =    {Wiley Online Library}
}

@article{gerds2008performance,
  title =    {The performance of risk prediction models},
  author =   {Gerds, T.A. and Cai, T. and Schumacher, M.},
  journal =  {Biometrical Journal},
  volume =   50,
  number =   4,
  pages =    {457--479},
  year =     2008,
  publisher =    {Wiley Online Library}
}

@Article{mogensen2012pec,
  title =    {Evaluating random forests for survival analysis
                  using prediction error curves},
  author =   {Mogensen, U B and Ishwaran, H. and Gerds, T A},
  journal =  {Journal of Statistical Software},
  year =     2012,
  volume =   50,
  number =   11
}

@article{Blanche2013statmed,
  title =    "{Estimating and comparing time-dependent areas under
                  receiver operating characteristic curves for
                  censored event times with competing risks}",
  author =   {Blanche, P. and Dartigues, J-F and Jacqmin-Gadda,
                  H.},
  journal =  {Statistics in Medicine},
  volume =   32,
  number =   30,
  pages =    {5381--5397},
  year =     2013
}

@article{blanche2015,
  title =    {Quantifying and comparing dynamic predictive
                  accuracy of joint models for longitudinal marker and
                  time-to-event in presence of censoring and competing
                  risks},
  author =   {Blanche, Paul and Proust-Lima, C{\'e}cile and
                  Loub{\`e}re, Lucie and Berr, Claudine and Dartigues,
                  Jean-Fran{\c{c}}ois and Jacqmin-Gadda,
                  H{\'e}l{\`e}ne},
  journal =  {Biometrics},
  volume =   71,
  number =   1,
  pages =    {102--113},
  year =     2015,
  publisher =    {Wiley Online Library}
}

@article{ozenne2017,
  title =    {riskRegression: Predicting the Risk of an Event
                using Cox Regression Modelss},
  author =   {Ozenne, Brice and Sørensen, Anne Lyngholm 
                and Scheike, Thomas and Torp-Pedersen, Christian
                and Gerds, Thomas Alexander},
  journal =  {The R Journal},
  volume =   9,
  number =   2,
  pages =    {440--460},
  year =     2017
}
```

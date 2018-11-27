ifelse(!require (betareg), install.packages("betareg"), library(betareg))
ifelse(!require (tidyverse), install.packages("tidyverse"), library(tidyverse))
ifelse(!require (nlme), install.packages("nlme"), library(nlme))
ifelse(!require (lavaan), install.packages("lavaan"), library(lavaan))
ifelse(!require (lavaan.survey), install.packages("lavaan.survey"), library(lavaan.survey))
ifelse(!require (piecewiseSEM), install.packages("piecewiseSEM"), library(piecewiseSEM))
ifelse(!require (brms), install.packages("brms"), library(brms))


### Read in data
plot_dat <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/Grace2016/master/Grace2016_Plot_Data.csv",
                       sep = ",",
                       header = TRUE)
### rename colnames as "." can not be accepted by the piecewiseSEM functions
colnames(plot_dat) <- 
  c("psitecode", "prich", "ln_prich", "site_rich", "ln_site_rich", "ptotmass", "ln_ptotmass", "ln_site_totmass",
    "ln_pshade", "pprod", "ln_pprod", "ln_site_prod", "SoilSuitability", "SoilWithShade", "pn", "sqr_pn", "ln_pn", 
    "ln_pc", "ln_pp")

### Reproducing Grace et al. 2016's results
### The following code are identical from the suppliment information of Grace et al. 2016
plotmod_1 <- 'ln_prich ~ ln_site_rich + ln_pshade + SoilSuitability 
              ln_pshade   ~ ln_ptotmass + SoilWithShade 
              ln_ptotmass ~ ln_site_totmass + ln_pprod  
              ln_pprod ~  ln_site_prod + ln_prich                                     
              ln_pprod ~~ ln_ptotmass #controlling for error correlation
             '
survey_designplot <- svydesign(ids = ~psitecode, nest = TRUE, data = plot_dat)
survey_fit.plot <- lavaan.survey(lavaan.fit = plotmod_fit_1, survey.design = survey_designplot)
param_Grace <- parameterEstimates(survey_fit.plot)[1:9, c("lhs", "rhs", "est", "pvalue")]

### Plain SEM
plotmod_1 <- 'ln_prich ~ ln_site_rich + ln_pshade + SoilSuitability 
              ln_pshade   ~ ln_ptotmass + SoilWithShade 
              ln_ptotmass ~ ln_site_totmass + ln_pprod  
              ln_pprod ~  ln_site_prod + ln_prich                                     
              ln_pprod ~~ ln_ptotmass #controlling for error correlation
             '
plotmod_fit_1 <- sem(plotmod_1, data = plot_dat, meanstructure = TRUE)
param0 <- parameterEstimates(plotmod_fit_1)[1:9, c("lhs", "rhs", "est", "pvalue")]

### hierarchical SEM (site variables are considered)  with piecewiseSEM package
mod_Site <- 
  list(
    lme(ln_prich ~ ln_pshade + SoilSuitability + ln_site_rich, 
        random = ~ 1 | psitecode, 
        data = plot_dat, 
        control =  lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000))),
    lme(ln_pshade ~ ln_ptotmass + SoilWithShade, 
        random = ~ 1 | psitecode, 
        data = plot_dat, control = lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000))),
    lme(ln_ptotmass ~ ln_pprod + ln_site_totmass, 
        random = ~ 1 | psitecode, 
        data = plot_dat, control = lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000))),
    lme(ln_pprod ~ ln_prich + ln_site_prod, 
        random = ~ 1 | psitecode, 
        data = plot_dat, control = lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000)))
  )  

sem.coefs(modelList = mod_Site, data = plot_dat, corr.errors = "ln_pprod ~~ ln_ptotmass")
sem.fit(modelList = mod_Site, data = plot_dat, corr.errors = "ln_pprod ~~ ln_ptotmass",
        model.control = list(lmeControl(maxIter = 1000, msMaxIter = 1000, opt = "optim")),
        conditional = TRUE)
sem.model.fits(modelList = mod_Site)

### hierarchical SEM (site variables are NOT considered)  with piecewiseSEM package
mod_NoSite <- 
  list(
    lme(ln_prich ~ ln_pshade + SoilSuitability, 
        random = ~ 1 | psitecode, 
        data = plot_dat, 
        control =  lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000))),
    lme(ln_pshade ~ ln_ptotmass + SoilWithShade, 
        random = ~ 1 | psitecode, 
        data = plot_dat, control = lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000))),
    lme(ln_ptotmass ~  ln_pprod, 
        random = ~ 1 | psitecode, 
        data = plot_dat, control = lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000))),
    lme(ln_pprod ~ ln_prich, 
        random = ~ 1 | psitecode, 
        data = plot_dat, control = lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000)))
  ) 

sem.coefs(modelList = mod_NoSite, data = plot_dat, corr.errors = "ln_pprod ~~ ln_ptotmass")
sem.fit(modelList = mod_NoSite, data = plot_dat, corr.errors = "ln_pprod ~~ ln_ptotmass",
        model.control = list(lmeControl(maxIter = 1000, msMaxIter = 1000, opt = "optim")),
        conditional = TRUE)
sem.model.fits(modelList = mod_NoSite)

### Additional analyses, using brms package to do hierarchical SEM with Bayesian method
### (site varaiable considered)
bf_prich_Site <- brmsformula(ln_prich ~ ln_pshade + SoilSuitability + ln_site_rich + (1|ID|psitecode))
bf_pshade_Site <- brmsformula(ln_pshade ~ ln_ptotmass + SoilWithShade + (1|ID|psitecode))
bf_ptotmass_Site <- brmsformula(ln_ptotmass ~ ln_pprod + ln_site_totmass + (1|ID|psitecode))
bf_pprod_Site <- brmsformula(ln_pprod ~ ln_prich + ln_site_prod + (1|ID|psitecode))
bform_Site <- bf_prich_Site + bf_pshade_Site + bf_ptotmass_Site + bf_pprod_Site + set_rescor(FALSE)
brms_fit_Site <- brm(bform_Site, data = plot_dat, cores = 4, refresh = 0, chains = 4)
summary(brms_fit_Site, waic = TRUE)

### Additional analyses, using brms package to do hierarchical SEM with Bayesian method
### (site varaiable NOT considered)
bf_prich_NoSite <- brmsformula(ln_prich ~ ln_pshade + SoilSuitability + (1|ID|psitecode))
bf_pshade_NoSite <- brmsformula(ln_pshade ~ ln_ptotmass + SoilWithShade + (1|ID|psitecode))
bf_ptotmass_NoSite <- brmsformula(ln_ptotmass ~  ln_pprod + (1|ID|psitecode))
bf_pprod_NoSite <- brmsformula(ln_pprod ~ ln_prich + (1|ID|psitecode))
bform_NoSite <- bf_prich_NoSite + bf_pshade_NoSite + bf_ptotmass_NoSite + bf_pprod_NoSite + set_rescor(FALSE)
brms_fit_NoSite <- brm(bform_NoSite, data = plot_dat, cores = 4, refresh = 0, chains = 4)
summary(brms_fit_NoSite, waic = TRUE)



cor(plot_dat[,"ln_prich"], plot_dat[, "ln_site_rich"])








library(tidyverse)
library(nlme)
library(lavaan)
library(lavaan.survey)
library(piecewiseSEM)
library(brms)

plot_dat <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/Grace2016/master/Grace2016_Plot_Data.csv",
                       sep = ",",
                       header = TRUE)
colnames(plot_dat) <- 
  c("psitecode", "prich", "ln_prich", "site_rich", "ln_site_rich", "ptotmass", "ln_ptotmass", "ln_site_totmass",
    "ln_pshade", "pprod", "ln_pprod", "ln_site_prod", "SoilSuitability", "SoilWithShade", "pn", "sqr_pn", "ln_pn", 
    "ln_pc", "ln_pp")
### Using Grace's code to estimate effects of each variable ################

plotmod_1 <- 'ln_prich ~ ln_site_rich + ln_pshade + SoilSuitability 
              ln_pshade   ~ ln_ptotmass + SoilWithShade 
              ln_ptotmass ~ ln_site_totmass + ln_pprod  
              ln_pprod ~  ln_site_prod + ln_prich                                     
              ln_pprod ~~ ln_ptotmass #controlling for error correlation
             '

plotmod_fit_1 <- sem(plotmod_1, data = plot_dat, meanstructure = TRUE)
summary(plotmod_fit_1)
# adjusting for nested design
survey_designplot <- svydesign(ids = ~psitecode, nest = TRUE, data = plot_dat)
survey_fit.plot <- lavaan.survey(lavaan.fit = plotmod_fit_1, survey.design = survey_designplot)
survey_fit.plot  # gives chi-square
summary(survey_fit_plot)
AIC(survey_fit_plot)
standardizedSolution(survey_fit_plot)
### Using Grace's code to estimate effects of each variable ################

### Using piecewiseSEM (with random effects) to estimate effects of each variable ################

mod_list_1 <- 
  list(
    lme(ln_prich ~ ln_pshade + SoilSuitability, 
        random = ~ 1 | psitecode, 
        data = plot_dat, 
        control =  lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000))),
    lme(ln_pshade ~ ln_ptotmass + SoilWithShade, 
        random = ~ 1 | psitecode, 
        data = plot_dat, control = lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000))),
    lme(ln_ptotmass ~ ln_pprod, 
        random = ~ 1 | psitecode, 
        data = plot_dat, control = lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000))),
    lme(ln_pprod ~ ln_prich, 
        random = ~ 1 | psitecode, 
        data = plot_dat, control = lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000)))
  ) 

sem.coefs(modelList = mod_list_1, data = plot_dat, corr.errors = "ln_pprod ~~ ln_ptotmass")
sem.fit(modelList = mod_list_1, data = plot_dat, corr.errors = "ln_pprod ~~ ln_ptotmass",
        model.control = list(lmeControl(maxIter = 1000, msMaxIter = 1000, opt = "optim")),
        conditional = TRUE)
sem.model.fits(modelList = mod_list_1)
sem.plot(modelList = mod_list_1, data = plot_dat, corr.errors = "ln_pprod ~~ ln_ptotmass")
sem.lavaan(modelList = mod_list_1, data = plot_dat, corr.errors = "ln_pprod ~~ ln_ptotmass")


mod_list <- 
list(
  lme(ln_prich ~ ln_pshade + SoilSuitability, 
      random = ~ ln_pshade | psitecode, 
      data = plot_dat, 
      control =  lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000))),
  lme(ln_pshade ~ ln_ptotmass + SoilWithShade, 
      random = ~ 1 + ln_ptotmass + SoilWithShade | psitecode, 
      data = plot_dat, control = lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000))),
  lme(ln_ptotmass ~ ln_pprod, 
      random = ~ 1 + ln_pprod | psitecode, 
      data = plot_dat, control = lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000))),
  lme(ln_pprod ~ ln_prich, 
      random = ~ 1 + ln_prich | psitecode, 
      data = plot_dat, control = lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000)))
) 

sem.coefs(modelList = mod_list, data = plot_dat, corr.errors = "ln_pprod ~~ ln_ptotmass")
sem.fit(modelList = mod_list, data = plot_dat, corr.errors = "ln_pprod ~~ ln_ptotmass",
        model.control = list(lmeControl(maxIter = 1000, msMaxIter = 1000, opt = "optim")),
        conditional = TRUE)
sem.model.fits(modelList = mod_list)
sem.plot(modelList = mod_list, data = plot_dat, corr.errors = "ln_pprod ~~ ln_ptotmass")
sem.lavaan(modelList = mod_list, data = plot_dat, corr.errors = "ln_pprod ~~ ln_ptotmass")



bf_prich <- brmsformula(ln_prich ~ ln_pshade + SoilSuitability + (1|ID|psitecode))
bf_pshade <- brmsformula(ln_pshade ~ ln_ptotmass + SoilWithShade + (1|ID|psitecode))
bf_ptotmass <- brmsformula(ln_ptotmass ~ ln_pprod + (1|ID|psitecode))
bf_pprod <- brmsformula(ln_pprod ~ ln_prich + (1|ID|psitecode))
bform <- bf_prich + bf_pshade + bf_ptotmass + bf_pprod + set_rescor(FALSE)
brms_fit <- brm(bform, data = plot_dat, cores = 4, refresh = 0, chains = 4)
summary(brms_fit, waic = TRUE)


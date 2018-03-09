library(tidyverse)
library(nlme)
library(lavaan)
library(lavaan.survey)
library(piecewiseSEM)
library(brms)

plot.dat <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/Grace2016/master/Grace2016_Plot_Data.csv",
                       sep = ",",
                       header = TRUE)

### Using Grace's code to estimate effects of each variable ################

plotmod_1 <- 'ln.prich ~ ln.site.rich + ln.pshade + SoilSuitability 
              ln.pshade   ~ ln.ptotmass + SoilWithShade 
              ln.ptotmass ~ ln.site.totmass + ln.pprod  
              ln.pprod ~  ln.site.prod + ln.prich                                     
              ln.pprod ~~ ln.ptotmass #controlling for error correlation
             '

plotmod_fit_1 <- sem(plotmod_1, data = plot.dat, meanstructure = TRUE)
summary(plotmod_fit_1)
# adjusting for nested design
survey.designplot <- svydesign(ids = ~psitecode, nest = TRUE, data = plot.dat)
survey.fit.plot <- lavaan.survey(lavaan.fit = plotmod_fit_1, survey.design = survey.designplot)
survey.fit.plot  # gives chi-square
summary(survey.fit.plot)
AIC(survey.fit.plot)
standardizedSolution(survey.fit.plot)
### Using Grace's code to estimate effects of each variable ################

### Using piecewiseSEM (with random effects) to estimate effects of each variable ################

mod.list <- 
list(
  lme(ln.prich ~ ln.pshade + SoilSuitability, 
      random = ~ ln.pshade | psitecode, 
      data = plot.dat, 
      control =  lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000))),
  lme(ln.pshade ~ ln.ptotmass + SoilWithShade, 
      random = ~ 1 + ln.ptotmass + SoilWithShade | psitecode, 
      data = plot.dat, control = lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000))),
  lme(ln.ptotmass ~ ln.pprod, 
      random = ~ 1 + ln.pprod | psitecode, 
      data = plot.dat, control = lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000))),
  lme(ln.pprod ~ ln.prich, 
      random = ~ 1 + ln.prich | psitecode, 
      data = plot.dat, control = lmeControl(lmeControl(maxIter = 1000, msMaxIter = 1000)))
) 

sem.coefs(modelList = mod.list, data = plot.dat, corr.errors = "ln.pprod ~~ ln.ptotmass")
sem.fit(modelList = mod.list, data = plot.dat, corr.errors = "ln.pprod ~~ ln.ptotmass",
        model.control = list(lmeControl(maxIter = 1000, msMaxIter = 1000, opt = "optim")),
        conditional = TRUE)
sem.model.fits(modelList = mod.list)
sem.plot(modelList = mod.list, data = plot.dat, corr.errors = "ln.pprod ~~ ln.ptotmass")
sem.lavaan(modelList = mod.list, data = plot.dat, corr.errors = "ln.pprod ~~ ln.ptotmass")




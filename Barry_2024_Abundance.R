#######################
#  Abundance GLMM's for 
#  both GLSA and TATO
#   
#   1/10/2024 BRB
########################
library(dplyr)
rm(list=ls())
# read in the data files
Ncovs <- read.csv("./doc/final_files/Ncovs_master_final.csv")

#####################
# scale optimization
####################

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}


model.scores <- NULL
p.values <- NULL
estimate <- NULL
dispersion <- NULL

for (i in 16:(ncol(Ncovs))){

  # Fit with either abundance of each species N_neto or N_glor
  fm <- lme4::glmer.nb(N_neto~Ncovs[,i]+ (1|YrSt) +(1|Stand)+Study, data = Ncovs) 
  print(round(overdisp_fun(fm),2))
  mod <- summary(fm)

  model.scores <-rbind(model.scores,mod$AICtab[1])        
  p.values <- rbind(p.values,mod$coefficients[2,4])        
  estimate <- rbind(estimate,mod$coefficients[2,1])
  dispersion <- rbind(dispersion,round(overdisp_fun(fm),2))
}

scale <- c("1","9","11","13","15","17","19","21","23","25","27","29","31","33","35","37","39")
params <- c("asp","elev","tpi","slop","cvr","fhd","pa20","pa40","pa10","rh50","rh75","rh98")
scale.output <- as.data.frame(cbind(model.scores[1:(length(model.scores)-4)],p.values[1:(length(p.values)-4)],estimate[1:(length(estimate)-4)]))

scale.output <- cbind(scale.output,scale=factor(rep(scale,length(params)), levels=scale),
                      params <- factor(rep(params,each=length(scale)),levels=params))

colnames(scale.output) <- c("AIC","p.value","Coefficient","Scale","Params")

# plot and examine for best parameters for each species


##############
# GLOR Models     
##############
# glor optim variables. Check for corrilations and those above 0.7 can't be in the same model
round(cor(Ncovs[,c("asp_13", "elev_39", "tpi_11", "slop_13", "cvr_13", "fhd_17", "pa20_27", "pa40_39", "pa10_13","rh50_27","rh75_15","rh98_15")]),2)

# Null
fm.null<- lme4::glmer.nb(N_glor~(1|YrSt)+(1|Stand)+Study, data = Ncovs)

# Abiotic
fm.a.elev.tmax.asp <- lme4::glmer.nb(N_glor~elev_39+tmax+asp_13+(1|YrSt)+(1|Stand)+Study, data = Ncovs) 
fm.a.elev.tmin.ppt <- lme4::glmer.nb(N_glor~elev_39+ppt+tmin+(1|YrSt)+(1|Stand)+Study, data = Ncovs) 
fm.a.tpi.tmin.asp <- lme4::glmer.nb(N_glor~tpi_11+tmin+asp_13+(1|YrSt)+(1|Stand)+Study, data = Ncovs) 
fm.a.slop.asp.tmin <- lme4::glmer.nb(N_glor~slop_13+asp_13+tmin+(1|YrSt)+(1|Stand)+Study, data = Ncovs)
fm.a.slop.asp.ppt <- lme4::glmer.nb(N_glor~slop_13+asp_13+ppt+(1|YrSt)+(1|Stand)+Study, data = Ncovs) 

# GEDI only
fm.g.fhd.pa10.cvr <- lme4::glmer.nb(N_glor~fhd_17+pa10_13+cvr_13+(1|YrSt)+(1|Stand)+Study, data = Ncovs) 
fm.g.pa20.rh98 <- lme4::glmer.nb(N_glor~pa20_27+rh98_15+(1|YrSt)+(1|Stand)+Study, data = Ncovs) 
fm.g.fhd.pa20.cvr <- lme4::glmer.nb(N_glor~fhd_17+pa20_27+cvr_13+(1|YrSt)+(1|Stand)+Study, data = Ncovs) 
fm.g.fhd <- lme4::glmer.nb(N_glor~fhd_17+(1|YrSt)+(1|Stand)+Study, data = Ncovs) 
fm.g.cvr.pa40.pa10 <- lme4::glmer.nb(N_glor~cvr_13+pa40_39+pa10_13+(1|YrSt)+(1|Stand)+Study, data = Ncovs)

# Hypothesis driven models, 
fm.h.elev.tmin.pa10 <- lme4::glmer.nb(N_glor~elev_39+tmin+pa10_13+(1|YrSt)+(1|Stand)+Study, data = Ncovs)
fm.h.tpi.fhd.pa10 <- lme4::glmer.nb(N_glor~tpi_11+fhd_17+pa10_13+(1|YrSt)+(1|Stand)+Study, data = Ncovs) 
fm.h.asp.rh98.cvr <- lme4::glmer.nb(N_glor~asp_13+rh98_15+cvr_13+(1|YrSt)+(1|Stand)+Study, data = Ncovs) 
fm.h.tmax.pa20.fhd <- lme4::glmer.nb(N_glor~tmax+pa20_27+fhd_17+(1|YrSt)+(1|Stand)+Study, data = Ncovs)
fm.h.tmax.ppt.pa40<- lme4::glmer.nb(N_glor~tmax+ppt+pa40_39+(1|YrSt)+(1|Stand)+Study, data = Ncovs) 

###################################################################################### 

glor.fm <- MuMIn::model.sel(fm.null,
                            
                            fm.a.elev.tmax.asp,fm.a.elev.tmin.ppt,fm.a.tpi.tmin.asp,fm.a.slop.asp.tmin,fm.a.slop.asp.ppt,
                            
                            fm.g.fhd.pa10.cvr,fm.g.pa20.rh98,fm.g.fhd.pa20.cvr,fm.g.fhd,fm.g.cvr.pa40.pa10,
                            
                            fm.h.elev.tmin.pa10,fm.h.tpi.fhd.pa10,fm.h.asp.rh98.cvr,fm.h.tmax.pa20.fhd,fm.h.tmax.ppt.pa40)

# GLOR model selection table
glor.fm

################
# NETO Models  #
################

# check for neto correlations 0.7 is the cutoff
round(cor(Ncovs[,c("asp_29", "elev_11", "tpi_11", "slop_39", "cvr_11","fhd_13","pa20_39","pa40_11","pa10_39", "rh50_13","rh75_13","rh98_15")]),2) 

# Null
fm.n.null<- lme4::glmer.nb(N_neto~(1|YrSt)+(1|Stand)+Study, data = Ncovs) 

# Abiotic
fm.a.n.elev.slop.asp <- lme4::glmer.nb(N_neto~elev_11+slop_39+asp_29+(1|YrSt)+(1|Stand)+Study, data = Ncovs)
fm.a.n.elev.asp.ppt <- lme4::glmer.nb(N_neto~elev_11+asp_29+ppt+(1|YrSt)+(1|Stand)+Study+Study, data = Ncovs)
fm.a.n.elev.tmin <- lme4::glmer.nb(N_neto~elev_11+tmin+(1|YrSt)+(1|Stand)+Study, data = Ncovs)
fm.a.n.tpi.slop.asp <- lme4::glmer.nb(N_neto~tpi_11+slop_39+asp_29+(1|YrSt)+(1|Stand)+Study, data = Ncovs)
fm.a.n.tpi.tmin.asp <- lme4::glmer.nb(N_neto~tpi_11+tmin+asp_29+(1|Stand)+(1|Year)+Study, data = Ncovs)

# GEDI only
fm.g.n.pa20.pa40 <- lme4::glmer.nb(N_neto~pa20_39+pa40_11+(1|YrSt)+(1|Stand)+Study, data = Ncovs) 
fm.g.n.fhd.cvr.pa20 <- lme4::glmer.nb(N_neto~fhd_13+pa20_39+cvr_11+(1|YrSt)+(1|Stand)+Study, data = Ncovs)
fm.g.n.rh98.pa20.pa10<- lme4::glmer.nb(N_neto~rh98_15+pa20_39+pa10_39+(1|YrSt)+(1|Stand)+Study, data = Ncovs) 
fm.g.n.pa10.cvr <- lme4::glmer.nb(N_neto~pa10_39+cvr_11+(1|YrSt)+(1|Stand)+Study, data = Ncovs)
fm.g.n.cvr.rh98.pa10 <- lme4::glmer.nb(N_neto~cvr_11+rh98_15+pa10_39+(1|YrSt)+(1|Stand)+Study, data = Ncovs)

# Hypothesis driven models
fm.h.n.elev.cvr.pa10 <- lme4::glmer.nb(N_neto~elev_11+cvr_11+pa10_39+(1|YrSt)+(1|Stand)+Study, data = Ncovs)
fm.h.n.ppt.pa20.pa40 <- lme4::glmer.nb(N_neto~ppt+pa20_39+pa40_11+(1|YrSt)+(1|Stand)+Study, data = Ncovs)
fm.h.n.elev.tmin.pa40 <- lme4::glmer.nb(N_neto~elev_11+tmin+pa40_11+(1|YrSt)+(1|Stand)+Study, data = Ncovs)
fm.h.n.asp.rh98.cvr <- lme4::glmer.nb(N_neto~asp_29+rh98_15+cvr_11+(1|YrSt)+(1|Stand)+Study, data = Ncovs)
fm.h.n.elev.pa20.fhd <- lme4::glmer.nb(N_neto~elev_11+pa20_39+fhd_13+(1|YrSt)+(1|Stand)+Study, data = Ncovs)
#######################################################################
neto.fm <- MuMIn::model.sel(fm.n.null,
                            fm.a.n.elev.slop.asp,fm.a.n.elev.asp.ppt,fm.a.n.elev.tmin,fm.a.n.tpi.slop.asp,fm.a.n.tpi.tmin.asp,
                            fm.g.n.pa20.pa40, fm.g.n.fhd.cvr.pa20,fm.g.n.rh98.pa20.pa10, fm.g.n.pa10.cvr, fm.g.n.cvr.rh98.pa10,
                            fm.h.n.elev.cvr.pa10,fm.h.n.ppt.pa20.pa40, fm.h.n.elev.tmin.pa40,fm.h.n.asp.rh98.cvr,fm.h.n.elev.pa20.fhd)

# NETO model selection table
neto.fm


# Build a function for overdispersion testing for both glor.fm and neto.fm
mods <- row.names(glor.fm) # swap for neto.fm  
d.test <- NULL
for(i in 1:length(mods)){
  fm=get(mods[i])
  prs <- overdisp_fun(fm)
  d.test <- c(d.test,prs)
}
round(d.test, 2)

##############################
# Uncertainty with Bootstrap #
# and upper/lower estimates  #
##############################

# Best model for GLOR or NETO. The process is the same just swap out the appropriate model
fm.mod <- fm.h.tpi.fhd.pa10
summary(fm.mod)
performance::performance(fm.mod) # fit values

# Run models at the upper and lower N estimates to determine how the model output changes
Ncovs$lcl_glor <- round(Ncovs$lcl_glor)
Ncovs$ucl_glor <- round(Ncovs$ucl_glor)

fm.mod.ucl <- lme4::glmer.nb(ucl_glor~tpi_11+fhd_17+pa10_13+(1|YrSt)+(1|Stand)+Study, data = Ncovs) # upper CI rerun
fm.mod.lcl <- lme4::glmer.nb(lcl_glor~tpi_11+fhd_17+pa10_13+(1|YrSt)+(1|Stand)+Study, data = Ncovs) # lower CI rerun


# fixed effects bootstrap
feEx <- rbind(merTools::FEsim(fm.mod, 5000),
              merTools::FEsim(fm.mod.ucl, 5000),
              merTools::FEsim(fm.mod.lcl, 5000))

# random effects
reEx <- rbind(merTools::REsim(fm.mod),
              merTools::REsim(fm.mod.ucl),
              merTools::REsim(fm.mod.lcl))
# then plot as needed

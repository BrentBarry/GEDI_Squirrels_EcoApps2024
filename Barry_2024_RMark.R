#############################
# Abundance Estimates       #
# in rmark using Huggins RD #
# models                    #
#                           #
# BRB 1/10/2024             #
#############################

# load packages
library(RMark)
library(ggplot2)
library(dplyr)
library(knitr)


rm(list=ls())

# read in the data either glor or neto
dater <- read.csv("./doc/final_files/data_glor.csv")
dater<-subset(dater,ch!="000000000000000000000000000000000000000000000000000000000000000000000000") # remove blanks
intervals<-c(rep(c(0,0,0,0,0,0,0,0,0,0,0,1),5),c(0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0)) # 1 indicates that one interval of time (e.g., year) has passed between primary occasions, 0's are n-1 secondary occasions indicating closure periods.  

data.processed=process.data(dater,model="RDHuggins",time.intervals=intervals,groups=c("area","grid"), begin.time=2014)
data.processed$group.covariates$area<-relevel(data.processed$group.covariates$area,ref='HJA') 


# Parameter Fix Function --------------------------------------------------
# The example shown is for GLOR but neto has slightly different trapping periods explained in the methods of the main text
fix.params=function(data.ddl)
{
  # Fix survival estimates for years where no trapping occurred including when the HJA grids drop down in size
  data.ddl$S$fix=NA
  data.ddl$S$fix[data.ddl$S$area=="UMP" & data.ddl$S$time%in%2016:2021]=0 # survival fix the years appropriately
  data.ddl$S$fix[data.ddl$S$area=="SIU" & data.ddl$S$time%in% 2016:2021]=0 # Siuslaw only trapped 2014-2016
  
  # These grids were not trapped in later years
  data.ddl$S$fix[data.ddl$S$area=="HJA" & data.ddl$S$grid==1 & data.ddl$S$time%in%2018:2021]=0  
  data.ddl$S$fix[data.ddl$S$area=="HJA" & data.ddl$S$grid==3 & data.ddl$S$time%in%2018:2021]=0
  data.ddl$S$fix[data.ddl$S$area=="HJA" & data.ddl$S$grid==4 & data.ddl$S$time%in%2018:2021]=0
  data.ddl$S$fix[data.ddl$S$area=="HJA" & data.ddl$S$grid==6 & data.ddl$S$time%in%2018:2021]=0
  data.ddl$S$fix[data.ddl$S$area=="HJA" & data.ddl$S$grid==7 & data.ddl$S$time%in%2018:2021]=0
  data.ddl$S$fix[data.ddl$S$area=="HJA" & data.ddl$S$grid==9 & data.ddl$S$time%in%2018:2021]=0
  
  # Gamma fixes according to trapped years
  data.ddl$GammaDoublePrime$fix=NA
  data.ddl$GammaDoublePrime$fix[data.ddl$GammaDoublePrime$area=="UMP" & data.ddl$GammaDoublePrime$time%in%2016:2021]=0
  data.ddl$GammaDoublePrime$fix[data.ddl$GammaDoublePrime$area=="SIU" & data.ddl$GammaDoublePrime$time%in%2016:2021]=0
  data.ddl$GammaDoublePrime$fix[data.ddl$GammaDoublePrime$area=="HJA" & data.ddl$GammaDoublePrime$grid==1 & data.ddl$GammaDoublePrime$time%in%2018:2021]=0 
  data.ddl$GammaDoublePrime$fix[data.ddl$GammaDoublePrime$area=="HJA" & data.ddl$GammaDoublePrime$grid==3 & data.ddl$GammaDoublePrime$time%in%2018:202001]=0
  data.ddl$GammaDoublePrime$fix[data.ddl$GammaDoublePrime$area=="HJA" & data.ddl$GammaDoublePrime$grid==4 & data.ddl$GammaDoublePrime$time%in%2018:2021]=0 
  data.ddl$GammaDoublePrime$fix[data.ddl$GammaDoublePrime$area=="HJA" & data.ddl$GammaDoublePrime$grid==6 & data.ddl$GammaDoublePrime$time%in%2018:2021]=0 
  data.ddl$GammaDoublePrime$fix[data.ddl$GammaDoublePrime$area=="HJA" & data.ddl$GammaDoublePrime$grid==7& data.ddl$GammaDoublePrime$time%in%2018:2021]=0 
  data.ddl$GammaDoublePrime$fix[data.ddl$GammaDoublePrime$area=="HJA" & data.ddl$GammaDoublePrime$grid==9 & data.ddl$GammaDoublePrime$time%in%2018:2021]=0 
  
  data.ddl$GammaPrime$fix=NA
  data.ddl$GammaPrime$fix[data.ddl$GammaPrime$area=="UMP" & data.ddl$GammaPrime$time%in%2016:2021]=0
  data.ddl$GammaPrime$fix[data.ddl$GammaPrime$area=="SIU" & data.ddl$GammaPrime$time%in%2016:2021]=0
  data.ddl$GammaPrime$fix[data.ddl$GammaPrime$area=="HJA" & data.ddl$GammaPrime$grid==1 & data.ddl$GammaPrime$time%in%2018:2021]=0
  data.ddl$GammaPrime$fix[data.ddl$GammaPrime$area=="HJA" & data.ddl$GammaPrime$grid==3 & data.ddl$GammaPrime$time%in%2018:2021]=0
  data.ddl$GammaPrime$fix[data.ddl$GammaPrime$area=="HJA" & data.ddl$GammaPrime$grid==4 & data.ddl$GammaPrime$time%in%2018:2021]=0
  data.ddl$GammaPrime$fix[data.ddl$GammaPrime$area=="HJA" & data.ddl$GammaPrime$grid==6 & data.ddl$GammaPrime$time%in%2018:2021]=0
  data.ddl$GammaPrime$fix[data.ddl$GammaPrime$area=="HJA" & data.ddl$GammaPrime$grid==7 & data.ddl$GammaPrime$time%in%2018:2021]=0
  data.ddl$GammaPrime$fix[data.ddl$GammaPrime$area=="HJA" & data.ddl$GammaPrime$grid==9 & data.ddl$GammaPrime$time%in%2018:2021]=0
  
  # recapture probability
  data.ddl$p$fix=NA
  data.ddl$p$fix[data.ddl$p$area=="UMP" & data.ddl$p$session%in%2017:2021]=0  
  data.ddl$p$fix[data.ddl$p$area=="SIU" & data.ddl$p$session%in%2017:2021]=0
  data.ddl$p$fix[data.ddl$p$area=="HJA" & data.ddl$p$session%in%2017:2021 & data.ddl$p$time%in%9:12]=0 # switch to 8 day capture period
  data.ddl$p$fix[data.ddl$p$area=="HJA" & data.ddl$p$grid == 1 & data.ddl$p$session%in%2019:2021]=0 
  data.ddl$p$fix[data.ddl$p$area=="HJA" & data.ddl$p$grid == 3 & data.ddl$p$session%in%2019:2021]=0 
  data.ddl$p$fix[data.ddl$p$area=="HJA" & data.ddl$p$grid == 4 & data.ddl$p$session%in%2019:2021]=0 
  data.ddl$p$fix[data.ddl$p$area=="HJA" & data.ddl$p$grid == 6 & data.ddl$p$session%in%2019:2021]=0 
  data.ddl$p$fix[data.ddl$p$area=="HJA" & data.ddl$p$grid == 7 & data.ddl$p$session%in%2019:2021]=0 
  data.ddl$p$fix[data.ddl$p$area=="HJA" & data.ddl$p$grid == 9 & data.ddl$p$session%in%2019:2021]=0 

  # capture probability
  data.ddl$c$fix=NA
  data.ddl$c$fix[data.ddl$c$area=="UMP" & data.ddl$c$session%in%2017:2021]=0  
  data.ddl$c$fix[data.ddl$c$area=="SIU" & data.ddl$c$session%in%2017:2021]=0
  data.ddl$c$fix[data.ddl$c$area=="HJA" & data.ddl$c$session%in%2017:2021 & data.ddl$c$time%in%9:12]=0 # switch to 8 day capture period
  data.ddl$c$fix[data.ddl$c$area=="HJA" & data.ddl$c$grid == 1 & data.ddl$c$session%in%2019:2021]=0 
  data.ddl$c$fix[data.ddl$c$area=="HJA" & data.ddl$c$grid == 3 & data.ddl$c$session%in%2019:2021]=0 
  data.ddl$c$fix[data.ddl$c$area=="HJA" & data.ddl$c$grid == 4 & data.ddl$c$session%in%2019:2021]=0 
  data.ddl$c$fix[data.ddl$c$area=="HJA" & data.ddl$c$grid == 6 & data.ddl$c$session%in%2019:2021]=0 
  data.ddl$c$fix[data.ddl$c$area=="HJA" & data.ddl$c$grid == 7 & data.ddl$c$session%in%2019:2021]=0 
  data.ddl$c$fix[data.ddl$c$area=="HJA" & data.ddl$c$grid == 9 & data.ddl$c$session%in%2019:2021]=0 
  
  return(data.ddl)
}


##################
# GAMMA Models   #
##################
Gamma.models=function(data.processed)
{
  
  data.ddl=make.design.data(data.processed)
  data.ddl<-fix.params(data.ddl)
  #
  # Define parameter models
  #
  c.global=list(formula=~area*time)
  p.global=list(formula=~area*time)
  S.global=list(formula=~area*time)
  
  Gamma.Prime=list(formula=~1)
  Gamma.Double.Prime=list(formula=~1)
  Gammas.Zero=list(formula=~1,fixed=0)
  Gammas.Double.Zero=list(formula=~1,fixed=0)
  Gammas.random=list(formula=~1, share=TRUE)
  
  Gamma.Prime.GRID=list(formula=~grid)
  Gamma.Double.Prime.GRID=list(formula=~grid)
  Gammas.shared.GRID=list(formula=~grid, share=TRUE)
  
  Gamma.Prime.Area=list(formula=~area)
  Gamma.Double.Prime.Area=list(formula=~area)
  Gammas.shared.Area=list(formula=~area, share=TRUE)
  
  
  # Candidate Models
  mod_separate.gammas<-mark(data.processed,data.ddl,model.parameters=list(c=c.global, p=p.global,S=S.global,GammaDoublePrime=Gamma.Double.Prime,GammaPrime=Gamma.Prime))
  mod_gammas.random<-mark(data.processed,data.ddl,model.parameters=list(c=c.global, p=p.global,S=S.global,GammaDoublePrime=Gammas.random))
  mod_gammas.zero<-mark(data.processed,data.ddl,model.parameters=list(c=c.global, p=p.global,S=S.global,GammaDoublePrime=Gammas.Double.Zero,GammaPrime=Gammas.Zero))
  
  mod_Gamma.Double.GRID<-mark(data.processed,data.ddl,model.parameters=list(c=c.global, p=p.global,S=S.global,GammaDoublePrime=Gamma.Double.Prime.GRID,GammaPrime=Gamma.Prime))
  mod_Gamma.GRID<-mark(data.processed,data.ddl,model.parameters=list(c=c.global, p=p.global,S=S.global,GammaDoublePrime=Gamma.Double.Prime,GammaPrime=Gamma.Prime.GRID))
  mod_Gamma.Shared.GRID<-mark(data.processed,data.ddl,model.parameters=list(c=c.global, p=p.global,S=S.global,GammaDoublePrime=Gammas.shared.GRID))
  mod_Gammas.GRID<-mark(data.processed,data.ddl,model.parameters=list(c=c.global, p=p.global,S=S.global,GammaDoublePrime=Gamma.Double.Prime.GRID,GammaPrime=Gamma.Prime.GRID))
  
  mod_Gamma.Double.Area<-mark(data.processed,data.ddl,model.parameters=list(c=c.global, p=p.global,S=S.global,GammaDoublePrime=Gamma.Double.Prime.Area,GammaPrime=Gamma.Prime))
  mod_Gamma.Area<-mark(data.processed,data.ddl,model.parameters=list(c=c.global, p=p.global,S=S.global,GammaDoublePrime=Gamma.Double.Prime,GammaPrime=Gamma.Prime.Area))
  mod_Gamma.Shared.Area<-mark(data.processed,data.ddl,model.parameters=list(c=c.global, p=p.global,S=S.global,GammaDoublePrime=Gammas.shared.Area))
  mod_Gammas.Area<-mark(data.processed,data.ddl,model.parameters=list(c=c.global, p=p.global,S=S.global,GammaDoublePrime=Gamma.Double.Prime.Area,GammaPrime=Gamma.Prime.Area))
  
  return(collect.models())
  
}
gamma.results=Gamma.models(data.processed)

(gamma.table<-kable(model.table(model.list = gamma.results, sort = TRUE, adjust = TRUE,
                                pf = 1, use.lnl = FALSE, use.AIC = FALSE,
                                model.name = FALSE)))

# Fit behavior models
behavior.models=function(data.processed)
{
  data.ddl=make.design.data(data.processed)
  data.ddl<-fix.params(data.ddl)
  #
  # Define parameter models
  #
  
  S.global=list(formula=~area*time)
  gammas.best=list(formula=~1,fixed=0)
  
  c.dot=list(formula=~1)
  p.dot=list(formula=~1)
  c.Time=list(formula=~time)
  p.Time=list(formula=~time)
  p.c.share=list(formula=~1, share=TRUE)
  p.c.Time.share=list(formula=~time, share=TRUE)
  
  # Candidate Models
  mod_behavior<-mark(data.processed,data.ddl,model.parameters=list(c=c.dot, p=p.dot,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_p.c.share<-mark(data.processed,data.ddl,model.parameters=list(p=p.c.share,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_p.c.Time.share<-mark(data.processed,data.ddl,model.parameters=list(p=p.c.Time.share,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_p.Time<-mark(data.processed,data.ddl,model.parameters=list(c=c.dot, p=p.Time,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_c.Time<-mark(data.processed,data.ddl,model.parameters=list(c=c.Time, p=p.dot,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_p.c.Time<-mark(data.processed,data.ddl,model.parameters=list(c=c.Time, p=p.Time,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  
  return(collect.models())
  
}
behavior.results=behavior.models(data.processed)

(behavior.table<-kable(model.table(model.list = behavior.results, sort = TRUE, adjust = TRUE,
                                   pf = 1, use.lnl = FALSE, use.AIC = FALSE,
                                   model.name = FALSE)))

# Capture with behavioral models
 c.model=function(data.processed)
{
  data.ddl=make.design.data(data.processed)
  data.ddl<-fix.params(data.ddl)
  #
  # Define constant parameter models
  c.global=list(formula=~area*time)
  p.global=list(formula=~area*time)
  S.global=list(formula=~area*time)
  
  # Best models 
  gammas.best=list(formula=~1,fixed=0)
  
  # univariate
  c.dot=list(formula=~1)
  c.grid=list(formula=~grid)
  c.time=list(formula=~time)
  c.Time=list(formula=~Time)
  c.area=list(formula=~area)
  c.session=list(formula=~session)

  #multivariate
  c.session.grid =list(formula=~session+grid)
  c.session.area =list(formula=~session+area)
  c.session.time =list(formula=~session+time)
  c.session.Time =list(formula=~session+Time)
  c.grid.Time =list(formula=~grid+Time)
  c.grid.time =list(formula=~grid+time)
  c.grid.area =list(formula=~grid+area)
  c.area.time =list(formula=~area+time)
  c.area.Time =list(formula=~area+Time)
  
  # Run models
  mod_null <- mark(data.processed,model.parameters=list(c=list(formula=~1), p=list(formula=~1),S=list(formula=~1),GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  
  mod_c.dot <- mark(data.processed,model.parameters=list(p=p.global, c=c.dot,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_c.grid <- mark(data.processed,model.parameters=list(p=p.global, c=c.grid,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_c.time <- mark(data.processed,model.parameters=list(p=p.global, c=c.time,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_c.Time <- mark(data.processed,model.parameters=list(p=p.global, c=c.Time,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_c.area <- mark(data.processed,model.parameters=list(p=p.global, c=c.area,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_c.session <- mark(data.processed,model.parameters=list(p=p.global, c=c.session,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  
  # multivariate
  mod_c.session.grid <- mark(data.processed,model.parameters=list(p=p.global,c=c.session.grid,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_c.session.area <- mark(data.processed,model.parameters=list(p=p.global,c=c.session.area,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_c.session.time <- mark(data.processed,model.parameters=list(p=p.global,c=c.session.time,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_c.session.Time <- mark(data.processed,model.parameters=list(p=p.global,c=c.session.Time,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_c.grid.Time <- mark(data.processed,model.parameters=list(p=p.global,c=c.grid.Time,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_c.grid.time <- mark(data.processed,model.parameters=list(p=p.global,c=c.grid.time,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_c.grid.area <- mark(data.processed,model.parameters=list(p=p.global,c=c.grid.area,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_c.area.time <- mark(data.processed,model.parameters=list(p=p.global,c=c.area.time,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_c.area.Time <- mark(data.processed,model.parameters=list(p=p.global,c=c.area.Time,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  
  return(collect.models())
  
}
c.results=c.model(data.processed)

(c.table <- kable(model.table(model.list = c.results, sort = TRUE,adjust = TRUE,
                                        pf = 1, use.lnl = FALSE, use.AIC = FALSE,
                                        model.name = FALSE), digits=2, 
                            col.names=c("S", "GammaPR","GammaDPR","p","c","model","npar","AICc", "DeltaAICc","wi","Deviance"),
                            caption="Table S1:Probability of capture models using null values for capture and survival, while holding gamma and double gamma fixed at zero"))



# Capture with behavioral models
p.model=function(data.processed)
{
  data.ddl=make.design.data(data.processed)
  data.ddl<-fix.params(data.ddl)
  #
  # Define constant parameter models
  p.global=list(formula=~area*time)
  S.global=list(formula=~area*time)
  
  # Define best model 
  gammas.best=list(formula=~1,fixed=0)
  c.best=list(formula=~session+time)

    # univariate
  p.dot=list(formula=~1)
  p.grid=list(formula=~grid)
  p.time=list(formula=~time)
  p.Time=list(formula=~Time)
  p.area=list(formula=~area)
  p.session=list(formula=~session)
  
  #multivariate
  p.session.grid =list(formula=~session+grid)
  p.session.area =list(formula=~session+area)
  p.session.time =list(formula=~session+time)
  p.session.Time =list(formula=~session+Time)
  p.grid.Time =list(formula=~grid+Time)
  p.grid.time =list(formula=~grid+time)
  p.grid.area =list(formula=~grid+area)
  p.area.time =list(formula=~area+time)
  p.area.Time =list(formula=~area+Time)
  
  # Run models
  mod_null <- mark(data.processed,model.parameters=list(c=list(formula=~1), p=list(formula=~1),S=list(formula=~1),GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_p.dot <- mark(data.processed,model.parameters=list(c=c.best, p=p.dot,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_p.grid <- mark(data.processed,model.parameters=list(c=c.best, p=p.grid,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_p.time <- mark(data.processed,model.parameters=list(c=c.best, p=p.time,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_p.Time <- mark(data.processed,model.parameters=list(c=c.best, p=p.Time,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_p.area <- mark(data.processed,model.parameters=list(c=c.best, p=p.area,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_p.session <- mark(data.processed,model.parameters=list(c=c.best, p=p.session,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  
  #multivariate
  mod_p.session.grid <- mark(data.processed,model.parameters=list(c=c.best, p=p.session.grid,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_p.session.area <- mark(data.processed,model.parameters=list(c=c.best, p=p.session.area,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_p.session.time <- mark(data.processed,model.parameters=list(c=c.best, p=p.session.time,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_p.session.Time <- mark(data.processed,model.parameters=list(c=c.best, p=p.session.Time,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_p.grid.Time <- mark(data.processed,model.parameters=list(c=c.best, p=p.grid.Time,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_p.grid.time <- mark(data.processed,model.parameters=list(c=c.best, p=p.grid.time,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_p.grid.area <- mark(data.processed,model.parameters=list(c=c.best, p=p.grid.area,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_p.area.time <- mark(data.processed,model.parameters=list(c=c.best, p=p.area.time,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  mod_p.area.Time <- mark(data.processed,model.parameters=list(c=c.best, p=p.area.Time,S=S.global,GammaDoublePrime=gammas.best,GammaPrime=gammas.best))
  return(collect.models())
  
}
p.results=p.model(data.processed)

(p.table <- kable(model.table(model.list = p.results, sort = TRUE,adjust = TRUE,
                                                 pf = 1, use.lnl = FALSE, use.AIC = FALSE,
                                                 model.name = FALSE), digits=2, 
                                     col.names=c("S", "GammaPR","GammaDPR","p","c","model","npar","AICc", "DeltaAICc","wi","Deviance"),
                                     caption="Table S1:Probability of recapture models using null values for capture and survival, while holding gamma and double gamma fixed at zero"))


s.model=function(data.processed)
{
  data.ddl=make.design.data(data.processed)
  data.ddl<-fix.params(data.ddl)
  #
  # Define constant parameter models
  c.null=list(formula=~1)
  p.null=list(formula=~1)
  S.null=list(formula=~1)
  gamma.null=list(formula=~1)
  
  # Define best models
  gammas.best =list(formula = ~ 1, fixed=0)
  c.best=list(formula=~session+time)
  p.best=list(formula=~area+time)
  
  # Nuisance
  S.grid =list(formula=~grid)
  S.area =list(formula=~area)
  S.time =list(formula=~time)
  S.Time =list(formula=~Time)
  S.grid.time =list(formula=~grid+time)
  S.grid.area =list(formula=~grid+area)
  S.area.time =list(formula=~area+time)
  
  # S GEDI variants    
  
  S.cvr =list(formula = ~ cvr_11)
  S.fhd =list(formula = ~ fhd_11)
  S.pa10 =list(formula = ~ pa10_11)
  S.pa20 =list(formula = ~ pa20_11)
  S.pa40 =list(formula = ~ pa40_11)
  S.rh50 =list(formula = ~ rh50_11)
  S.rh98 =list(formula = ~ rh98_11)
  
  S.cvr.rh98 =list(formula = ~ cvr_11+rh98_11)
  S.cvr.pa40 =list(formula = ~ cvr_11+pa40_11)
  S.cvr.pa10 =list(formula = ~ cvr_11+pa10_11)
  S.cvr.fhd =list(formula = ~ cvr_11+fhd_11)
  S.fhd.pa10 =list(formula = ~ fhd_11+pa10_11)
  S.fhd.pa20 =list(formula = ~ fhd_11+pa20_11)
  S.pa10.pa40 =list(formula = ~ pa10_11+pa40_11)
  S.pa10.pa20 =list(formula = ~ pa10_11+pa20_11)
  S.pa10.rh50 =list(formula = ~ pa10_11+rh50_11)
  S.pa10.rh98 =list(formula = ~ pa10_11+rh98_11)

  # S GEDI + 1 cov
  S.area.cvr =list(formula = ~ area+cvr_11)
  S.area.fhd =list(formula = ~ area+fhd_11)
  S.area.pa10 =list(formula = ~ area+pa10_11)
  S.area.pa20 =list(formula = ~ area+pa20_11)
  S.area.pa40 =list(formula = ~ area+pa40_11)
  S.area.rh50 =list(formula = ~ area+rh50_11)
  S.area.rh98 =list(formula = ~ area+rh98_11)
  
  S.grid.cvr =list(formula = ~ grid+cvr_11)
  S.grid.fhd =list(formula = ~ grid+fhd_11)
  S.grid.pa10 =list(formula = ~ grid+pa10_11)
  S.grid.pa20 =list(formula = ~ grid+pa20_11)
  S.grid.pa40 =list(formula = ~ grid+pa40_11)
  S.grid.rh50 =list(formula = ~ grid+rh50_11)
  S.grid.rh98 =list(formula = ~ grid+rh98_11)
  
  S.time.cvr =list(formula = ~ time+cvr_11)
  S.time.fhd =list(formula = ~ time+fhd_11)
  S.time.pa10 =list(formula = ~ time+pa10_11)
  S.time.rh50 = list(formula = ~ time+rh50_11)
  S.time.pa20 =list(formula = ~ time+pa20_11)
  S.time.pa40 =list(formula = ~ time+pa40_11)
  S.time.rh98 =list(formula = ~ time+rh98_11)
  
  # Base models
  mod_s.null <-mark(data.processed,data.ddl,model.parameters=list(p= p.null,c=c.null,S=S.null,GammaDoublePrime=gamma.null, GammaPrime=gamma.null))
  mod_s.dot <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,S=S.null,GammaDoublePrime=gammas.best, GammaPrime= gammas.best))
  mod_s.grid <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.grid, GammaPrime= gammas.best))
  mod_s.area <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.area, GammaPrime= gammas.best))
  mod_s.time <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.time, GammaPrime= gammas.best))
  mod_s.Time <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.Time, GammaPrime= gammas.best))
  mod_s.area.time <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.area.time, GammaPrime= gammas.best))
  mod_s.grid.time <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.grid.time, GammaPrime= gammas.best))
  mod_s.grid.area <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.grid.area, GammaPrime= gammas.best))
  
  
  # GEDI univariate
  mod_s.cvr <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,S=S.cvr,GammaDoublePrime=gammas.best, GammaPrime= gammas.best))
  mod_s.fhd <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,S=S.fhd,GammaDoublePrime=gammas.best, GammaPrime= gammas.best))
  mod_s.pa20 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,S=S.pa20,GammaDoublePrime=gammas.best, GammaPrime= gammas.best))
  mod_s.pa40 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,S=S.pa40,GammaDoublePrime=gammas.best, GammaPrime= gammas.best))
  mod_s.pa10 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,S=S.pa10,GammaDoublePrime=gammas.best, GammaPrime= gammas.best))
  mod_s.rh50 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,S=S.rh50,GammaDoublePrime=gammas.best, GammaPrime= gammas.best))
  mod_s.rh98 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,S=S.rh98,GammaDoublePrime=gammas.best, GammaPrime= gammas.best))

  # GEDI multivariate
  mod_s.cvr.rh98 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,S=S.cvr.rh98,GammaDoublePrime=gammas.best, GammaPrime= gammas.best))
  mod_s.cvr.pa40 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,S=S.cvr.pa40,GammaDoublePrime=gammas.best, GammaPrime= gammas.best))
  mod_s.cvr.pa10 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,S=S.cvr.pa10,GammaDoublePrime=gammas.best, GammaPrime= gammas.best))
  mod_s.cvr.fhd <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,S=S.cvr.fhd,GammaDoublePrime=gammas.best, GammaPrime= gammas.best))
  mod_s.fhd.pa10 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,S=S.fhd.pa10,GammaDoublePrime=gammas.best, GammaPrime= gammas.best))
  mod_s.fhd.pa20 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,S=S.fhd.pa20,GammaDoublePrime=gammas.best, GammaPrime= gammas.best))
  mod_s.pa10.pa40 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,S=S.pa10.pa40,GammaDoublePrime=gammas.best, GammaPrime= gammas.best))
  mod_s.pa10.pa20 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,S=S.pa10.pa20,GammaDoublePrime=gammas.best, GammaPrime= gammas.best))
  mod_s.pa10.rh50 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,S=S.pa10.rh50,GammaDoublePrime=gammas.best, GammaPrime= gammas.best))
  mod_s.pa10.rh98 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,S=S.pa10.rh98,GammaDoublePrime=gammas.best, GammaPrime= gammas.best))
  
  
  # GEDI and 1 cov
  mod_s.area.cvr <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.area.cvr, GammaPrime= gammas.best))
  mod_s.area.fhd <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.area.fhd, GammaPrime= gammas.best))
  mod_s.area.pa20 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.area.pa20, GammaPrime= gammas.best))
  mod_s.area.pa40 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.area.pa40, GammaPrime= gammas.best))
  mod_s.area.pa10 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.area.pa10, GammaPrime= gammas.best))
  mod_s.area.rh50 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.area.rh50, GammaPrime= gammas.best))
  mod_s.area.rh98 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.area.rh98, GammaPrime= gammas.best))

#  time + cov
  mod_s.time.cvr <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.time.cvr, GammaPrime= gammas.best))
  mod_s.time.fhd <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.time.fhd, GammaPrime= gammas.best))
  mod_s.time.pa20 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.time.pa20, GammaPrime= gammas.best))
  mod_s.time.pa40 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.time.pa40, GammaPrime= gammas.best))
  mod_s.time.pa10 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.time.pa10, GammaPrime= gammas.best))
  mod_s.time.rh50 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.time.rh50, GammaPrime= gammas.best))
  mod_s.time.rh98 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.time.rh98, GammaPrime= gammas.best))
  
  # grid + cov
  mod_s.grid.cvr <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.grid.cvr, GammaPrime= gammas.best))
  mod_s.grid.fhd <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.grid.fhd, GammaPrime= gammas.best))
  mod_s.grid.pa20 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.grid.pa20, GammaPrime= gammas.best))
  mod_s.grid.pa40 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.grid.pa40, GammaPrime= gammas.best))
  mod_s.grid.pa10 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.grid.pa10, GammaPrime= gammas.best))
  mod_s.grid.rh50 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.grid.rh50, GammaPrime= gammas.best))
  mod_s.grid.rh98 <-mark(data.processed,data.ddl,model.parameters=list(p=p.best,c=c.best,GammaDoublePrime=gammas.best,S=S.grid.rh98, GammaPrime= gammas.best))

  
  return(collect.models())
  
}

s.results=s.model(data.processed)

(s.table<-kable(model.table(model.list = s.results, sort = TRUE, adjust = TRUE,
                                 pf = 1, use.lnl = FALSE, use.AIC = FALSE,
                                 model.name = FALSE), digits=2,
                     col.names=c("S", "GammaPR","GammaDPR","p","c","model","npar","AICc", "DeltaAICc","wi","Deviance"),
                     caption="Table S4: Apparent annual survival models using gedi covariates and the most supported parameters for gamma, c and p"))

rm(list=ls())
library(devtools)
library(deSolve)
library(rootSolve)
library(ReacTran)
library(rodeo)
library(readxl)
library(FME)
library(coda)
library(ggplot2)
library(ggpubr)



####################################################
# SETTING INITIAL GUESS VALUES FOR INLET AND OUTLET#
####################################################


st_sim<-unclass(as.POSIXlt(Sys.time()))
parm_2 <- as.data.frame(read_excel("parameters.xlsx", 
                                   sheet = "pars",
                                   col_names = TRUE,
                                   range="A1:F2"))

obs_data_out<- as.data.frame(read_excel(path= "data_exp.xlsx",
                                sheet="exp_3", 
                                col_names=TRUE))
obs_data_in<- as.data.frame(read_excel(path= "data_exp.xlsx",
                                        sheet="exp3_in", 
                                        col_names=TRUE))

###########################################################
# DETERMINING DECAY RATE USING INLET CONCENTRATION CURVES #
###########################################################

parmguess_in <- c(mu_d=0.02)
Model_fun<-function(p){
  derivs<-function(t,state,p){
    with(as.list(c(state,p)),{
      dC<- -mu_d*C
      return(list(c(dC)))
    })
  }
  state<- c(C=0.03)
  times <- seq(from=15,to=1500,by=5)
  t<-Data$time
  out<- ode(y=state,times=times,func = derivs,parms = p)
  return(out)
}
Model_Cost<-function(p){
  out_df<-Model_fun(p)
  return(modCost(model=out_df,obs=obs_data_in))
}
Fit1_in<-modFit(f = Model_Cost, p = parmguess_in ,method = 'bobyqa')




######################################################################
######  FUNCTIONS AND INITIALIZATION FOR OUTLET CONCENTRATION ########
######################################################################

Model_fun<-function(p){
  asMatrix<- function(df){
    matrix(unlist(df[,2:ncol(df)]), 
           nrow=nrow(df),
           ncol= ncol(df)-1,
           dimnames =list(df[,1], names(df)[2:ncol(df)]))
  }
  L <- 20
  nLayers <- 200
  dx <- L/nLayers
  
  ############### TIME DEFINITION ###########
  ###### INITIALIZATION OF ALL OBJECTS#######
  vars <- as.data.frame(read_excel("parameters.xlsx", sheet = "vars", col_names = TRUE,  skip = 0))
  pars <- as.data.frame(read_excel("parameters.xlsx", sheet = "pars", col_names = TRUE,  skip = 0))
  funs <- as.data.frame(read_excel("parameters.xlsx", sheet = "funcs",col_names = TRUE,  skip = 0))
  pros <- as.data.frame(read_excel("parameters.xlsx", sheet = "pros", col_names = TRUE,  skip = 0))
  stoi <- as.data.frame(read_excel("parameters.xlsx", sheet = "stoi", col_names = TRUE,  skip = 0))
  
  model <- rodeo$new(
    vars = vars[1:3],
    pars = pars[1:3],
    funs= funs,
    pros = pros,
    stoi = asMatrix(stoi),
    asMatrix=TRUE,
    dim= c(nLayers)
  )
  
  model$compile("boundary.f95")
  tmp_var<- vars$values
  names(tmp_var)<- vars$name
  tmp_var <- matrix (rep(tmp_var,
                         each=nLayers),
                     nrow=nLayers, 
                     ncol=length(tmp_var),
                     dimnames=list(NULL,
                                   names(tmp_var)))
  
  model$setVars(tmp_var)
  tmp<-c(as.numeric(p),pars$values_3[2:10])
  names(tmp)<-pars$name
  tmp <- matrix(rep(tmp,
                    each=nLayers), 
                nrow=nLayers, 
                ncol=length(tmp),
                dimnames=list(NULL,names(tmp)))
  tmp[1,"leftmost"]<-1
  model$setPars(tmp)
  # t<-Data$time
   t<-seq(from=1,to=250,by=0.1)
  out1<- model$dynamics(times=t)
  out_df<-(cbind(out1[,"time"],out1[,"C.165"]))
  colnames(out_df)=c("time","conc")
  return(out_df)
  }

 Fit0_1<-Model_fun(parm_guess)

Model_cost<-function(p){
  
  out_df<-Model_fun(p)
  Data<-as.data.frame((read_excel(path= "data_exp.xlsx",
                                  sheet="exp_3", 
                                  col_names=TRUE))
                       )
  cost_fn<-modCost(model=out_df, obs=Data) 
}

Fit <- modFit(f = Model_cost,
                 p = parm_guess,
                 method = 'bobyqa',
                 lower = rep(0,1),
                 upper = rep(1,1)
)
Fit1<-Model_fun(Fit$par)
summary(Fit)
##################################################
#       CALCUALTING SEDIMENT RETAINED DATA      #
##################################################


Model_fun_soil<-function(p){
  asMatrix<- function(df){
    matrix(unlist(df[,2:ncol(df)]), 
           nrow=nrow(df),
           ncol= ncol(df)-1,
           dimnames =list(df[,1], names(df)[2:ncol(df)]))
  }
  L <- 30
  nLayers <- 400
  dx <- L/nLayers
  
  ############### TIME DEFINITION ###########
  ###### INITIALIZATION OF ALL OBJECTS#######
  vars <- as.data.frame(read_excel("parameters.xlsx", sheet = "vars", col_names = TRUE,  skip = 0))
  pars <- as.data.frame(read_excel("parameters.xlsx", sheet = "pars", col_names = TRUE,  skip = 0))
  funs <- as.data.frame(read_excel("parameters.xlsx", sheet = "funcs",col_names = TRUE,  skip = 0))
  pros <- as.data.frame(read_excel("parameters.xlsx", sheet = "pros", col_names = TRUE,  skip = 0))
  stoi <- as.data.frame(read_excel("parameters.xlsx", sheet = "stoi", col_names = TRUE,  skip = 0))
  
  model <- rodeo$new(
    vars = vars[1:3],
    pars = pars[1:3],
    funs= funs,
    pros = pros,
    stoi = asMatrix(stoi),
    asMatrix=TRUE,
    dim= c(nLayers)
  )
  
  model$compile("boundary.f95")
  tmp_var<- vars$values
  names(tmp_var)<- vars$name
  tmp_var <- matrix (rep(tmp_var,
                         each=nLayers),
                     nrow=nLayers, 
                     ncol=length(tmp_var),
                     dimnames=list(NULL,
                                   names(tmp_var)))
  
  model$setVars(tmp_var)
  tmp<-c(as.numeric(p),pars$values_3[2:10])
  names(tmp)<-pars$name
  tmp <- matrix(rep(tmp,
                    each=nLayers), 
                nrow=nLayers, 
                ncol=length(tmp),
                dimnames=list(NULL,names(tmp)))
  tmp[1,"leftmost"]<-1
  model$setPars(tmp)
  t<-seq(from=0,to=250,by=0.1)
  out1<- model$dynamics(times=t)
 return(out1)
}

Fit_soil<-Model_fun_soil(Fit1_1$par)
Fit_soil_end<-as.data.frame(Fit_soil[2280,402:801])
C_1<-mean(Fit_soil_end$`Fit_soil[2280, 402:801]`[1:60])
C_2<-mean(Fit_soil_end$`Fit_soil[2280, 402:801]`[60:120])
C_3<-mean(Fit_soil_end$`Fit_soil[2280, 402:801]`[120:180])
C_4<-mean(Fit_soil_end$`Fit_soil[2280, 402:801]`[180:240])
C_5<-mean(Fit_soil_end$`Fit_soil[2280, 402:801]`[240:300])


inl<-ggplot(obs_data_in, aes(x=time1, y=conc1)) + 
  geom_errorbar(aes(ymin= conc1-sd1, ymax=conc1+sd1), width=0.01) +
  geom_point()+
  geom_line(data = data.frame(Fit1_in), aes(x=time,y=conc))+
  xlab("Time (minutes)") +
  ylab("C/Co(-)")+
  theme_bw() +
  coord_cartesian(xlim = c(0.0, 250.0),ylim = c(0,1))

out<-ggplot(data = obs_data_out, aes(x=time, y=conc)) + 
  geom_errorbar(aes(ymin=conc-err, ymax=conc+err), width=0.01) +
  geom_point()+
  geom_line(data = data.frame(Fit1),aes(x=time,y=conc))+
  xlab("Time (minutes)") +
  ylab("C/Co(-)")+
  theme_bw()+
  coord_cartesian(xlim = c(0, 250),ylim = c(0,1))

comb_fig<- ggarrange ( inl, out,
                       ncol = 2, nrow = 1,
                       labels= c ("a.","b.")
)
comb_fig

en_sim<-unclass(as.POSIXlt(Sys.time()))

print(paste0("simulation time= ",
             (((en_sim$hour*60)+(en_sim$min)+(en_sim$sec/60))-
                ((st_sim$hour*60)+(st_sim$min)+(st_sim$sec/60))),
             "minutes"))
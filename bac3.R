rm(list=ls())
setwd("C:\\Users\\achandra\\Desktop\\APARNA1\\codes\\Experiments\\Exp_bac1")
library(devtools)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
library(deSolve)
library(rootSolve)
library(ReacTran)
library(rodeo)
library(readxl)
library(FME)
library(coda)


st_sim<-unclass(as.POSIXlt(Sys.time()))
parm_2 <- as.data.frame(read_excel("definitions.xlsx", 
                                   sheet = "pars",
                                   col_names = TRUE,
                                   range="A1:F2"))

parm_guess<-c(10^-3)

vars <- as.data.frame(read_excel("definitions.xlsx", 
                                 sheet = "vars", 
                                 col_names = TRUE,  
                                 skip = 0))
pars <- as.data.frame(read_excel("definitions.xlsx", 
                                 sheet = "pars", 
                                 col_names = TRUE,  
                                 skip = 0))
funs <- as.data.frame(read_excel("definitions.xlsx", 
                                 sheet = "funcs",
                                 col_names = TRUE,  
                                 skip = 0))
pros <- as.data.frame(read_excel("definitions.xlsx", 
                                 sheet = "pros", 
                                 col_names = TRUE,  
                                 skip = 0))
stoi <- as.data.frame(read_excel("definitions.xlsx", 
                                 sheet = "stoi", 
                                 col_names = TRUE,  
                                 skip = 0))

Data<- as.data.frame(read_excel(path= "data_exp.xlsx",
                                sheet="exp_3", 
                                col_names=TRUE
                                )
                     )
##########################################
######FUNCTIONS AND INITIALIZATION########
##########################################

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
  vars <- as.data.frame(read_excel("definitions.xlsx", sheet = "vars", col_names = TRUE,  skip = 0))
  pars <- as.data.frame(read_excel("definitions.xlsx", sheet = "pars", col_names = TRUE,  skip = 0))
  funs <- as.data.frame(read_excel("definitions.xlsx", sheet = "funcs",col_names = TRUE,  skip = 0))
  pros <- as.data.frame(read_excel("definitions.xlsx", sheet = "pros", col_names = TRUE,  skip = 0))
  stoi <- as.data.frame(read_excel("definitions.xlsx", sheet = "stoi", col_names = TRUE,  skip = 0))
  
  model <- rodeo$new(
    vars = vars[1:3],
    pars = pars[1:3],
    funs= funs,
    pros = pros,
    stoi = asMatrix(stoi),
    asMatrix=TRUE,
    dim= c(nLayers)
  )
  
  model$compile("ADE3.f95")
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
Fit<-nls.lm(par=parm_guess, lower=0, upper=5, fn=Model_cost)
Fit<-bobyqa(x0=parm_guess, fn=Model_fun, lower = 0, upper = Inf)
Fit1_1 <- modFit(f = Model_cost,
                 p = parm_guess,
                 method = 'bobyqa',
                 lower = rep(0,1),
                 upper = rep(1,1)
)
Fit1_1plot<-Model_fun(Fit1_1$par)

MC_results3<- modMCMC(f=Model_cost,
                      p=parm_guess,
                      updatecov = 15,
                      lower = rep(0,2),
                      upper = rep(1,2),
                      niter = 5000)
save.image(file = "D:\\APARNA1\\codes\\Experiments\\results\\fast\\noDOnoN_Ffast__combineddata_con(0,1)_withsoil_AM15DR20_26thFeb(1-5000).RData")
MC_results3<- modMCMC(f=Model_cost,
                      p=MC_results3$bestpar,
                      updatecov = 15,
                      ntrydr = 5,
                      lower = rep(0,2),
                      upper = rep(1,2),
                      niter = 5000)
save.image(file = "D:\\APARNA1\\codes\\Experiments\\results\\fast\\noDOnoN_Ffast__combineddata_con(0,1)_withsoil_AM15DR20_26thFeb(5000-10000).RData")

Fit2<-Model_fun(MC_results3$bestpar)

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
  vars <- as.data.frame(read_excel("definitions.xlsx", sheet = "vars", col_names = TRUE,  skip = 0))
  pars <- as.data.frame(read_excel("definitions.xlsx", sheet = "pars", col_names = TRUE,  skip = 0))
  funs <- as.data.frame(read_excel("definitions.xlsx", sheet = "funcs",col_names = TRUE,  skip = 0))
  pros <- as.data.frame(read_excel("definitions.xlsx", sheet = "pros", col_names = TRUE,  skip = 0))
  stoi <- as.data.frame(read_excel("definitions.xlsx", sheet = "stoi", col_names = TRUE,  skip = 0))
  
  model <- rodeo$new(
    vars = vars[1:3],
    pars = pars[1:3],
    funs= funs,
    pros = pros,
    stoi = asMatrix(stoi),
    asMatrix=TRUE,
    dim= c(nLayers)
  )
  
  model$compile("ADE3.f95")
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
write.csv(x=Fit_soil,file="D:\\APARNA1\\codes\\Experiments\\results\\fast\\soil_noDOnoN_t1_unconfined")
en_sim<-unclass(as.POSIXlt(Sys.time()))

print(paste0("simulation time= ",
             (((en_sim$hour*60)+(en_sim$min)+(en_sim$sec/60))-
                ((st_sim$hour*60)+(st_sim$min)+(st_sim$sec/60))),
             "minutes")
)
Data1<- as.data.frame(read_excel(path= "data_exp.xlsx",
                                 sheet="exp_3", 
                                 col_names=TRUE))
# Data1<-cbind(Data1$time,Data1$val)
# colnames(Data1)<-c("time","conc")

SSR1<-modCost(obs=Data1,model = Fit2_1)

plot(Data1,pch="*",ylim=c(0,0.4),xlab="time (mins)",ylab="C/Co(-)")
lines(Fit1_1plot,lwd="3")
lines(Fit0_1,col="red")


save.image(file = "D:\\APARNA1\\codes\\Experiments\\results\\fast\\noDOnoN_Ffast__combineddata_con(0,1)_withsoil_AM15DR20_21thFeb.RData")


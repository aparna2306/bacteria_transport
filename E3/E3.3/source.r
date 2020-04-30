rm(list=ls())
library(ReacTran)
library(deSolve)
library(rootSolve)
library(readxl)
library("rodeo")
library(FME)
library(coda)
library(ggplot2)
library(ggpubr)
#==================================================#
#       INITIALIZE TIME AND LOAD DATA FUNCTIONS    #
#==================================================#

st_sim<-unclass(as.POSIXlt(Sys.time()))

parm_guess1<-c(0.00001)
obs_data_out<-read_excel("observation_data.xlsx",
                 sheet = "F1",
                 col_names = TRUE,
                 skip = 0)

obs_data_in<- read_excel("observation_data.xlsx",
                         sheet = "F1_in",
                         col_names = TRUE,
                         skip = 0)
###########################################################
# DETERMINING DECAY RATE USING INLET CONCENTRATION CURVES #
###########################################################

parmguess_in <- c(mu_r=0.02)
Model_fun<-function(p){
  derivs<-function(t,state,p){
    with(as.list(c(state,p)),{
      dC<- -mu_r*C
      return(list(c(dC)))
    })
  }
  state<- c(C=1)
  times <- seq(from=0,to=300,by=5)
  t<-Data$time
  out<- ode(y=state,times=times,func = derivs,parms = p)
  return(out)
}
Model_Cost<-function(p){
  out_df<-Model_fun(p)
  return(modCost(model=out_df,obs=obs_data_in))
}
Fit1_in<-modFit(f = Model_Cost, p = parmguess_in ,method = 'bobyqa')

#LOAD DATA
parm <- as.data.frame(read_excel("dparameters.xlsx", 
                                 sheet = "pars",
                                 range = "A1:F3",
                                 col_names = TRUE))

parm_guess<-parm$values
#====================================#
# RODEO SOLUTION- FUNCTIONS OUTLET   #
#====================================#

Model_fun1<-function(p1){
asMatrix<- function(df){
  matrix(unlist(df[,2:ncol(df)]), 
         nrow=nrow(df), 
         ncol= ncol(df)-1,
         dimnames =list(df[,1], names(df)[2:ncol(df)]))
}
#define grid
L <- 20
nLayers <- 200
dx <- L/nLayers
#Extract data for rodeo
vars <- as.data.frame(read_excel("parameters.xlsx", sheet = "vars", col_names = TRUE,  skip = 0))
pars <- as.data.frame(read_excel("parameters.xlsx", sheet = "pars", col_names = TRUE,  skip = 0))
funs <- as.data.frame(read_excel("parameters.xlsx", sheet = "funcs", col_names = TRUE,  skip = 0))
pros <- as.data.frame(read_excel("parameters.xlsx", sheet = "pros", col_names = TRUE,  skip = 0))
stoi <- as.data.frame(read_excel("parameters.xlsx", sheet = "stoi", col_names = TRUE,  skip = 0))

model <- rodeo$new(
  vars = vars[1:3],
  pars = pars[1:3],
  funs = funs,
  pros = pros,
  stoi = asMatrix(stoi),
  asMatrix=TRUE,
  dim= c(nLayers)
)

#Generate function to return derivatives
model$compile("boundary.f95")

#duplicating the variables to grid
tmp_var<-vars$values
names(tmp_var)<- vars$name
tmp_var <- matrix (rep(tmp_var,
                       each=nLayers),
                   nrow=nLayers, 
                   ncol=length(tmp_var),
                   dimnames=list(NULL,
                                 names(tmp_var)))

model$setVars(tmp_var)

#Duplicating parameters to the grid
tmp_par<-c(as.numeric(p1),pars$values_1[2:10])
names(tmp_par)<-pars$name
tmp_par <- matrix(rep(tmp_par,
                  each=nLayers), 
              nrow=nLayers, 
              ncol=length(tmp_par),
              dimnames=list(NULL,names(tmp_par)))
tmp_par[1,"leftmost"]<-1
model$setPars(tmp_par)
times<-seq(from=0,to=300,by=1)
out<- model$dynamics(times=times)
out_df<-(cbind(out[,"time"],out[,"C.165"]))
colnames(out_df)=c("time","conc")
return(out_df)
}

#====================================#
# GENERATING MODEL COST FUNCTION     #
#====================================#

Model_cost1<-function(p1){
  out_df<-Model_fun1(p1)
  Data<-as.data.frame(read_excel("observation_data.xlsx",
                                  sheet = "F1",
                                  col_names = TRUE,
                                  skip = 0))
  cost_fn<-modCost(model=out_df, obs=Data)
return(cost_fn)
}


#====================================#
# MODEL FITTING USING MODFIT         #
#====================================#
Fit<- modFit(f= Model_cost1,
              p = parm_guess1,
              lower = rep(0,1),
              upper = rep(2,1),
              method = 'bobyqa'
              )

Fit1<-Model_fun1(Fit1$par)
summary(Fit1)
#############################################
######sediment retained data concentration###
#############################################

Model_fun_soil<-function(p1){
  asMatrix<- function(df){
    matrix(unlist(df[,2:ncol(df)]),
           nrow=nrow(df),
           ncol= ncol(df)-1,
           dimnames =list(df[,1], names(df)[2:ncol(df)]))
  }
  #define grid
  L <- 20
  nLayers <- 300
  dx <- L/nLayers
  #Extract data for rodeo
  vars <- as.data.frame(read_excel("parameters.xlsx", sheet = "vars", col_names = TRUE,  skip = 0))
  pars <- as.data.frame(read_excel("parameters.xlsx", sheet = "pars", col_names = TRUE,  skip = 0))
  funs <- as.data.frame(read_excel("parameters.xlsx", sheet = "funcs", col_names = TRUE,  skip = 0))
  pros <- as.data.frame(read_excel("parameters.xlsx", sheet = "pros", col_names = TRUE,  skip = 0))
  stoi <- as.data.frame(read_excel("parameters.xlsx", sheet = "stoi", col_names = TRUE,  skip = 0))

  model <- rodeo$new(
    vars = vars[1:3],
    pars = pars[1:3],
    funs = funs,
    pros = pros,
    stoi = asMatrix(stoi),
    asMatrix=TRUE,
    dim= c(nLayers)
  )

  #Generate function to return derivatives
  model$compile("boundary.f95")

  #duplicating the variables to grid
  tmp_var<-vars$values
  names(tmp_var)<- vars$name
  tmp_var <- matrix (rep(tmp_var,
                         each=nLayers),
                     nrow=nLayers,
                     ncol=length(tmp_var),
                     dimnames=list(NULL,
                                   names(tmp_var)))

  model$setVars(tmp_var)

  #Duplicating parameters to the grid
  tmp_par<-c(as.numeric(p1),pars$values_1[2:10])
  names(tmp_par)<-pars$name
  tmp_par <- matrix(rep(tmp_par,
                        each=nLayers),
                    nrow=nLayers,
                    ncol=length(tmp_par),
                    dimnames=list(NULL,names(tmp_par)))
  tmp_par[1,"leftmost"]<-1
  model$setPars(tmp_par)

  #define time and compile
  times <- seq(0,250,by=0.1)
  out<- model$dynamics(times=times)
  return(out)
 }

Fit_soil<-Model_fun_soil(Fit1$par)

Fit_soil_end<-as.data.frame(Fit_soil[2300,302:601])
C_1<-mean(Fit_soil_end$`Fit_soil[2300, 302:601]`[1:60])
C_2<-mean(Fit_soil_end$`Fit_soil[2300, 302:601]`[60:120])
C_3<-mean(Fit_soil_end$`Fit_soil[2300, 302:601]`[120:180])
C_4<-mean(Fit_soil_end$`Fit_soil[2300, 302:601]`[180:240])
C_5<-mean(Fit_soil_end$`Fit_soil[2300, 302:601]`[240:300])

# #====================================#
# # PLOTTING                           #
# #====================================#

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
             (
               ((en_sim$hour*60)+(en_sim$min)+(en_sim$sec/60))-
                 ((st_sim$hour*60)+(st_sim$min)+st_sim$sec/60)
             ),
             "minutes"))

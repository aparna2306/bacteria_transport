#Solution to the Advective Dispersive equation by the rodeo method and compare this with the analytical solution
rm(list=ls())
setwd("C:/Users/achandra/Desktop/APARNA1/codes/Experiments/batch")
library(ReacTran)
library(deSolve)
library(rootSolve)
library(readxl)
library("rodeo")
library(FME)

#==================================================#
#       INITIALIZE TIME AND LOAD DATA FUNCTIONS    #
#==================================================#
st_sim<-unclass(as.POSIXlt(Sys.time()))

#LOAD DATA
# sol1<-read_excel("paper_data.xlsx",
#                 sheet = "F1",
#                 col_names = TRUE,
#                 skip = 0)

#LOAD PARAMTER GUESS, EXCLUDING BOUNDARY CONDITION
# parm <- as.data.frame(read_excel("definitions_col1_o2.xlsx", 
#                                  sheet = "pars",
#                                  range = "A1:F3",
#                                  col_names = TRUE))
parms <- c(mu_d=0.02)
parms1<-c(mu_r=0.01)
Data<-as.data.frame(read_excel(path = 'inlet_noDOnoN.xlsx',
                               sheet = 'Sheet3',
                               col_names=TRUE,
                               skip=0)
)
#====================================#
# RODEO SOLUTION- FUNCTIONS          #
#====================================#

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
 plot<-(Model_fun(parms))
# Model_fun1<-function(p){
#   derivs<-function(t,state,p){
#     with(as.list(c(state,p)),{
#       dC<- -mu_r*C
#       return(list(c(dC)))
#     })
#   }
#   state<- c(C=1)
#   times <- seq(0,1200,by=0.1)
#   out<- ode(y=state,times=times,func = derivs,parms = p)
#   return(out)
# }
Data<-as.data.frame(read_excel(path = 'C:/Users/achandra/Desktop/APARNA1/codes/Experiments/batch/inlet_noDOnoN.xlsx',
                               sheet = 'Sheet3',
                               col_names=TRUE,
                               skip=0)
)
# Data1<-as.data.frame(read_excel(path = 'C:/Users/achandra/Desktop/APARNA1/codes/Experiments/batch/inlet_noDOnoN.xlsx',
#                                sheet = 'Sheet2',
#                                 col_names=TRUE,
#                                skip=0)
# )

Model_Cost<-function(p){
 out_df<-Model_fun(p)
 return(modCost(model=out_df,obs=Data))
}

# Model_Cost1<-function(p){
#   out_df<-Model_fun1(p)
#   return(modCost(model = out_df,obs = Data1))
# }
#====================================#
# MODEL FITTING USING MODFIT         #
#====================================#

Fit<-modFit(f = Model_Cost, p = parms,method = 'bobyqa')
MC_fit<-modMCMC(f = Model_Cost, updatecov = 150, p = Fit$par, niter = 10000)

# Fit1 <- modFit(f = Model_Cost1,p = parms1)
# MC_fit1<-modMCMC(f = Model_Cost1,p = Fit1$par, niter = 10000)
# Fit0<-Model_fun(parms)
# Fit01 <-Model_fun1(parms1)
Fit_1<-Model_fun(Fit$par)
# Fit1_1<-Model_fun1(Fit1$par)
Fit_2<-Model_fun(MC_fit$bestpar)
plot(MC_fit)
plot(Data)
lines(Fit_1)
lines(Fit0,col='red')

plot(Data1,ylim=c(0,2))
lines(Fit1_1)
lines(Fit01,col='red')

print(Fit)
print(Fit1)

# Snsbac<-sensFun(func=Model_fun1,parms = pars, sensvar = 'C')
# cor(Snsbac)
# Coll<-collin(Snsbac)


# asMatrix<- function(df){
#   matrix(unlist(df[,2:ncol(df)]), 
#          nrow=nrow(df), 
#          ncol= ncol(df)-1,
#          dimnames =list(df[,1], names(df)[2:ncol(df)]))
# }
#define grid
# L <- 20
# nLayers <- 400
# dx <- L/nLayers
#Extract data for rodeo
# vars <- as.data.frame(read_excel("definitions_col1_o2.xlsx", sheet = "vars", col_names = TRUE,  skip = 0))
# pars <- as.data.frame(read_excel("definitions_col1_o2.xlsx", sheet = "pars", col_names = TRUE,  skip = 0))
# funs <- as.data.frame(read_excel("definitions_col1_o2.xlsx", sheet = "funcs", col_names = TRUE,  skip = 0))
# pros <- as.data.frame(read_excel("definitions_col1_o2.xlsx", sheet = "pros", col_names = TRUE,  skip = 0))
# stoi <- as.data.frame(read_excel("definitions_col1_o2.xlsx", sheet = "stoi", col_names = TRUE,  skip = 0))

# model <- rodeo$new(
#   vars = vars[1:3],
#   pars = pars[1:3],
#   funs = funs,
#   pros = pros,
#   stoi = asMatrix(stoi),
#   asMatrix=TRUE,
# )

#Generate function to return derivatives
# model$compile("ADE1.f95")

#duplicating the variables to grid
# tmp_var<-vars$values
# names(tmp_var)<- vars$name
# tmp_var <- matrix (rep(tmp_var,
                   #     each=nLayers),
                   # nrow=nLayers, 
                   # ncol=length(tmp_var),
                   # dimnames=list(NULL,
                   #               names(tmp_var)))

# model$setVars(tmp_var)

#Duplicating parameters to the grid
# tmp_par<-c(as.numeric(p1))
# names(tmp_par)<-pars$name
# tmp_par <- matrix(rep(tmp_par,
              #     each=nLayers), 
              # nrow=nLayers, 
              # ncol=length(tmp_par),
              # dimnames=list(NULL,names(tmp_par)))
# tmp_par[1,"leftmost"]<-1
# model$setPars(tmp_par)
#define time and compile
  

# out_df<-(cbind(out[,"time"],out[,"C.200"]))
# colnames(out_df)=c("time","conc")


#====================================#
# INITIAL FITTING                    #
#====================================#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
# STOP HERE FOR FITTING JUST MODEL TO DATA WITHOUT INVERSE PARAMTERIZATION  #
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
# Fit0_1<-Model_fun1(parm_guess1)


#====================================#
# GENERATING MODEL COST FUNCTION     #
#====================================#

# Model_cost1<-function(p1){
# 
#   out_df<-Model_fun1(p1)
#   Data<- as.data.frame(read_excel("paper_data.xlsx",
#                                   sheet = "F1",
#                                   col_names = TRUE,
#                                   skip = 0))
#   cost_fn<-modCost(model=out_df,obs=Data)
# }
# 
# #====================================#
# # MODEL FITTING USING MODFIT         #
# #====================================#
# Fit1 <- modFit(f = Model_cost1,
#               p = parm_guess1,
#               lower = rep(0,2),
#               upper = rep(Inf,2))
# 
# Fit1_1<-Model_fun1(Fit1$par)
# 
# #====================================#
# # MODEL FITTING USING MODMCMC        #
# #====================================#
# MC_results1<- modMCMC(f=Model_cost1,
#                      p=Fit1$par,
#                      lower = rep(0,2),
#                      upper = rep(Inf,2),
#                      niter = 20000)
# Fit2_1<-Model_fun1(MC_results1$bestpar)
# write.csv(Fit2_1,file = "F:\\codes\\Experiments\\results\\fast\\withDOwithN_fast_conf_t1")
# 
# Model_fun_soil<-function(p1){
#   asMatrix<- function(df){
#     matrix(unlist(df[,2:ncol(df)]), 
#            nrow=nrow(df), 
#            ncol= ncol(df)-1,
#            dimnames =list(df[,1], names(df)[2:ncol(df)]))
#   }
#   #define grid
#   L <- 20
#   nLayers <- 400
#   dx <- L/nLayers
#   #Extract data for rodeo
#   vars <- as.data.frame(read_excel("definitions_attdet_col1.xlsx", sheet = "vars", col_names = TRUE,  skip = 0))
#   pars <- as.data.frame(read_excel("definitions_attdet_col1.xlsx", sheet = "pars", col_names = TRUE,  skip = 0))
#   funs <- as.data.frame(read_excel("definitions_attdet_col1.xlsx", sheet = "funcs", col_names = TRUE,  skip = 0))
#   pros <- as.data.frame(read_excel("definitions_attdet_col1.xlsx", sheet = "pros", col_names = TRUE,  skip = 0))
#   stoi <- as.data.frame(read_excel("definitions_attdet_col1.xlsx", sheet = "stoi", col_names = TRUE,  skip = 0))
#   
#   model <- rodeo$new(
#     vars = vars[1:3],
#     pars = pars[1:3],
#     funs = funs,
#     pros = pros,
#     stoi = asMatrix(stoi),
#     asMatrix=TRUE,
#     dim= c(nLayers)
#   )
#   
#   #Generate function to return derivatives
#   model$compile("ADE1.f95")
#   
#   #duplicating the variables to grid
#   tmp_var<-vars$values
#   names(tmp_var)<- vars$name
#   tmp_var <- matrix (rep(tmp_var,
#                          each=nLayers),
#                      nrow=nLayers, 
#                      ncol=length(tmp_var),
#                      dimnames=list(NULL,
#                                    names(tmp_var)))
#   
#   model$setVars(tmp_var)
#   
#   #Duplicating parameters to the grid
#   tmp_par<-c(as.numeric(p1),pars$values_1[4:10])
#   names(tmp_par)<-pars$name
#   tmp_par <- matrix(rep(tmp_par,
#                         each=nLayers), 
#                     nrow=nLayers, 
#                     ncol=length(tmp_par),
#                     dimnames=list(NULL,names(tmp_par)))
#   tmp_par[1,"leftmost"]<-1
#   model$setPars(tmp_par)
#   
#   #define time and compile
#   times <- seq(0,350,by=0.1)
#   out<- model$dynamics(times=times)
#   return(out)
# }
# 
# Fit_soil<-Model_fun_soil(MC_results1$bestpar)
# Fit_soil_end<-as.data.frame(Fit_soil[2308,401:801])
# C_1<-mean(Fit_soil_end$`Fit_soil[2308, 401:801]`[1:60])
# C_2<-mean(Fit_soil_end$`Fit_soil[2308, 401:801]`[60:120])
# C_3<-mean(Fit_soil_end$`Fit_soil[2308, 401:801]`[120:180])
# C_4<-mean(Fit_soil_end$`Fit_soil[2308, 401:801]`[180:240])
# C_5<-mean(Fit_soil_end$`Fit_soil[2308, 401:801]`[240:300])
# 
# 
# 
# 
# #====================================#
# # SENSITIVITY ANALYSIS               #
# #====================================#
# 
# #sR <- sensRange(func = Model_fun, parms = MC_results$bestpar, parInput = MC_results$pars)
# 
# #====================================#
# # PLOTTING                           #
# #====================================#
# plot(sol1,pch="*")
# lines(Fit0_1,col="dark red")
# lines(Fit1_1,col="red")
# lines(Fit2_1,col="black")
# SSR<-modCost(model = Fit2_1,obs=sol1)
# en_sim<-unclass(as.POSIXlt(Sys.time()))
# print(paste0("simulation time= ",
#              (
#                ((en_sim$hour*60)+(en_sim$min)+(en_sim$sec/60))-
#                  ((st_sim$hour*60)+(st_sim$min)+st_sim$sec/60)
#              ),
#              "minutes"))
# save.image(file="F:/codes/Experiments/results/combined data_withDOwithN_Ffast_t1.RData")
# 

#====================================#
# ANALYTICAL SOLUTION-IMPORT VALUES  #
#====================================#


#====================================#
# ERRORS                             #
#====================================#


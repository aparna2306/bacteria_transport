rm(list=ls())
library(ReacTran)
library(deSolve)
library(rootSolve)
library(readxl)
library("rodeo")
library(FME)
library(ggplot2)
#==================================================#
#       INITIALIZE TIME AND LOAD DATA FUNCTIONS    #
#==================================================#
st_sim<-unclass(as.POSIXlt(Sys.time()))

#LOAD DATA from experimental observations 
obs_data_out<-read_excel("obs_data.xlsx",
                sheet = "obs_data",
                col_names = TRUE,
                skip = 0)

#LOAD PARAMTER GUESS, EXCLUDING BOUNDARY CONDITION
parm <- as.data.frame(read_excel("parameters.xlsx", 
                                 sheet = "pars",
                                 range = "A1:D3",
                                 col_names = TRUE))
parm_guess<-(parm$values)
parm_guess_tracer<-c(0,1)
names(parm_guess)<-parm$name

#====================================#
# RODEO SOLUTION- FUNCTIONS          #
#====================================#

Model_fun<-function(p){
asMatrix<- function(df){
  matrix(unlist(df[,2:ncol(df)]), 
         nrow=nrow(df), 
         ncol= ncol(df)-1,
         dimnames =list(df[,1], names(df)[2:ncol(df)]))
}
#define grid
L <- 15*2
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
tmp_par<-c(as.numeric(p),as.numeric(pars$values[3:7]))
names(tmp_par)<-pars$name
tmp_par <- matrix(rep(tmp_par,
                  each=nLayers), 
              nrow=nLayers, 
              ncol=length(tmp_par),
              dimnames=list(NULL,names(tmp_par)))
tmp_par[1,"leftmost"]<-1
model$setPars(tmp_par)

#define time and compile
times <- seq(0,7,by=0.01)
out<- model$dynamics(times=times)
out_df<-(cbind(out[,"time"],out[,"C.150"]))
colnames(out_df)=c("time","conc")
return(out_df)
}

#====================================#
# FITTING USING INITIAL GUESS        #
#====================================#

Fit_0<-Model_fun(parm_guess)


#====================================#
# GENERATING MODEL COST FUNCTION     #
#====================================#

Model_cost<-function(p){

  out_df<-Model_fun(p)
  Data<- as.data.frame(read_excel("obs_data.xlsx",
                                  sheet = "obs_data",
                                  col_names = TRUE,
                                  skip = 0))
  cost_fn<-modCost(model=out_df,obs=Data)
}
#====================================#
# MODEL FITTING USING MODFIT         #
#====================================#
Fit <- modFit(f = Model_cost,
              p = parm_guess,
              lower = rep(0,2),
              upper = rep(Inf,2))

Fit_1<-Model_fun(Fit$par)

#====================================#
# MODEL FITTING USING MODMCMC        #
#====================================#
MC_results<- modMCMC(f=Model_cost,
                     p=Fit$par,
                     lower = rep(0,2),
                     upper=c(1,10),
                     updatecov = 100,
                     niter = 10000)

Fit_2<-Model_fun(MC_results$bestpar)

#====================================#
# PLOTTING                           #
#====================================#
out1<-ggplot(data = obs_data_out, aes(x=time, y=conc)) + 
  geom_point()+
  geom_line(data = data.frame(Fit2),aes(x=time,y=conc))+
  xlab("Time (minutes)") +
  ylab("C/Co(-)")+
  theme_bw()+
  coord_cartesian(xlim = c(0, 250),ylim = c(0,1))

en_sim<-unclass(as.POSIXlt(Sys.time()))
print(paste0("simulation time= ",
             (
               ((en_sim$hour*60)+(en_sim$min)+(en_sim$sec/60))-
                 ((st_sim$hour*60)+(st_sim$min)+st_sim$sec/60)
             ),
             "minutes"))
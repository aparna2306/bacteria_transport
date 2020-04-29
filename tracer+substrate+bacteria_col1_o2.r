#Solution to the Advective Dispersive equation by the rodeo method and compare this with the analytical solution
rm(list=ls())
setwd("C:/Users/achandra/Desktop/APARNA1/codes/Experiments/column/tracer+substrate+bacteria_column/col1_o2")
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
sol1<-read_excel("paper_data.xlsx",
                sheet = "F1",
                col_names = TRUE,
                skip = 0)

#LOAD PARAMTER GUESS, EXCLUDING BOUNDARY CONDITION
parm <- as.data.frame(read_excel("definitions_col1_o2.xlsx", 
                                 sheet = "pars",
                                 range = "A1:F3",
                                 col_names = TRUE))
parm_guess1<-parm$values_1

#====================================#
# RODEO SOLUTION- FUNCTIONS          #
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
nLayers <- 400
dx <- L/nLayers
#Extract data for rodeo
vars <- as.data.frame(read_excel("definitions_col1_o2.xlsx", sheet = "vars", col_names = TRUE,  skip = 0))
pars <- as.data.frame(read_excel("definitions_col1_o2.xlsx", sheet = "pars", col_names = TRUE,  skip = 0))
funs <- as.data.frame(read_excel("definitions_col1_o2.xlsx", sheet = "funcs", col_names = TRUE,  skip = 0))
pros <- as.data.frame(read_excel("definitions_col1_o2.xlsx", sheet = "pros", col_names = TRUE,  skip = 0))
stoi <- as.data.frame(read_excel("definitions_col1_o2.xlsx", sheet = "stoi", col_names = TRUE,  skip = 0))

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
model$compile("ADE1.f95")

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
tmp_par<-c(as.numeric(p1),pars$values_1[3:12])
names(tmp_par)<-pars$name
tmp_par <- matrix(rep(tmp_par,
                  each=nLayers), 
              nrow=nLayers, 
              ncol=length(tmp_par),
              dimnames=list(NULL,names(tmp_par)))
tmp_par[1,"leftmost"]<-1
model$setPars(tmp_par)
#define time and compile
times <- seq(0,240,by=0.5)
out<- model$dynamics(times=times)
out_df<-(cbind(out[,"time"],out[,"C.200"]))
colnames(out_df)=c("time","conc")
return(out_df)
}

#====================================#
# INITIAL FITTING                    #
#====================================#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
# STOP HERE FOR FITTING JUST MODEL TO DATA WITHOUT INVERSE PARAMTERIZATION  #
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
Fit0_1<-Model_fun1(parm_guess1)


#====================================#
# GENERATING MODEL COST FUNCTION     #
#====================================#

Model_cost1<-function(p1){

  out_df<-Model_fun1(p1)
  Data<- as.data.frame(read_excel("paper_data.xlsx",
                                  sheet = "F1",
                                  col_names = TRUE,
                                  skip = 0))
  cost_fn<-modCost(model=out_df,obs=Data)
}

#====================================#
# MODEL FITTING USING MODFIT         #
#====================================#
Fit1 <- modFit(f = Model_cost1,
              p = parm_guess1,
              method = 'bobyqa',
              lower = rep(0,2),
              upper = rep(20,2))

Fit1_1<-Model_fun1(Fit1$par)

#====================================#
# MODEL FITTING USING MODMCMC        #
#====================================#
MC_results1<- modMCMC(f=Model_cost1,
                     p=Fit1$par,
                     lower = rep(0,2),
                     upper = rep(Inf,2),
                     niter = 20000)
Fit2_1<-Model_fun1(MC_results1$bestpar)
write.csv(Fit2_1,file = "F:\\codes\\Experiments\\results\\fast\\withDOwithN_fast_conf_t1")

Model_fun_soil<-function(p1){
  asMatrix<- function(df){
    matrix(unlist(df[,2:ncol(df)]), 
           nrow=nrow(df), 
           ncol= ncol(df)-1,
           dimnames =list(df[,1], names(df)[2:ncol(df)]))
  }
  #define grid
  L <- 20
  nLayers <- 400
  dx <- L/nLayers
  #Extract data for rodeo
  vars <- as.data.frame(read_excel("definitions_col1_o2.xlsx", sheet = "vars", col_names = TRUE,  skip = 0))
  pars <- as.data.frame(read_excel("definitions_col1_o2.xlsx", sheet = "pars", col_names = TRUE,  skip = 0))
  funs <- as.data.frame(read_excel("definitions_col1_o2.xlsx", sheet = "funcs", col_names = TRUE,  skip = 0))
  pros <- as.data.frame(read_excel("definitions_col1_o2.xlsx", sheet = "pros", col_names = TRUE,  skip = 0))
  stoi <- as.data.frame(read_excel("definitions_col1_o2.xlsx", sheet = "stoi", col_names = TRUE,  skip = 0))
  
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
  model$compile("ADE1.f95")
  
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
  tmp_par<-c(c(0.00848,0),pars$values_1[3:12])
  names(tmp_par)<-pars$name
  tmp_par <- matrix(rep(tmp_par,
                        each=nLayers), 
                    nrow=nLayers, 
                    ncol=length(tmp_par),
                    dimnames=list(NULL,names(tmp_par)))
  tmp_par[1,"leftmost"]<-1
  model$setPars(tmp_par)
  
  #define time and compile
  times <- seq(0,350,by=0.1)
  out<- model$dynamics(times=times)
  return(out)
}

Fit_soil<-Model_fun_soil(Fit1$par)
Fit_soil_end<-as.data.frame(Fit_soil[2308,402:801])
C_1<-mean(Fit_soil_end$`Fit_soil[2308, 402:801]`[1:60])
C_2<-mean(Fit_soil_end$`Fit_soil[2308, 402:801]`[60:120])
C_3<-mean(Fit_soil_end$`Fit_soil[2308, 402:801]`[120:180])
C_4<-mean(Fit_soil_end$`Fit_soil[2308, 402:801]`[180:240])
C_5<-mean(Fit_soil_end$`Fit_soil[2308, 402:801]`[240:300])




#====================================#
# SENSITIVITY ANALYSIS               #
#====================================#

#sR <- sensRange(func = Model_fun, parms = MC_results$bestpar, parInput = MC_results$pars)

#====================================#
# PLOTTING                           #
#====================================#
plot(sol1,pch="*")
lines(Fit0_1,col="dark red")
lines(Fit1_1,col="red")
lines(Fit2_1,col="black")
SSR<-modCost(model = Fit2_1,obs=sol1)
en_sim<-unclass(as.POSIXlt(Sys.time()))
print(paste0("simulation time= ",
             (
               ((en_sim$hour*60)+(en_sim$min)+(en_sim$sec/60))-
                 ((st_sim$hour*60)+(st_sim$min)+st_sim$sec/60)
             ),
             "minutes"))
save.image(file="F:/codes/Experiments/results/combined data_withDOwithN_Ffast_t1.RData")


#====================================#
# ANALYTICAL SOLUTION-IMPORT VALUES  #
#====================================#


#====================================#
# ERRORS                             #
#====================================#


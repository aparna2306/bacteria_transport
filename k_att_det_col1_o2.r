#Solution to the Advective Dispersive equation by the rodeo method and compare this with the analytical solution
rm(list=ls())
# Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
setwd("C:/users/achandra/Desktop//APARNA1/codes/Experiments/column/bacteria_attachement_detachment/Col1_o2")
library(ReacTran)
library(deSolve)
library(rootSolve)
library(readxl)
library("rodeo")
library(FME)
library(coda)
#==================================================#
#       INITIALIZE TIME AND LOAD DATA FUNCTIONS    #
#==================================================#
st_sim<-unclass(as.POSIXlt(Sys.time()))

#LOAD DATA
#LOAD PARAMTER GUESS, EXCLUDING BOUNDARY CONDITION
parm <- as.data.frame(read_excel("definitions_attdet_col1.xlsx", 
                                 sheet = "pars",
                                 range = "A1:F3",
                                 col_names = TRUE))
parm_guess1<-c(0.00001)
sol1<-read_excel("observation_data.xlsx",
                 sheet = "F1",
                 col_names = TRUE,
                 skip = 0)
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
nLayers <- 200
dx <- L/nLayers
#Extract data for rodeo
vars <- as.data.frame(read_excel("definitions_attdet_col1.xlsx", sheet = "vars", col_names = TRUE,  skip = 0))
pars <- as.data.frame(read_excel("definitions_attdet_col1.xlsx", sheet = "pars", col_names = TRUE,  skip = 0))
funs <- as.data.frame(read_excel("definitions_attdet_col1.xlsx", sheet = "funcs", col_names = TRUE,  skip = 0))
pros <- as.data.frame(read_excel("definitions_attdet_col1.xlsx", sheet = "pros", col_names = TRUE,  skip = 0))
stoi <- as.data.frame(read_excel("definitions_attdet_col1.xlsx", sheet = "stoi", col_names = TRUE,  skip = 0))

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
# times <- seq(0,250,by=0.5)
times<-sol1$time
out<- model$dynamics(times=times)
out_df<-(cbind(out[,"time"],out[,"C.165"]))
colnames(out_df)=c("time","conc")
return(out_df)
}

#====================================#
# INITIAL FITTING                    #
#====================================#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
# STOP HERE FOR FITTING JUST MODEL TO DATA WITHOUT INVERSE PARAMTERIZATION  #
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
fit<-Model_fun1(parm_guess1)

# sol1<-read_excel("observation_data.xlsx",
#                  sheet = "F1",
#                  col_names = TRUE,
#                  skip = 0)
#====================================#
# GENERATING MODEL COST FUNCTION     #
#====================================#

Model_cost1<-function(p1){
  
  out_df<-Model_fun1(p1)
  
  # Fit_soil_end<-as.data.frame(out[,202:401])
  # C_1<-rowMeans(Fit_soil_end[,1:33])
  # C_2<-rowMeans(Fit_soil_end[,33:66])
  # C_3<-rowMeans(Fit_soil_end[,66:99])
  # C_4<-rowMeans(Fit_soil_end[,99:132])
  # C_5<-rowMeans(Fit_soil_end[,132:165])
 
  # out_df<-as.data.frame(cbind(out[,"time"],out[,"C.165"],C_1,C_2,C_3,C_4,C_5))
  # out_df<-as.data.frame(cbind(out[,"time"],out[,"C.165"]))
  # colnames(out_df)<-c("time","conc")
 
  Data<-as.data.frame(read_excel("observation_data.xlsx",
                                  sheet = "F1",
                                  col_names = TRUE,
                                  skip = 0))
  # colnames(Data)<-c("time","conc",' C1',' C2',' C3',' C4',' C5')
  
  # Data_long<-cross2long(Data,x='time')
  
  # Data_long<-cbind(Data_l,Data2$std)
  # colnames(Data_long)<-c('name','time','val')
  cost_fn<-modCost(model=out_df, obs=Data)
 
return(cost_fn)
}


#====================================#
# MODEL FITTING USING MODFIT         #
#====================================#
Fit1<- modFit(f= Model_cost1,
              p = parm_guess1,
              lower = rep(0,1),
              upper = rep(2,1),
              method = 'bobyqa'
              )

Fit1_1<-Model_fun1(Fit1$par)

# sol1<-read_excel("observation_data.xlsx",
#                  sheet = "F1",
#                  col_names = TRUE,
#                  skip = 0)

# write.csv(Fit2_1,file = "D:\\APARNA1\\codes\\Experiments\\results\\fast\\DOnoN_fast_unconf_t1")

# SSR1<-modCost(obs=sol1,model=Fit1_1)



############################################################################################################################################

#====================================#
# MODEL FITTING USING MODMCMC        #
#====================================#
MC_results1<- modMCMC(f=Model_cost1,
                     p=parm_guess1,
                     updatecov = 15,
                     lower = rep(0,1),
                     upper = rep(Inf,1),
                     niter = 5000)
Fit2_1<-Model_fun1(MC_results1$bestpar)
############################################################################################################################################

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
  vars <- as.data.frame(read_excel("definitions_attdet_col1.xlsx", sheet = "vars", col_names = TRUE,  skip = 0))
  pars <- as.data.frame(read_excel("definitions_attdet_col1.xlsx", sheet = "pars", col_names = TRUE,  skip = 0))
  funs <- as.data.frame(read_excel("definitions_attdet_col1.xlsx", sheet = "funcs", col_names = TRUE,  skip = 0))
  pros <- as.data.frame(read_excel("definitions_attdet_col1.xlsx", sheet = "pros", col_names = TRUE,  skip = 0))
  stoi <- as.data.frame(read_excel("definitions_attdet_col1.xlsx", sheet = "stoi", col_names = TRUE,  skip = 0))

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
plot(x=sol1$time,y=sol1$conc,ylab="C/Co(-)",xlab="time (mins)",pch=16,col="red",ylim = c(0,1.5))
lines(Fit1_1,col="red",lwd="3")
lines(Fit1_0)
legend("topright", legend=c("observation", "model"),lty=c(0,1),lwd=c(0,3),pch=c(16,NA)) 

en_sim<-unclass(as.POSIXlt(Sys.time()))
print(paste0("simulation time= ",
             (
               ((en_sim$hour*60)+(en_sim$min)+(en_sim$sec/60))-
                 ((st_sim$hour*60)+(st_sim$min)+st_sim$sec/60)
             ),
             "minutes"))
save.image(file = "D:\\APARNA1\\codes\\Experiments\\results\\withDOnoN_fast_combineddata_modMCMC_iter20000_conf'(0,1)_AM20DR20_21Feb_AMDR.RData")

rm(list=ls())
setwd("C:/Users/achandra/Desktop/APARNA1/codes/Experiments/column/tracer+substrate+bacteria_column/col1")
library(deSolve)
library(rootSolve)
library(ReacTran)
library(rodeo)
library(readxl)
library(FME)
library(coda)

st_sim<-unclass(as.POSIXlt(Sys.time()))

parm <- as.data.frame(read_excel("definitions.xlsx", 
                                 sheet = "pars",
                                 col_names = TRUE,
                                 range="A1:F3"))
parm_guess1<-c(8.308736e-1,1.690628e-02)

##########################################
######FUNCTIONS AND INITIALIZATION########
##########################################
######## COL1 ########
Model_fun1<-function(p){
  
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
  
  model$compile("ADE1.f95")
  tmp_var<- vars$values
  names(tmp_var)<- vars$name
  tmp_var <- matrix (rep(tmp_var,
                         each=nLayers),
                     nrow=nLayers, 
                     ncol=length(tmp_var),
                     dimnames=list(NULL,
                                   names(tmp_var)))
  
  model$setVars(tmp_var)
  
  tmp<-c(as.numeric(p),pars$values_1[3:12])
  names(tmp)<-pars$name
  tmp <- matrix(rep(tmp,
                    each=nLayers), 
                nrow=nLayers, 
                ncol=length(tmp),
                dimnames=list(NULL,names(tmp)))
  tmp[1,"leftmost"]<-1
  model$setPars(tmp)
  t<-seq(from=0,to=300,length.out =100)
  out1<- model$dynamics(times=t)
  out_df<-(cbind(out1[,"time"],out1[,"C.165"]))
  colnames(out_df)=c("time","conc")
  return(out_df)
}


Fit0_1<-Model_fun1(parm_guess1)


Model_cost1<-function(p){
  
  out_df<-Model_fun1(p)
  Data<- as.data.frame(read_excel(path= "paper_data.xlsx",
                                  sheet="F1", 
                                  col_names=TRUE))
  cost_fn<-modCost(model=out_df,obs=Data) 
}


Fit1_1 <- modFit(f = Model_cost1, 
                 p = parm_guess1,
                 method = 'bobyqa',
                 lower = rep(0,2),
                 upper = rep(10,2)
)


MC_results1<- modMCMC(f=Model_cost1,
                      p=Fit1_1$par,
                      lower = rep(0,2),
                      upper = rep(Inf,2),
                      niter = 20000)


Fit2_1<-Model_fun1(Fit1_1$par)

write.csv(Fit2_1,file = "F:\\codes\\Experiments\\results\\fast\\noDOwithN_fast_conf_t2")
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
  tmp_par<-c(as.numeric(p1),pars$values_1[4:10])
  names(tmp_par)<-pars$name
  tmp_par <- matrix(rep(tmp_par,
                        each=nLayers), 
                    nrow=nLayers, 
                    ncol=length(tmp_par),
                    dimnames=list(NULL,names(tmp_par)))
  tmp_par[1,"leftmost"]<-1
  model$setPars(tmp_par)
  
  #define time and compile
  times <- seq(0,300,by=0.1)
  out<- model$dynamics(times=times)
  return(out)
}

Fit_soil<-Model_fun_soil(MC_results1$bestpar)
Fit_soil_end<-as.data.frame(Fit_soil[2280,401:801])
C_1<-mean(Fit_soil_end$`Fit_soil[2280, 401:801]`[1:60])
C_2<-mean(Fit_soil_end$`Fit_soil[2280, 401:801]`[60:120])
C_3<-mean(Fit_soil_end$`Fit_soil[2280, 401:801]`[120:180])
C_4<-mean(Fit_soil_end$`Fit_soil[2280, 401:801]`[180:240])
C_5<-mean(Fit_soil_end$`Fit_soil[2280, 401:801]`[240:300])

Data1<- as.data.frame(read_excel(path= "F:/codes/Experiments/column/tracer+substrate+bacteria_column/col1/paper_data.xlsx",
                                sheet="F1", 
                                col_names=TRUE))


SSR1<-modCost(obs=Data1,model = Fit2_1)

###################################################################################################################
en_sim<-unclass(as.POSIXlt(Sys.time()))

print(paste0("simulation time= ",
             (((en_sim$hour*60)+(en_sim$min)+(en_sim$sec/60))-
                ((st_sim$hour*60)+(st_sim$min)+(st_sim$sec/60))),
             "minutes")
)


plot(Data1,ylim=c(0,1),ylab="C/Co(-)",xlab="time (mins)",pch=16,col="red")
lines(Fit2_1,col="black",lwd="3")

# legend("topright", legend=c("observation", "model"),lty=c(0,1),lwd=c(0,3),pch=c(16,NA)) 


save.image(file="F:/codes/Experiments/results/noDowithN_fast_combdata_trial2.RData")
#sR <- sensRange(func = Model_fun, parms = MC_results$bestpar, parInput = MC_results$pars)

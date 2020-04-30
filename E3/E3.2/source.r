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
obs_data_out<- as.data.frame(read_excel(path= "paper_data.xlsx",
                                sheet="F1", 
                                col_names=TRUE))

obs_data_in<- as.data.frame(read_excel(path= "paper_data.xlsx",
                                sheet="F1_in", 
                                col_names=TRUE))

parm <- as.data.frame(read_excel("parameters.xlsx", 
                                 sheet = "pars",
                                 col_names = TRUE,
                                 range="A1:F3"))

parm_guess<-parm$values

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

Fit1<-Model_fun1(Fit1_1$par)

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

###################################################################################################################
en_sim<-unclass(as.POSIXlt(Sys.time()))

print(paste0("simulation time= ",
             (((en_sim$hour*60)+(en_sim$min)+(en_sim$sec/60))-
                ((st_sim$hour*60)+(st_sim$min)+(st_sim$sec/60))),
             "minutes")
)


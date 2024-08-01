######################################################################
######################################################################

##R Code for TerrorDecay-DiseaseControl Model with discrete time steps


#######################################################################
##Control starts once death toll reaches certain level
library(deSolve)
library(ggplot2)
library(scales)
library(pheatmap)
library(reshape)
library(here)
library(ggpubr)


###############################
Timesteps = 365 ##timestep after the first control

thres = 3  ###threshold value of death when control could start

CD = 30 ## control duration time for one control period

K_M = 3500 ###mosquito carrying capacity

#########step 1 create a vector for control on mosquito: C_M

C_M <- rep(0, Timesteps)


#####determine the first time when D_H > thres. 


S_H <- c()

I_H <- c()

R_H <- c()

D_H <- c()

S_M <- c()

I_M <- c()


S_H[1] <- 5000

I_H[1] <- 10

R_H[1] <- 0

D_H[1] <- 0

S_M[1] <- 1000

I_M[1] <- 0


b_H = 2.4657534e-05*2000/700 #human birth rate
beta_H = 0.00005 #human transmission rate
mu_H = 2.35616e-05*2000/700#natural mortality rate of human
r = 0.037 #recovery rate of human from infection 
delta_H = 3*5.468913e-05 ##rate how infected human to die
eta_M = 5 #egg laying rate - determine mosquito birth 
#K_M = 3500 #mosquito carrying capacity
beta_M = 0.0003#mosquito transmission rate
mu_M = 1/13#natural mortality of mosquito
t_int <- c() ###track all time time points when death > thres

for(t in 1: Timesteps)
{
  S_H[t+1] <- max(S_H[t] + b_H * (S_H[t] + I_H[t] + R_H[t]) - beta_H * I_M[t] * S_H[t] - mu_H *  S_H[t], 0)
  I_H[t+1] <- max(I_H[t] + beta_H * I_M[t] * S_H[t] - r * I_H[t] - mu_H * I_H[t] - delta_H * I_H[t], 0)
  R_H[t+1] <- max(R_H[t] + r * I_H[t] - mu_H * R_H[t], 0)
  D_H[t+1] <- max(D_H[t] + delta_H * I_H[t], 0)
  S_M[t+1] <- max(S_M[t] + eta_M * (S_M[t] + I_M[t]) * (1 -  eta_M * (S_M[t] + I_M[t])/ K_M) - beta_M *  I_H[t] * S_M[t] - mu_M * S_M[t] - C_M[t] * S_M[t], 0)
  I_M[t+1] <- max(I_M[t] + beta_M * I_H[t] * S_M[t] - mu_M * I_M[t] - C_M[t] * I_M[t], 0)
  
  ##first period release
  if(D_H[t] > thres)
  {
    t_int <- c(t_int, t)
  }
}

start_t <- min(t_int) #determine the first time point when death > thres and government starts to report its death toll


###update total time

Timesteps_new =  Timesteps + start_t  ###total time = time before control + time when control already starts

###case without control: serve as comparison group
sita = 0  ###theta: initial fear

tao = 0  ###tau: fear decay

S_H <- c() ##define one vector to save the time series output of susceptible humans

I_H <- c() ##define one vector to save the time series output of infected humans

R_H <- c() ##define one vector to save the time series output of recovered humans

D_H <- c() ##define one vector to save the time series output of death cases in human

S_M <- c() ##define one vector to save the time series output of susceptible mosquitoes

I_M <- c() ##define one vector to save the time series output of infected mosquitoes

C_M <- c() ##define one vector to save the time series of control efforts


S_H[1] <- 5000  ##initial value for the simulation

I_H[1] <- 10

R_H[1] <- 0

D_H[1] <- 0

S_M[1] <- 1000

I_M[1] <- 0

C_M[1] <- 0


b_H = 2.4657534e-05*2000/700 #human birth rate
beta_H = 0.00005 #human transmission rate
mu_H = 2.35616e-05*2000/700#natural mortality rate of human
r = 0.037 #recovery rate of human from infection 
delta_H = 3*5.468913e-05 ##rate how infected human to die
eta_M = 5 #egg laying rate - determine mosquito birth 
beta_M = 0.0003#mosquito transmission rate
mu_M = 1/13#natural mortality of mosquito


for(t in 1: Timesteps_new)
{
  FD <- function(x) {sita * exp(- x * tao)} ##fear of death equation
  S_H[t+1] <- max(S_H[t] + b_H * (S_H[t] + I_H[t] + R_H[t]) - beta_H * I_M[t] * S_H[t] - mu_H *  S_H[t], 0)
  I_H[t+1] <- max(I_H[t] + beta_H * I_M[t] * S_H[t] - r * I_H[t] - mu_H * I_H[t] - delta_H * I_H[t], 0)
  R_H[t+1] <- max(R_H[t] + r * I_H[t] - mu_H * R_H[t], 0)
  D_H[t+1] <- max(D_H[t] + delta_H * I_H[t], 0)
  S_M[t+1] <- max(S_M[t] + eta_M * (S_M[t] + I_M[t]) * (1 -  eta_M * (S_M[t] + I_M[t])/ K_M) - beta_M *  I_H[t] * S_M[t] - mu_M * S_M[t] - C_M[t] * S_M[t], 0)
  I_M[t+1] <- max(I_M[t] + beta_M * I_H[t] * S_M[t] - mu_M * I_M[t] - C_M[t] * I_M[t], 0)
  C_M[t+1] <- max(C_M[t] + FD(t) * D_H[t], 0)   
  
}

S_H_NoControl <- S_H

I_H_NoControl <- I_H

R_H_NoControl <- R_H

D_H_NoControl <- D_H

S_M_NoControl <- S_M

I_M_NoControl <- I_M

C_M_NoControl <- C_M

###############################################################################

###Fig. 5  disease dynamics under control periods = 1, 6, 365 with intermediate initial fear without decay and CD = 3


##############################################################################################

#################level 1 of theta

Timesteps = 365

C_int = 0.1 #

tau = 0 

period_number = 12  ##death toll released times per year

period_time = floor(365/period_number) ##first release when death toll number > thres

CD_period = 3

CD = round(min(CD_period, period_time), 0)

con_start <- c()

con_start[1] <- 1

for(kk in 1:period_number)
{
  con_start <- c(con_start, (start_t + (kk - 1) * round(period_time, 0)))
}

con_start <- c(con_start, (Timesteps_new + 1))

sita = C_int
tao = tau

C_M = rep(0, (Timesteps_new + 1))

S_H <- c()

I_H <- c()

R_H <- c()

D_H <- c()

S_M <- c()

I_M <- c()


S_H[1] <- 5000

I_H[1] <- 10

R_H[1] <- 0

D_H[1] <- 0

S_M[1] <- 1000

I_M[1] <- 0


b_H = 2.4657534e-05*2000/700 #human birth rate
beta_H = 0.00005 #human transmission rate
mu_H = 2.35616e-05*2000/700#natural mortality rate of human
r = 0.037 #recovery rate of human from infection 
delta_H = 3*5.468913e-05 ##rate how infected human to die
eta_M = 5 #egg laying rate - determine mosquito birth 
#K_M = 3500 #mosquito carrying capacity
beta_M = 0.0003#mosquito transmission rate
mu_M = 1/13#natural mortality of mosquito
t_int <- c()
#jk = 0

for(j in 1:(length(con_start)-1))
{
  
  jk <- 0
  
  
  for(t in con_start[j]:(con_start[j+1]-1))
  {
    
    FD <- function(x) {sita * exp(- x * tao)}
    S_H[t+1] <- max(S_H[t] + b_H * (S_H[t] + I_H[t] + R_H[t]) - beta_H * I_M[t] * S_H[t] - mu_H *  S_H[t], 0)
    I_H[t+1] <- max(I_H[t] + beta_H * I_M[t] * S_H[t] - r * I_H[t] - mu_H * I_H[t] - delta_H * I_H[t], 0)
    R_H[t+1] <- max(R_H[t] + r * I_H[t] - mu_H * R_H[t], 0)
    D_H[t+1] <- max(D_H[t] + delta_H * I_H[t], 0)
    S_M[t+1] <- max(S_M[t] + eta_M * (S_M[t] + I_M[t]) * (1 -  eta_M * (S_M[t] + I_M[t])/ K_M) - beta_M *  I_H[t] * S_M[t] - mu_M * S_M[t] - C_M[t] * S_M[t], 0)
    I_M[t+1] <- max(I_M[t] + beta_M * I_H[t] * S_M[t] - mu_M * I_M[t] - C_M[t] * I_M[t], 0)
    
    
    if(((t ==  con_start[j]) & (D_H[t] > thres)) | ((jk < CD) & (jk > 0)))
    {
      C_M[con_start[j]: (con_start[j] +  CD)] = FD(0:CD) * D_H[con_start[j]]/D_H[con_start[2]] ###every period, control always starts at time 0, lasting for CD days
      
      jk = jk + 1
      
      print(t)
      print(j)
    }
  }
}


C_M_period12_sitalevel1_nodecay = C_M

S_H_period12_sitalevel1_nodecay = S_H

I_H_period12_sitalevel1_nodecay = I_H

R_H_period12_sitalevel1_nodecay = R_H

D_H_period12_sitalevel1_nodecay = D_H

S_M_period12_sitalevel1_nodecay = S_M

I_M_period12_sitalevel1_nodecay = I_M

###calculate control efficacy

con_eff_period12_sitalevel1_nodecay <- (round(sum(D_H_NoControl)) - round(sum(D_H_period12_sitalevel1_nodecay)))/(sum(C_M_period12_sitalevel1_nodecay) + 0.0000000000001)


###########period 6

period_number = 14  ##death toll released times per year

period_time = floor(365/period_number) ##first release when death toll number > thres

CD_period = 3

CD = round(min(CD_period, period_time), 0)

con_start <- c()

con_start[1] <- 1

for(kk in 1:period_number)
{
  con_start <- c(con_start, (start_t + (kk - 1) * round(period_time, 0)))
}

con_start <- c(con_start, (Timesteps_new + 1))

sita = C_int
tao = tau

C_M = rep(0, (Timesteps_new + 1))

S_H <- c()

I_H <- c()

R_H <- c()

D_H <- c()

S_M <- c()

I_M <- c()


S_H[1] <- 5000

I_H[1] <- 10

R_H[1] <- 0

D_H[1] <- 0

S_M[1] <- 1000

I_M[1] <- 0


b_H = 2.4657534e-05*2000/700 #human birth rate
beta_H = 0.00005 #human transmission rate
mu_H = 2.35616e-05*2000/700#natural mortality rate of human
r = 0.037 #recovery rate of human from infection 
delta_H = 3*5.468913e-05 ##rate how infected human to die
eta_M = 5 #egg laying rate - determine mosquito birth 
#K_M = 3500 #mosquito carrying capacity
beta_M = 0.0003#mosquito transmission rate
mu_M = 1/13#natural mortality of mosquito
t_int <- c()
#jk = 0

for(j in 1:(length(con_start)-1))
{
  
  jk <- 0
  
  
  for(t in con_start[j]:(con_start[j+1]-1))
  {
    
    FD <- function(x) {sita * exp(- x * tao)}
    S_H[t+1] <- max(S_H[t] + b_H * (S_H[t] + I_H[t] + R_H[t]) - beta_H * I_M[t] * S_H[t] - mu_H *  S_H[t], 0)
    I_H[t+1] <- max(I_H[t] + beta_H * I_M[t] * S_H[t] - r * I_H[t] - mu_H * I_H[t] - delta_H * I_H[t], 0)
    R_H[t+1] <- max(R_H[t] + r * I_H[t] - mu_H * R_H[t], 0)
    D_H[t+1] <- max(D_H[t] + delta_H * I_H[t], 0)
    S_M[t+1] <- max(S_M[t] + eta_M * (S_M[t] + I_M[t]) * (1 -  eta_M * (S_M[t] + I_M[t])/ K_M) - beta_M *  I_H[t] * S_M[t] - mu_M * S_M[t] - C_M[t] * S_M[t], 0)
    I_M[t+1] <- max(I_M[t] + beta_M * I_H[t] * S_M[t] - mu_M * I_M[t] - C_M[t] * I_M[t], 0)
    
    
    if(((t ==  con_start[j]) & (D_H[t] > thres)) | ((jk < CD) & (jk > 0)))
    {
      C_M[con_start[j]: (con_start[j] +  CD)] = FD(0:CD) * D_H[con_start[j]]/D_H[con_start[2]] ###every period, control always starts at time 0, lasting for CD days
      
      jk = jk + 1
      
      print(t)
      print(j)
    }
  }
}

C_M_period14_sitalevel1_nodecay = C_M

S_H_period14_sitalevel1_nodecay = S_H

I_H_period14_sitalevel1_nodecay = I_H

R_H_period14_sitalevel1_nodecay = R_H

D_H_period14_sitalevel1_nodecay = D_H

S_M_period14_sitalevel1_nodecay = S_M

I_M_period14_sitalevel1_nodecay = I_M

###calculate control efficacy

con_eff_period14_sitalevel1_nodecay <- (round(sum(D_H_NoControl)) - round(sum(D_H_period14_sitalevel1_nodecay)))/(sum(C_M_period14_sitalevel1_nodecay) + 0.0000000000001)



###########################################################################################
#######period 365

period_number = 365  ##death toll released times per year

period_time = floor(365/period_number) ##first release when death toll number > thres

CD_period = 3

CD = round(min(CD_period, period_time), 0)

con_start <- c()

con_start[1] <- 1

for(kk in 1:period_number)
{
  con_start <- c(con_start, (start_t + (kk - 1) * round(period_time, 0)))
}

con_start <- c(con_start, (Timesteps_new + 1))

sita = C_int
tao = tau

C_M = rep(0, (Timesteps_new + 1))

S_H <- c()

I_H <- c()

R_H <- c()

D_H <- c()

S_M <- c()

I_M <- c()


S_H[1] <- 5000

I_H[1] <- 10

R_H[1] <- 0

D_H[1] <- 0

S_M[1] <- 1000

I_M[1] <- 0


b_H = 2.4657534e-05*2000/700 #human birth rate
beta_H = 0.00005 #human transmission rate
mu_H = 2.35616e-05*2000/700#natural mortality rate of human
r = 0.037 #recovery rate of human from infection 
delta_H = 3*5.468913e-05 ##rate how infected human to die
eta_M = 5 #egg laying rate - determine mosquito birth 
#K_M = 3500 #mosquito carrying capacity
beta_M = 0.0003#mosquito transmission rate
mu_M = 1/13#natural mortality of mosquito
t_int <- c()
#jk = 0

for(j in 1:(length(con_start)-1))
{
  
  jk <- 0
  
  
  for(t in con_start[j]:(con_start[j+1]-1))
  {
    
    FD <- function(x) {sita * exp(- x * tao)}
    S_H[t+1] <- max(S_H[t] + b_H * (S_H[t] + I_H[t] + R_H[t]) - beta_H * I_M[t] * S_H[t] - mu_H *  S_H[t], 0)
    I_H[t+1] <- max(I_H[t] + beta_H * I_M[t] * S_H[t] - r * I_H[t] - mu_H * I_H[t] - delta_H * I_H[t], 0)
    R_H[t+1] <- max(R_H[t] + r * I_H[t] - mu_H * R_H[t], 0)
    D_H[t+1] <- max(D_H[t] + delta_H * I_H[t], 0)
    S_M[t+1] <- max(S_M[t] + eta_M * (S_M[t] + I_M[t]) * (1 -  eta_M * (S_M[t] + I_M[t])/ K_M) - beta_M *  I_H[t] * S_M[t] - mu_M * S_M[t] - C_M[t] * S_M[t], 0)
    I_M[t+1] <- max(I_M[t] + beta_M * I_H[t] * S_M[t] - mu_M * I_M[t] - C_M[t] * I_M[t], 0)
    
    
    if(((t ==  con_start[j]) & (D_H[t] > thres)) | ((jk < CD) & (jk > 0)))
    {
      C_M[con_start[j]: (con_start[j] +  CD)] = FD(0:CD) * D_H[con_start[j]]/D_H[con_start[2]] ###every period, control always starts at time 0, lasting for CD days
      
      jk = jk + 1
      
      print(t)
      print(j)
    }
  }
}


C_M_period365_sitalevel1_nodecay = C_M

S_H_period365_sitalevel1_nodecay = S_H

I_H_period365_sitalevel1_nodecay = I_H

R_H_period365_sitalevel1_nodecay = R_H

D_H_period365_sitalevel1_nodecay = D_H

S_M_period365_sitalevel1_nodecay = S_M

I_M_period365_sitalevel1_nodecay = I_M

###calculate control efficacy

con_eff_period365_sitalevel1_nodecay <- (round(sum(D_H_NoControl)) - round(sum(D_H_period365_sitalevel1_nodecay)))/(sum(C_M_period365_sitalevel1_nodecay) + 0.0000000000001)

##########level 2 of theta

Timesteps = 365

C_int = 1 #

tau = 0 

period_number = 12  ##death toll released times per year

period_time = floor(365/period_number) ##first release when death toll number > thres

CD_period = 3

CD = round(min(CD_period, period_time), 0)

con_start <- c()

con_start[1] <- 1

for(kk in 1:period_number)
{
  con_start <- c(con_start, (start_t + (kk - 1) * round(period_time, 0)))
}

con_start <- c(con_start, (Timesteps_new + 1))

sita = C_int
tao = tau

C_M = rep(0, (Timesteps_new + 1))

S_H <- c()

I_H <- c()

R_H <- c()

D_H <- c()

S_M <- c()

I_M <- c()


S_H[1] <- 5000

I_H[1] <- 10

R_H[1] <- 0

D_H[1] <- 0

S_M[1] <- 1000

I_M[1] <- 0


b_H = 2.4657534e-05*2000/700 #human birth rate
beta_H = 0.00005 #human transmission rate
mu_H = 2.35616e-05*2000/700#natural mortality rate of human
r = 0.037 #recovery rate of human from infection 
delta_H = 3*5.468913e-05 ##rate how infected human to die
eta_M = 5 #egg laying rate - determine mosquito birth 
#K_M = 3500 #mosquito carrying capacity
beta_M = 0.0003#mosquito transmission rate
mu_M = 1/13#natural mortality of mosquito
t_int <- c()
#jk = 0

for(j in 1:(length(con_start)-1))
{
  
  jk <- 0
  
  
  for(t in con_start[j]:(con_start[j+1]-1))
  {
    
    FD <- function(x) {sita * exp(- x * tao)}
    S_H[t+1] <- max(S_H[t] + b_H * (S_H[t] + I_H[t] + R_H[t]) - beta_H * I_M[t] * S_H[t] - mu_H *  S_H[t], 0)
    I_H[t+1] <- max(I_H[t] + beta_H * I_M[t] * S_H[t] - r * I_H[t] - mu_H * I_H[t] - delta_H * I_H[t], 0)
    R_H[t+1] <- max(R_H[t] + r * I_H[t] - mu_H * R_H[t], 0)
    D_H[t+1] <- max(D_H[t] + delta_H * I_H[t], 0)
    S_M[t+1] <- max(S_M[t] + eta_M * (S_M[t] + I_M[t]) * (1 -  eta_M * (S_M[t] + I_M[t])/ K_M) - beta_M *  I_H[t] * S_M[t] - mu_M * S_M[t] - C_M[t] * S_M[t], 0)
    I_M[t+1] <- max(I_M[t] + beta_M * I_H[t] * S_M[t] - mu_M * I_M[t] - C_M[t] * I_M[t], 0)
    
    
    if(((t ==  con_start[j]) & (D_H[t] > thres)) | ((jk < CD) & (jk > 0)))
    {
      C_M[con_start[j]: (con_start[j] +  CD)] = FD(0:CD) * D_H[con_start[j]]/D_H[con_start[2]] ###every period, control always starts at time 0, lasting for CD days
      
      jk = jk + 1
      
      print(t)
      print(j)
    }
  }
}


C_M_period12_sitalevel2_nodecay = C_M

S_H_period12_sitalevel2_nodecay = S_H

I_H_period12_sitalevel2_nodecay = I_H

R_H_period12_sitalevel2_nodecay = R_H

D_H_period12_sitalevel2_nodecay = D_H

S_M_period12_sitalevel2_nodecay = S_M

I_M_period12_sitalevel2_nodecay = I_M

###calculate control efficacy

con_eff_period12_sitalevel2_nodecay <- (round(sum(D_H_NoControl)) - round(sum(D_H_period12_sitalevel2_nodecay)))/(sum(C_M_period12_sitalevel2_nodecay) + 0.0000000000001)


###########period 6

period_number = 14  ##death toll released times per year

period_time = floor(365/period_number) ##first release when death toll number > thres

CD_period = 3

CD = round(min(CD_period, period_time), 0)

con_start <- c()

con_start[1] <- 1

for(kk in 1:period_number)
{
  con_start <- c(con_start, (start_t + (kk - 1) * round(period_time, 0)))
}

con_start <- c(con_start, (Timesteps_new + 1))

sita = C_int
tao = tau

C_M = rep(0, (Timesteps_new + 1))

S_H <- c()

I_H <- c()

R_H <- c()

D_H <- c()

S_M <- c()

I_M <- c()


S_H[1] <- 5000

I_H[1] <- 10

R_H[1] <- 0

D_H[1] <- 0

S_M[1] <- 1000

I_M[1] <- 0


b_H = 2.4657534e-05*2000/700 #human birth rate
beta_H = 0.00005 #human transmission rate
mu_H = 2.35616e-05*2000/700#natural mortality rate of human
r = 0.037 #recovery rate of human from infection 
delta_H = 3*5.468913e-05 ##rate how infected human to die
eta_M = 5 #egg laying rate - determine mosquito birth 
#K_M = 3500 #mosquito carrying capacity
beta_M = 0.0003#mosquito transmission rate
mu_M = 1/13#natural mortality of mosquito
t_int <- c()
#jk = 0

for(j in 1:(length(con_start)-1))
{
  
  jk <- 0
  
  
  for(t in con_start[j]:(con_start[j+1]-1))
  {
    
    FD <- function(x) {sita * exp(- x * tao)}
    S_H[t+1] <- max(S_H[t] + b_H * (S_H[t] + I_H[t] + R_H[t]) - beta_H * I_M[t] * S_H[t] - mu_H *  S_H[t], 0)
    I_H[t+1] <- max(I_H[t] + beta_H * I_M[t] * S_H[t] - r * I_H[t] - mu_H * I_H[t] - delta_H * I_H[t], 0)
    R_H[t+1] <- max(R_H[t] + r * I_H[t] - mu_H * R_H[t], 0)
    D_H[t+1] <- max(D_H[t] + delta_H * I_H[t], 0)
    S_M[t+1] <- max(S_M[t] + eta_M * (S_M[t] + I_M[t]) * (1 -  eta_M * (S_M[t] + I_M[t])/ K_M) - beta_M *  I_H[t] * S_M[t] - mu_M * S_M[t] - C_M[t] * S_M[t], 0)
    I_M[t+1] <- max(I_M[t] + beta_M * I_H[t] * S_M[t] - mu_M * I_M[t] - C_M[t] * I_M[t], 0)
    
    
    if(((t ==  con_start[j]) & (D_H[t] > thres)) | ((jk < CD) & (jk > 0)))
    {
      C_M[con_start[j]: (con_start[j] +  CD)] = FD(0:CD) * D_H[con_start[j]]/D_H[con_start[2]] ###every period, control always starts at time 0, lasting for CD days
      
      jk = jk + 1
      
      print(t)
      print(j)
    }
  }
}

C_M_period14_sitalevel2_nodecay = C_M

S_H_period14_sitalevel2_nodecay = S_H

I_H_period14_sitalevel2_nodecay = I_H

R_H_period14_sitalevel2_nodecay = R_H

D_H_period14_sitalevel2_nodecay = D_H

S_M_period14_sitalevel2_nodecay = S_M

I_M_period14_sitalevel2_nodecay = I_M

###calculate control efficacy

con_eff_period14_sitalevel2_nodecay <- (round(sum(D_H_NoControl)) - round(sum(D_H_period14_sitalevel2_nodecay)))/(sum(C_M_period14_sitalevel2_nodecay) + 0.0000000000001)



###########################################################################################
#######period 365

period_number = 365  ##death toll released times per year

period_time = floor(365/period_number) ##first release when death toll number > thres

CD_period = 3

CD = round(min(CD_period, period_time), 0)

con_start <- c()

con_start[1] <- 1

for(kk in 1:period_number)
{
  con_start <- c(con_start, (start_t + (kk - 1) * round(period_time, 0)))
}

con_start <- c(con_start, (Timesteps_new + 1))

sita = C_int
tao = tau

C_M = rep(0, (Timesteps_new + 1))

S_H <- c()

I_H <- c()

R_H <- c()

D_H <- c()

S_M <- c()

I_M <- c()


S_H[1] <- 5000

I_H[1] <- 10

R_H[1] <- 0

D_H[1] <- 0

S_M[1] <- 1000

I_M[1] <- 0


b_H = 2.4657534e-05*2000/700 #human birth rate
beta_H = 0.00005 #human transmission rate
mu_H = 2.35616e-05*2000/700#natural mortality rate of human
r = 0.037 #recovery rate of human from infection 
delta_H = 3*5.468913e-05 ##rate how infected human to die
eta_M = 5 #egg laying rate - determine mosquito birth 
#K_M = 3500 #mosquito carrying capacity
beta_M = 0.0003#mosquito transmission rate
mu_M = 1/13#natural mortality of mosquito
t_int <- c()
#jk = 0

for(j in 1:(length(con_start)-1))
{
  
  jk <- 0
  
  
  for(t in con_start[j]:(con_start[j+1]-1))
  {
    
    FD <- function(x) {sita * exp(- x * tao)}
    S_H[t+1] <- max(S_H[t] + b_H * (S_H[t] + I_H[t] + R_H[t]) - beta_H * I_M[t] * S_H[t] - mu_H *  S_H[t], 0)
    I_H[t+1] <- max(I_H[t] + beta_H * I_M[t] * S_H[t] - r * I_H[t] - mu_H * I_H[t] - delta_H * I_H[t], 0)
    R_H[t+1] <- max(R_H[t] + r * I_H[t] - mu_H * R_H[t], 0)
    D_H[t+1] <- max(D_H[t] + delta_H * I_H[t], 0)
    S_M[t+1] <- max(S_M[t] + eta_M * (S_M[t] + I_M[t]) * (1 -  eta_M * (S_M[t] + I_M[t])/ K_M) - beta_M *  I_H[t] * S_M[t] - mu_M * S_M[t] - C_M[t] * S_M[t], 0)
    I_M[t+1] <- max(I_M[t] + beta_M * I_H[t] * S_M[t] - mu_M * I_M[t] - C_M[t] * I_M[t], 0)
    
    
    if(((t ==  con_start[j]) & (D_H[t] > thres)) | ((jk < CD) & (jk > 0)))
    {
      C_M[con_start[j]: (con_start[j] +  CD)] = FD(0:CD) * D_H[con_start[j]]/D_H[con_start[2]] ###every period, control always starts at time 0, lasting for CD days
      
      jk = jk + 1
      
      print(t)
      print(j)
    }
  }
}


C_M_period365_sitalevel2_nodecay = C_M

S_H_period365_sitalevel2_nodecay = S_H

I_H_period365_sitalevel2_nodecay = I_H

R_H_period365_sitalevel2_nodecay = R_H

D_H_period365_sitalevel2_nodecay = D_H

S_M_period365_sitalevel2_nodecay = S_M

I_M_period365_sitalevel2_nodecay = I_M

###calculate control efficacy

con_eff_period365_sitalevel2_nodecay <- (round(sum(D_H_NoControl)) - round(sum(D_H_period365_sitalevel2_nodecay)))/(sum(C_M_period365_sitalevel2_nodecay) + 0.0000000000001)

#########level 3 of theta

Timesteps = 365

C_int = 1.5 #

tau = 0 

period_number = 12  ##death toll released times per year

period_time = floor(365/period_number) ##first release when death toll number > thres

CD_period = 3

CD = round(min(CD_period, period_time), 0)

con_start <- c()

con_start[1] <- 1

for(kk in 1:period_number)
{
  con_start <- c(con_start, (start_t + (kk - 1) * round(period_time, 0)))
}

con_start <- c(con_start, (Timesteps_new + 1))

sita = C_int
tao = tau

C_M = rep(0, (Timesteps_new + 1))

S_H <- c()

I_H <- c()

R_H <- c()

D_H <- c()

S_M <- c()

I_M <- c()


S_H[1] <- 5000

I_H[1] <- 10

R_H[1] <- 0

D_H[1] <- 0

S_M[1] <- 1000

I_M[1] <- 0


b_H = 2.4657534e-05*2000/700 #human birth rate
beta_H = 0.00005 #human transmission rate
mu_H = 2.35616e-05*2000/700#natural mortality rate of human
r = 0.037 #recovery rate of human from infection 
delta_H = 3*5.468913e-05 ##rate how infected human to die
eta_M = 5 #egg laying rate - determine mosquito birth 
#K_M = 3500 #mosquito carrying capacity
beta_M = 0.0003#mosquito transmission rate
mu_M = 1/13#natural mortality of mosquito
t_int <- c()
#jk = 0

for(j in 1:(length(con_start)-1))
{
  
  jk <- 0
  
  
  for(t in con_start[j]:(con_start[j+1]-1))
  {
    
    FD <- function(x) {sita * exp(- x * tao)}
    S_H[t+1] <- max(S_H[t] + b_H * (S_H[t] + I_H[t] + R_H[t]) - beta_H * I_M[t] * S_H[t] - mu_H *  S_H[t], 0)
    I_H[t+1] <- max(I_H[t] + beta_H * I_M[t] * S_H[t] - r * I_H[t] - mu_H * I_H[t] - delta_H * I_H[t], 0)
    R_H[t+1] <- max(R_H[t] + r * I_H[t] - mu_H * R_H[t], 0)
    D_H[t+1] <- max(D_H[t] + delta_H * I_H[t], 0)
    S_M[t+1] <- max(S_M[t] + eta_M * (S_M[t] + I_M[t]) * (1 -  eta_M * (S_M[t] + I_M[t])/ K_M) - beta_M *  I_H[t] * S_M[t] - mu_M * S_M[t] - C_M[t] * S_M[t], 0)
    I_M[t+1] <- max(I_M[t] + beta_M * I_H[t] * S_M[t] - mu_M * I_M[t] - C_M[t] * I_M[t], 0)
    
    
    if(((t ==  con_start[j]) & (D_H[t] > thres)) | ((jk < CD) & (jk > 0)))
    {
      C_M[con_start[j]: (con_start[j] +  CD)] = FD(0:CD) * D_H[con_start[j]]/D_H[con_start[2]] ###every period, control always starts at time 0, lasting for CD days
      
      jk = jk + 1
      
      print(t)
      print(j)
    }
  }
}


C_M_period12_sitalevel3_nodecay = C_M

S_H_period12_sitalevel3_nodecay = S_H

I_H_period12_sitalevel3_nodecay = I_H

R_H_period12_sitalevel3_nodecay = R_H

D_H_period12_sitalevel3_nodecay = D_H

S_M_period12_sitalevel3_nodecay = S_M

I_M_period12_sitalevel3_nodecay = I_M

###calculate control efficacy

con_eff_period12_sitalevel3_nodecay <- (round(sum(D_H_NoControl)) - round(sum(D_H_period12_sitalevel3_nodecay)))/(sum(C_M_period12_sitalevel3_nodecay) + 0.0000000000001)


###########period 6

period_number = 14  ##death toll released times per year

period_time = floor(365/period_number) ##first release when death toll number > thres

CD_period = 3

CD = round(min(CD_period, period_time), 0)

con_start <- c()

con_start[1] <- 1

for(kk in 1:period_number)
{
  con_start <- c(con_start, (start_t + (kk - 1) * round(period_time, 0)))
}

con_start <- c(con_start, (Timesteps_new + 1))

sita = C_int
tao = tau

C_M = rep(0, (Timesteps_new + 1))

S_H <- c()

I_H <- c()

R_H <- c()

D_H <- c()

S_M <- c()

I_M <- c()


S_H[1] <- 5000

I_H[1] <- 10

R_H[1] <- 0

D_H[1] <- 0

S_M[1] <- 1000

I_M[1] <- 0


b_H = 2.4657534e-05*2000/700 #human birth rate
beta_H = 0.00005 #human transmission rate
mu_H = 2.35616e-05*2000/700#natural mortality rate of human
r = 0.037 #recovery rate of human from infection 
delta_H = 3*5.468913e-05 ##rate how infected human to die
eta_M = 5 #egg laying rate - determine mosquito birth 
#K_M = 3500 #mosquito carrying capacity
beta_M = 0.0003#mosquito transmission rate
mu_M = 1/13#natural mortality of mosquito
t_int <- c()
#jk = 0

for(j in 1:(length(con_start)-1))
{
  
  jk <- 0
  
  
  for(t in con_start[j]:(con_start[j+1]-1))
  {
    
    FD <- function(x) {sita * exp(- x * tao)}
    S_H[t+1] <- max(S_H[t] + b_H * (S_H[t] + I_H[t] + R_H[t]) - beta_H * I_M[t] * S_H[t] - mu_H *  S_H[t], 0)
    I_H[t+1] <- max(I_H[t] + beta_H * I_M[t] * S_H[t] - r * I_H[t] - mu_H * I_H[t] - delta_H * I_H[t], 0)
    R_H[t+1] <- max(R_H[t] + r * I_H[t] - mu_H * R_H[t], 0)
    D_H[t+1] <- max(D_H[t] + delta_H * I_H[t], 0)
    S_M[t+1] <- max(S_M[t] + eta_M * (S_M[t] + I_M[t]) * (1 -  eta_M * (S_M[t] + I_M[t])/ K_M) - beta_M *  I_H[t] * S_M[t] - mu_M * S_M[t] - C_M[t] * S_M[t], 0)
    I_M[t+1] <- max(I_M[t] + beta_M * I_H[t] * S_M[t] - mu_M * I_M[t] - C_M[t] * I_M[t], 0)
    
    
    if(((t ==  con_start[j]) & (D_H[t] > thres)) | ((jk < CD) & (jk > 0)))
    {
      C_M[con_start[j]: (con_start[j] +  CD)] = FD(0:CD) * D_H[con_start[j]]/D_H[con_start[2]] ###every period, control always starts at time 0, lasting for CD days
      
      jk = jk + 1
      
      print(t)
      print(j)
    }
  }
}

C_M_period14_sitalevel3_nodecay = C_M

S_H_period14_sitalevel3_nodecay = S_H

I_H_period14_sitalevel3_nodecay = I_H

R_H_period14_sitalevel3_nodecay = R_H

D_H_period14_sitalevel3_nodecay = D_H

S_M_period14_sitalevel3_nodecay = S_M

I_M_period14_sitalevel3_nodecay = I_M

###calculate control efficacy

con_eff_period14_sitalevel3_nodecay <- (round(sum(D_H_NoControl)) - round(sum(D_H_period14_sitalevel3_nodecay)))/(sum(C_M_period14_sitalevel3_nodecay) + 0.0000000000001)



###########################################################################################
#######period 365

period_number = 365  ##death toll released times per year

period_time = floor(365/period_number) ##first release when death toll number > thres

CD_period = 3

CD = round(min(CD_period, period_time), 0)

con_start <- c()

con_start[1] <- 1

for(kk in 1:period_number)
{
  con_start <- c(con_start, (start_t + (kk - 1) * round(period_time, 0)))
}

con_start <- c(con_start, (Timesteps_new + 1))

sita = C_int
tao = tau

C_M = rep(0, (Timesteps_new + 1))

S_H <- c()

I_H <- c()

R_H <- c()

D_H <- c()

S_M <- c()

I_M <- c()


S_H[1] <- 5000

I_H[1] <- 10

R_H[1] <- 0

D_H[1] <- 0

S_M[1] <- 1000

I_M[1] <- 0


b_H = 2.4657534e-05*2000/700 #human birth rate
beta_H = 0.00005 #human transmission rate
mu_H = 2.35616e-05*2000/700#natural mortality rate of human
r = 0.037 #recovery rate of human from infection 
delta_H = 3*5.468913e-05 ##rate how infected human to die
eta_M = 5 #egg laying rate - determine mosquito birth 
#K_M = 3500 #mosquito carrying capacity
beta_M = 0.0003#mosquito transmission rate
mu_M = 1/13#natural mortality of mosquito
t_int <- c()
#jk = 0

for(j in 1:(length(con_start)-1))
{
  
  jk <- 0
  
  
  for(t in con_start[j]:(con_start[j+1]-1))
  {
    
    FD <- function(x) {sita * exp(- x * tao)}
    S_H[t+1] <- max(S_H[t] + b_H * (S_H[t] + I_H[t] + R_H[t]) - beta_H * I_M[t] * S_H[t] - mu_H *  S_H[t], 0)
    I_H[t+1] <- max(I_H[t] + beta_H * I_M[t] * S_H[t] - r * I_H[t] - mu_H * I_H[t] - delta_H * I_H[t], 0)
    R_H[t+1] <- max(R_H[t] + r * I_H[t] - mu_H * R_H[t], 0)
    D_H[t+1] <- max(D_H[t] + delta_H * I_H[t], 0)
    S_M[t+1] <- max(S_M[t] + eta_M * (S_M[t] + I_M[t]) * (1 -  eta_M * (S_M[t] + I_M[t])/ K_M) - beta_M *  I_H[t] * S_M[t] - mu_M * S_M[t] - C_M[t] * S_M[t], 0)
    I_M[t+1] <- max(I_M[t] + beta_M * I_H[t] * S_M[t] - mu_M * I_M[t] - C_M[t] * I_M[t], 0)
    
    
    if(((t ==  con_start[j]) & (D_H[t] > thres)) | ((jk < CD) & (jk > 0)))
    {
      C_M[con_start[j]: (con_start[j] +  CD)] = FD(0:CD) * D_H[con_start[j]]/D_H[con_start[2]] ###every period, control always starts at time 0, lasting for CD days
      
      jk = jk + 1
      
      print(t)
      print(j)
    }
  }
}


C_M_period365_sitalevel3_nodecay = C_M

S_H_period365_sitalevel3_nodecay = S_H

I_H_period365_sitalevel3_nodecay = I_H

R_H_period365_sitalevel3_nodecay = R_H

D_H_period365_sitalevel3_nodecay = D_H

S_M_period365_sitalevel3_nodecay = S_M

I_M_period365_sitalevel3_nodecay = I_M

###calculate control efficacy

con_eff_period365_sitalevel3_nodecay <- (round(sum(D_H_NoControl)) - round(sum(D_H_period365_sitalevel3_nodecay)))/(sum(C_M_period365_sitalevel3_nodecay) + 0.0000000000001)


tiff(paste("Figure 5-MultiControl_nodecay_threesita_", "thres =",thres, "tau = ", tau, "CD =", CD_period, "K =", K_M, ".tiff"), width = 10, height = 10, units='in',res=600)

par(mfrow = c(1, 3))
par(mar = c(2,4,2,2))

plot(0:(Timesteps_new), C_M_period365_sitalevel2_nodecay, type = "l", ylim = c(0, 8.8),  ylab = "Control", lwd = 2, cex.lab = 1.3, col = "blue")
points(0:(Timesteps_new), C_M_period14_sitalevel2_nodecay, type = "l", col = "purple", lwd = 2)
points(0:(Timesteps_new), C_M_period12_sitalevel2_nodecay, type = "l", col = "red", lwd = 2)
legend("top", legend = expression(paste("Initial Fear", ~theta, ~"= 1")), cex = 1.3, bty = "n")


plot(0:(Timesteps_new), I_H_period365_sitalevel2_nodecay, type = "l", ylab = "Infected Human", cex.lab = 1.3, col = "blue", lwd = 2, ylim = c(0, 1500))
points(0:(Timesteps_new), I_H_period14_sitalevel2_nodecay, type = "l", col = "purple", lwd = 2)
points(0:(Timesteps_new), I_H_period12_sitalevel2_nodecay, type = "l", col = "red", lwd = 2)
legend("topright", legend =  c("log(Tot_Infec)", paste("=", round(log(sum(I_H_period12_sitalevel2_nodecay)), 3)), paste("=", round(log(sum(I_H_period14_sitalevel2_nodecay)), 3)), paste("=", round(log(sum(I_H_period365_sitalevel2_nodecay)), 3))), bty= "n", cex = 1.3, col = c("white","red", "purple", "blue"), lwd = c(2, 2, 2), lty = c(1, 1, 1, 1))

plot(0:(Timesteps_new), D_H_period14_sitalevel2_nodecay, type = "l", ylab = "Death Cases", cex.lab = 1.3, col = "purple", lwd = 2)
points(0:(Timesteps_new), D_H_period12_sitalevel2_nodecay, type = "l", col = "red", lwd = 2)
points(0:(Timesteps_new), D_H_period365_sitalevel2_nodecay, type = "l", col = "blue", lwd = 2)

legend("bottom", legend = c("Cont_Eff", paste(round(con_eff_period12_sitalevel2_nodecay, 2)), paste(round(con_eff_period14_sitalevel2_nodecay, 2)), paste(round(con_eff_period365_sitalevel2_nodecay, 2))), bty = "n", cex = 1.3, lwd = c(2, 2, 2, 2), lty = c(1, 1, 1, 1), col = c("white", "red", "purple", "blue"))

##################################################
dev.off()


tiff(paste("Figure S6-MultiControl_nodecay_threesita_susceptible_total_mos_human", "thres =",thres, "tau =", tau, "CD =", CD_period, "K =", K_M, ".tiff"), width = 10, height = 10, units='in',res=600)

par(mfrow = c(3, 3))
par(mar = c(2,4,2,2))

plot(0:(Timesteps_new), C_M_period365_sitalevel1_nodecay, type = "l", ylim = c(0, 8.8), ylab = "Control", lwd = 2, cex.lab = 1.3, col = "blue")
points(0:(Timesteps_new), C_M_period14_sitalevel1_nodecay, type = "l", col = "purple", lwd = 2)
points(0:(Timesteps_new), C_M_period12_sitalevel1_nodecay, type = "l", col = "red", lwd = 2)
legend("top", legend = expression(paste("Initial Fear", ~theta, ~"= 0.1")), cex = 1.3, bty = "n")
legend(100, 8, legend =  c("Cont_Period", " = 12", " = 14", " = 365"), bty = "n", cex = 1.25, col = c("white","red", "purple", "blue"), lwd = c(2, 2, 2), lty = c(1, 1, 1, 1))

plot(0:(Timesteps_new), S_H_period12_sitalevel1_nodecay, type = "l", xlab = "", ylim = c(0, 5000), cex.axis = 1, ylab = "Susceptible Human", cex.lab = 1.3, col = "red", lwd = 2)
points(0:(Timesteps_new), S_H_period14_sitalevel1_nodecay, type = "l", col = "purple", lwd = 2)
points(0:(Timesteps_new), S_H_period365_sitalevel1_nodecay, type = "l", col = "blue", lwd = 2)

plot(0:(Timesteps_new), S_M_period12_sitalevel1_nodecay + I_M_period12_sitalevel1_nodecay , type = "l", xlab = "", ylim = c(0, 1950), cex.axis = 1, ylab = "Mosquito Population", cex.lab = 1.3, col = "red", lwd = 2)
points(0:(Timesteps_new), S_M_period14_sitalevel1_nodecay + I_M_period12_sitalevel1_nodecay, type = "l", col = "purple", lwd = 2)
points(0:(Timesteps_new), S_M_period365_sitalevel1_nodecay + I_M_period365_sitalevel1_nodecay, type = "l", col = "blue", lwd = 2)
legend("topright", legend =  c("Mos_Pop", paste("=", round(sum(S_M_period12_sitalevel1_nodecay + I_M_period12_sitalevel1_nodecay), 2)), paste("=", round(sum(S_M_period14_sitalevel1_nodecay + I_M_period14_sitalevel1_nodecay), 2)), paste("=", round(sum(S_M_period365_sitalevel1_nodecay + I_M_period365_sitalevel1_nodecay), 2))), bty= "n", cex = 1.3, col = c("white","red", "purple", "blue"), lwd = c(2, 2, 2), lty = c(1, 1, 1, 1))

#######################################################################################
plot(0:(Timesteps_new), C_M_period365_sitalevel2_nodecay, type = "l", ylim = c(0, 8.8),  ylab = "Control", lwd = 2, cex.lab = 1.3, col = "blue")
points(0:(Timesteps_new), C_M_period14_sitalevel2_nodecay, type = "l", col = "purple", lwd = 2)
points(0:(Timesteps_new), C_M_period12_sitalevel2_nodecay, type = "l", col = "red", lwd = 2)
legend("top", legend = expression(paste("Initial Fear", ~theta, ~"= 1")), cex = 1.3, bty = "n")


plot(0:(Timesteps_new), S_H_period12_sitalevel2_nodecay, type = "l", xlab = "", cex.lab = 1.3 , cex.axis = 1, ylim = c(0, 5000), ylab = "Susceptible Human", col = "red", lwd = 2)
points(0:(Timesteps_new), S_H_period14_sitalevel2_nodecay, type = "l", col = "purple", lwd = 2)
points(0:(Timesteps_new), S_H_period365_sitalevel2_nodecay, type = "l", col = "blue", lwd = 2)

plot(0:(Timesteps_new), S_M_period12_sitalevel2_nodecay + I_M_period1_sitalevel2_nodecay , type = "l", xlab = "", ylim = c(0, 1950), cex.axis = 1, ylab = "Mosquito Population", cex.lab = 1.3, col = "red", lwd = 2)
points(0:(Timesteps_new), S_M_period14_sitalevel2_nodecay + I_M_period12_sitalevel2_nodecay, type = "l", col = "purple", lwd = 2)
points(0:(Timesteps_new), S_M_period365_sitalevel2_nodecay + I_M_period365_sitalevel2_nodecay, type = "l", col = "blue", lwd = 2)
legend("topright", legend =  c("Mos_Pop", paste("=", round(sum(S_M_period12_sitalevel2_nodecay + I_M_period12_sitalevel2_nodecay), 2)), paste("=", round(sum(S_M_period14_sitalevel2_nodecay + I_M_period14_sitalevel2_nodecay), 2)), paste("=", round(sum(S_M_period365_sitalevel2_nodecay + I_M_period365_sitalevel2_nodecay), 2))), bty= "n", cex = 1.3, col = c("white","red", "purple", "blue"), lwd = c(2, 2, 2), lty = c(1, 1, 1, 1))


#################################################################################
plot(0:(Timesteps_new), C_M_period365_sitalevel3_nodecay, type = "l", ylim = c(0, 8.8), ylab = "Control", lwd = 2, cex.lab = 1.3, col = "blue")
points(0:(Timesteps_new), C_M_period14_sitalevel3_nodecay, type = "l", col = "purple", lwd = 2)
points(0:(Timesteps_new), C_M_period12_sitalevel3_nodecay, type = "l", col = "red", lwd = 2)
legend("top", legend = expression(paste("Initial Fear", ~theta, ~"= 1.5")), cex = 1.3, bty = "n")


plot(0:(Timesteps_new), S_H_period12_sitalevel3_nodecay, type = "l", xlab = "", ylim = c(0, 5000), cex.axis = 1, ylab = "Susceptible Human", cex.lab = 1.3, col = "red", lwd = 2)
points(0:(Timesteps_new), S_H_period14_sitalevel3_nodecay, type = "l", col = "purple", lwd = 2)
points(0:(Timesteps_new), S_H_period365_sitalevel3_nodecay, type = "l", col = "blue", lwd = 2)

plot(0:(Timesteps_new), S_M_period12_sitalevel3_nodecay + I_M_period12_sitalevel3_nodecay , type = "l", xlab = "", ylim = c(0, 1950), cex.axis = 1, ylab = "Mosquito Population", cex.lab = 1.3, col = "red", lwd = 2)
points(0:(Timesteps_new), S_M_period14_sitalevel3_nodecay + I_M_period14_sitalevel3_nodecay, type = "l", col = "purple", lwd = 2)
points(0:(Timesteps_new), S_M_period365_sitalevel3_nodecay + I_M_period365_sitalevel3_nodecay, type = "l", col = "blue", lwd = 2)
legend("topright", legend =  c("Mos_Pop", paste("=", round(sum(S_M_period12_sitalevel3_nodecay + I_M_period12_sitalevel3_nodecay), 2)), paste("=", round(sum(S_M_period14_sitalevel3_nodecay + I_M_period14_sitalevel3_nodecay), 2)), paste("=", round(sum(S_M_period365_sitalevel3_nodecay + I_M_period365_sitalevel3_nodecay), 2))), bty= "n", cex = 1.3, col = c("white","red", "purple", "blue"), lwd = c(2, 2, 2), lty = c(1, 1, 1, 1))

dev.off()



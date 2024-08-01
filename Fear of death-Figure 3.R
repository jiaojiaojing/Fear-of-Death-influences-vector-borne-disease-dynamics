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

####Fig. 3 combined effect from both theta and tau


###############################

C_int_level = seq(0.001, 2.5, 0.05) ##set up levels of theta

tau_level = seq(0, 0.5, 0.25) ##set up levels of tau

period_number =  1  ##death toll released times per year-- one control period

period_time = floor(365/period_number) ###calculate the time interval among two controls, since only one control, this value would approach one year

con_start <- c() ###track all start time points for control efforts

CD_period = 30

CD = round(min(CD_period, (period_time - 1)), 0)  ###choose the CD as the minimum value of the setup duration and time interval between controls

for(kk in 1:period_number)
{
  con_start <- c(con_start, (start_t + (kk - 1) * round(period_time, 0)))
}


tot_infec_sitalevel1_nodecay <- matrix(NA, length(tau_level), length(C_int_level))

max_infec_sitalevel1_nodecay <- matrix(NA, length(tau_level), length(C_int_level))

tot_death_sitalevel1_nodecay <- matrix(NA, length(tau_level), length(C_int_level))

cont_eff_sitalevel1_nodecay <- matrix(NA, length(tau_level), length(C_int_level))

#########st

for(ii in 1:length(tau_level))
{
  for(jj in 1: length(C_int_level))
  {
    sita = C_int_level[jj]
    
    tao = tau_level[ii]
    
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
    beta_M = 0.0003#mosquito transmission rate
    mu_M = 1/13#natural mortality of mosquito
    t_int <- c()
    
    for(j in 1: length(con_start)) 
    {
      jk <- 0
      
      for(t in 1:Timesteps_new)
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
          C_M[con_start[j]: (con_start[j] +  CD)] = FD(0:CD) * D_H[t]/D_H[con_start[1]] ###every period, control always starts at time 0, lasting for CD days
          
          jk = jk + 1
          
          print(t)
          print(j)
        }
      }
    }
    
    tot_infec_sitalevel1_nodecay[ii, jj] <- sum(I_H)
    
    max_infec_sitalevel1_nodecay[ii, jj] <- max(I_H)
    
    tot_death_sitalevel1_nodecay[ii, jj] <- sum(D_H)
    
    cont_eff_sitalevel1_nodecay[ii,jj] <- (sum(D_H_NoControl) - sum(D_H))/(sum(C_M) + 0.0000000000001)
    
    
  }
}

x_axis <- C_int_level

tiff(paste("Figure 3-sita_tau-OneControl_thres =",thres,"CD =", CD, "K =", K_M, ".tiff"), width = 10, height = 10, units='in',res=600)

par(mfrow = c(2, 2))
par(mar = c(2,4.5,2,2))

plot(x_axis, log(tot_infec_sitalevel1_nodecay[1,]), type = "l", lwd = 2, ylab = "Total Infected Human (log)", xlab = "", cex.lab = 2, cex.axis = 2, col = "red")
points(x_axis, log(tot_infec_sitalevel1_nodecay[2,]), type = "l", lwd = 2, col = "purple")
points(x_axis, log(tot_infec_sitalevel1_nodecay[3,]), type = "l", lwd = 2, col = "blue")


plot(x_axis, log(max_infec_sitalevel1_nodecay[1,]), type = "l", lwd = 2, ylab = "Maximum Infected Human (log)", xlab = "", cex.lab = 2, cex.axis = 2, col = "red")
points(x_axis, log(max_infec_sitalevel1_nodecay[2,]), type = "l", lwd = 2, col = "purple")
points(x_axis, log(max_infec_sitalevel1_nodecay[3,]), type = "l", lwd = 2, col = "blue")
legend("topright", legend = c(expression(paste("Terror Decay (", tau, ")")), paste("=", tau_level[1]), paste("=",tau_level[2]), paste("=", tau_level[3])), lty = c(1, 1 , 1, 1), lwd = c(2, 2, 2, 2), col = c("white","red", "purple", "blue"), bty = "n", cex = 2)


plot(x_axis, log(tot_death_sitalevel1_nodecay[1,]), type = "l", lwd = 2, ylab = "Total Deaths (log)", xlab = "", cex.lab = 2, cex.axis = 2, col = "red")
points(x_axis, log(tot_death_sitalevel1_nodecay[2,]), type = "l", lwd = 2, col = "purple")
points(x_axis, log(tot_death_sitalevel1_nodecay[3,]), type = "l", lwd = 2, col = "blue")

plot(x_axis, log(cont_eff_sitalevel1_nodecay[1,]), type = "l", lwd = 2, ylab = "Control Efficacy (log)", xlab = "", cex.lab = 2, cex.axis = 2, col = "red", ylim = c(min(log(cont_eff_sitalevel1_nodecay)), max(log(cont_eff_sitalevel1_nodecay))))
points(x_axis, log(cont_eff_sitalevel1_nodecay[2,]), type = "l", lwd = 2, col = "purple")
points(x_axis, log(cont_eff_sitalevel1_nodecay[3,]), type = "l", lwd = 2, col = "blue")

dev.off()



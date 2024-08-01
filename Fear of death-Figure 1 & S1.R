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

################################Figure 1#######################################

###############################################################################

####case 1 for theta value

C_int <- 0.0  ### the value later theta would be assigned to

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


C_M <- 0

A1 <- c()  ##vector to keep control efforts

A1[1] <- C_M ##the first value of control effort time series

b_H = 2.4657534e-05*2000/700 #human birth rate
beta_H = 0.00005 #human transmission rate
mu_H = 2.35616e-05*2000/700#natural mortality rate of human
r = 0.037 #recovery rate of human from infection 
delta_H = 3*5.468913e-05 ##rate how infected human to die
eta_M = 5 #egg laying rate - determine mosquito birth 
beta_M = 0.0003#mosquito transmission rate
mu_M = 1/13#natural mortality of mosquito

jk1 = 0  ##tractor for CD conditions

for(t in 1: Timesteps_new)
{
  S_H[t+1] <- max(S_H[t] + b_H * (S_H[t] + I_H[t] + R_H[t]) - beta_H * I_M[t] * S_H[t] - mu_H *  S_H[t], 0)
  I_H[t+1] <- max(I_H[t] + beta_H * I_M[t] * S_H[t] - r * I_H[t] - mu_H * I_H[t] - delta_H * I_H[t], 0)
  R_H[t+1] <- max(R_H[t] + r * I_H[t] - mu_H * R_H[t], 0)
  D_H[t+1] <- max(D_H[t] + delta_H * I_H[t], 0)
  S_M[t+1] <- max(S_M[t] + eta_M * (S_M[t] + I_M[t]) * (1 -  eta_M * (S_M[t] + I_M[t])/ K_M) - beta_M *  I_H[t] * S_M[t] - mu_M * S_M[t] - C_M * S_M[t], 0)
  I_M[t+1] <- max(I_M[t] + beta_M * I_H[t] * S_M[t] - mu_M * I_M[t] - C_M * I_M[t], 0)
  
  if((D_H[t] > thres) & (jk1 < CD))
  {
    C_M = C_int
    
    jk1 = jk1 + 1
  }
  
  else 
  {
    C_M = 0
  }
  
  A1 <- c(A1, C_M)
}

C_M_sita_case1 <- A1

S_H_sita_case1 <- S_H

I_H_sita_case1 <- I_H

R_H_sita_case1 <- R_H

D_H_sita_case1 <- D_H

S_M_sita_case1 <- S_M

I_M_sita_case1 <- I_M


con_eff_sita_case1 <- (round(sum(D_H_NoControl)) - round(sum(D_H_sita_case1)))/sum(C_M_sita_case1)

#####case 2 for theta change

C_int <- 0.5

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


C_M <- 0

A1 <- c()

A1[1] <- C_M ##keep track of all control over the whole year

b_H = 2.4657534e-05*2000/700 #human birth rate
beta_H = 0.00005 #human transmission rate
mu_H = 2.35616e-05*2000/700#natural mortality rate of human
r = 0.037 #recovery rate of human from infection 
delta_H = 3*5.468913e-05 ##rate how infected human to die
eta_M = 5 #egg laying rate - determine mosquito birth 
#K_M = 3500 #mosquito carrying capacity
beta_M = 0.0003#mosquito transmission rate
mu_M = 1/13#natural mortality of mosquito

jk1 = 0

for(t in 1: Timesteps_new)
{
  S_H[t+1] <- max(S_H[t] + b_H * (S_H[t] + I_H[t] + R_H[t]) - beta_H * I_M[t] * S_H[t] - mu_H *  S_H[t], 0)
  I_H[t+1] <- max(I_H[t] + beta_H * I_M[t] * S_H[t] - r * I_H[t] - mu_H * I_H[t] - delta_H * I_H[t], 0)
  R_H[t+1] <- max(R_H[t] + r * I_H[t] - mu_H * R_H[t], 0)
  D_H[t+1] <- max(D_H[t] + delta_H * I_H[t], 0)
  S_M[t+1] <- max(S_M[t] + eta_M * (S_M[t] + I_M[t]) * (1 -  eta_M * (S_M[t] + I_M[t])/ K_M) - beta_M *  I_H[t] * S_M[t] - mu_M * S_M[t] - C_M * S_M[t], 0)
  I_M[t+1] <- max(I_M[t] + beta_M * I_H[t] * S_M[t] - mu_M * I_M[t] - C_M * I_M[t], 0)
  
  if((D_H[t] > thres) & (jk1 < CD))
  {
    C_M = C_int
    
    jk1 = jk1 + 1
  }
  
  else 
  {
    C_M = 0
  }
  
  A1 <- c(A1, C_M)
}

C_M_sita_case2 <- A1

S_H_sita_case2 <- S_H

I_H_sita_case2 <- I_H

R_H_sita_case2 <- R_H

D_H_sita_case2 <- D_H

S_M_sita_case2 <- S_M

I_M_sita_case2 <- I_M


con_eff_sita_case2 <- (round(sum(D_H_NoControl)) - round(sum(D_H_sita_case2)))/sum(C_M_sita_case2)


######case 3 for theta change

C_int <- 1.5

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


C_M <- 0

A1 <- c()

A1[1] <- C_M ##keep track of all control over the whole year

b_H = 2.4657534e-05*2000/700 #human birth rate
beta_H = 0.00005 #human transmission rate
mu_H = 2.35616e-05*2000/700#natural mortality rate of human
r = 0.037 #recovery rate of human from infection 
delta_H = 3*5.468913e-05 ##rate how infected human to die
eta_M = 5 #egg laying rate - determine mosquito birth 
#K_M = 3500 #mosquito carrying capacity
beta_M = 0.0003#mosquito transmission rate
mu_M = 1/13#natural mortality of mosquito

jk1 = 0

for(t in 1: Timesteps_new)
{
  S_H[t+1] <- max(S_H[t] + b_H * (S_H[t] + I_H[t] + R_H[t]) - beta_H * I_M[t] * S_H[t] - mu_H *  S_H[t], 0)
  I_H[t+1] <- max(I_H[t] + beta_H * I_M[t] * S_H[t] - r * I_H[t] - mu_H * I_H[t] - delta_H * I_H[t], 0)
  R_H[t+1] <- max(R_H[t] + r * I_H[t] - mu_H * R_H[t], 0)
  D_H[t+1] <- max(D_H[t] + delta_H * I_H[t], 0)
  S_M[t+1] <- max(S_M[t] + eta_M * (S_M[t] + I_M[t]) * (1 -  eta_M * (S_M[t] + I_M[t])/ K_M) - beta_M *  I_H[t] * S_M[t] - mu_M * S_M[t] - C_M * S_M[t], 0)
  I_M[t+1] <- max(I_M[t] + beta_M * I_H[t] * S_M[t] - mu_M * I_M[t] - C_M * I_M[t], 0)
  
  if((D_H[t] > thres) & (jk1 < CD))
  {
    C_M = C_int
    
    jk1 = jk1 + 1
  }
  
  else 
  {
    C_M = 0
  }
  
  A1 <- c(A1, C_M)
}

C_M_sita_case3 <- A1

S_H_sita_case3 <- S_H

I_H_sita_case3 <- I_H

R_H_sita_case3 <- R_H

D_H_sita_case3 <- D_H

S_M_sita_case3 <- S_M

I_M_sita_case3 <- I_M


con_eff_sita_case3 <- (round(sum(D_H_NoControl)) - round(sum(D_H_sita_case3)))/sum(C_M_sita_case3)


tiff(paste("Figure 1-Sita-OneControl_NoDecay_thres =",thres, "CD =", CD, "K =", K_M, ".tiff"), width = 10, height = 10, units='in',res=600)

par(mfrow = c(2, 2))
par(mar = c(2,4.5,2,2))

plot(0:(Timesteps_new), C_M_sita_case3, type = "l", lwd = 2, ylab = "Control", col = "blue", lty = 1,  ylim = c(min(C_M_sita_case1, C_M_sita_case2, C_M_sita_case3), max(C_M_sita_case1, C_M_sita_case2, C_M_sita_case3)), cex.lab = 2)
points(0:(Timesteps_new), C_M_sita_case2, type = "l", lwd = 2, col = "purple", cex.lab = 2)
points(0:(Timesteps_new), C_M_sita_case1, type = "l", lty = 5, lwd = 2, col = "red", cex.lab = 2)

plot(0:(Timesteps_new), S_H_sita_case1, type = "l", lwd = 2, ylab = "Susceptible Human", col = "red",lty = 5,  ylim = c(min(S_H_sita_case1, S_H_sita_case2, S_H_sita_case3), max(S_H_sita_case1, S_H_sita_case2, S_H_sita_case3)), cex.lab = 2)
points(0:(Timesteps_new), S_H_sita_case2, type = "l", lwd = 2, col = "purple", cex.lab = 2)
points(0:(Timesteps_new), S_H_sita_case3, type = "l", lwd = 2, col = "blue", cex.lab = 2)
legend("top", legend = c("Initial Fear", expression(paste(~theta, ~"= 0")), expression(paste(~theta, ~"= 1")), expression(paste(~theta, ~"= 1.5"))), lty = c(1, 5, 1 , 1), lwd = c(2, 2, 2, 2), col = c("white", "red", "purple", "blue"), bty = "n", cex = 2)

plot(0:(Timesteps_new), I_H_sita_case1, type = "l", lwd = 2, col = "red",lty = 5, ylim = c(min(I_H_sita_case1, I_H_sita_case2, I_H_sita_case3) - 80, max(I_H_sita_case1, I_H_sita_case2, I_H_sita_case3)), ylab = "Infected Human", cex.lab = 2)
points(0:(Timesteps_new), I_H_sita_case2, type = "l", lwd = 2, col = "purple", cex.lab = 2)
points(0:(Timesteps_new), I_H_sita_case3, type = "l", lwd = 2, col = "blue", cex.lab = 2)
legend("topright", legend = c("log(Tot_Infec)", paste("=",round(log(sum(I_H_sita_case1)), 2)), paste("=", round(log(sum(I_H_sita_case2)), 5)), paste("=", round(log(sum(I_H_sita_case3)), 2))), lty = c(1, 5, 1 , 1), lwd = c(2, 2, 2, 2), col = c("white", "red", "purple", "blue"), bty = "n", cex = 2)

plot(0:(Timesteps_new), D_H_sita_case1, type = "l", lwd = 2, ylab = "Death Cases", col = "red", lty = 5, ylim = c(min(D_H_sita_case1, D_H_sita_case2, D_H_sita_case3),max(D_H_sita_case1, D_H_sita_case2, D_H_sita_case3)), cex.lab = 2)
points(0:(Timesteps_new), D_H_sita_case2, type = "l", lwd = 2, col = "purple3", cex.lab = 2)
points(0:(Timesteps_new), D_H_sita_case3, type = "l", lwd = 2, col = "blue", cex.lab = 2)

dev.off()

#################################################CD change

###############################
Timesteps = 365 ##timestep after the first control

thres = 3  ###threshold value of death when control could start

CD = 3 ## control duration time for one control period

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

################################Figure 1#######################################

###############################################################################

####case 1 for theta value

C_int <- 0.0  ### the value later theta would be assigned to

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


C_M <- 0

A1 <- c()  ##vector to keep control efforts

A1[1] <- C_M ##the first value of control effort time series

b_H = 2.4657534e-05*2000/700 #human birth rate
beta_H = 0.00005 #human transmission rate
mu_H = 2.35616e-05*2000/700#natural mortality rate of human
r = 0.037 #recovery rate of human from infection 
delta_H = 3*5.468913e-05 ##rate how infected human to die
eta_M = 5 #egg laying rate - determine mosquito birth 
beta_M = 0.0003#mosquito transmission rate
mu_M = 1/13#natural mortality of mosquito

jk1 = 0  ##tractor for CD conditions

for(t in 1: Timesteps_new)
{
  S_H[t+1] <- max(S_H[t] + b_H * (S_H[t] + I_H[t] + R_H[t]) - beta_H * I_M[t] * S_H[t] - mu_H *  S_H[t], 0)
  I_H[t+1] <- max(I_H[t] + beta_H * I_M[t] * S_H[t] - r * I_H[t] - mu_H * I_H[t] - delta_H * I_H[t], 0)
  R_H[t+1] <- max(R_H[t] + r * I_H[t] - mu_H * R_H[t], 0)
  D_H[t+1] <- max(D_H[t] + delta_H * I_H[t], 0)
  S_M[t+1] <- max(S_M[t] + eta_M * (S_M[t] + I_M[t]) * (1 -  eta_M * (S_M[t] + I_M[t])/ K_M) - beta_M *  I_H[t] * S_M[t] - mu_M * S_M[t] - C_M * S_M[t], 0)
  I_M[t+1] <- max(I_M[t] + beta_M * I_H[t] * S_M[t] - mu_M * I_M[t] - C_M * I_M[t], 0)
  
  if((D_H[t] > thres) & (jk1 < CD))
  {
    C_M = C_int
    
    jk1 = jk1 + 1
  }
  
  else 
  {
    C_M = 0
  }
  
  A1 <- c(A1, C_M)
}

C_M_sita_case1 <- A1

S_H_sita_case1 <- S_H

I_H_sita_case1 <- I_H

R_H_sita_case1 <- R_H

D_H_sita_case1 <- D_H

S_M_sita_case1 <- S_M

I_M_sita_case1 <- I_M


con_eff_sita_case1 <- (round(sum(D_H_NoControl)) - round(sum(D_H_sita_case1)))/sum(C_M_sita_case1)

#####case 2 for theta change

C_int <- 0.5

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


C_M <- 0

A1 <- c()

A1[1] <- C_M ##keep track of all control over the whole year

b_H = 2.4657534e-05*2000/700 #human birth rate
beta_H = 0.00005 #human transmission rate
mu_H = 2.35616e-05*2000/700#natural mortality rate of human
r = 0.037 #recovery rate of human from infection 
delta_H = 3*5.468913e-05 ##rate how infected human to die
eta_M = 5 #egg laying rate - determine mosquito birth 
#K_M = 3500 #mosquito carrying capacity
beta_M = 0.0003#mosquito transmission rate
mu_M = 1/13#natural mortality of mosquito

jk1 = 0

for(t in 1: Timesteps_new)
{
  S_H[t+1] <- max(S_H[t] + b_H * (S_H[t] + I_H[t] + R_H[t]) - beta_H * I_M[t] * S_H[t] - mu_H *  S_H[t], 0)
  I_H[t+1] <- max(I_H[t] + beta_H * I_M[t] * S_H[t] - r * I_H[t] - mu_H * I_H[t] - delta_H * I_H[t], 0)
  R_H[t+1] <- max(R_H[t] + r * I_H[t] - mu_H * R_H[t], 0)
  D_H[t+1] <- max(D_H[t] + delta_H * I_H[t], 0)
  S_M[t+1] <- max(S_M[t] + eta_M * (S_M[t] + I_M[t]) * (1 -  eta_M * (S_M[t] + I_M[t])/ K_M) - beta_M *  I_H[t] * S_M[t] - mu_M * S_M[t] - C_M * S_M[t], 0)
  I_M[t+1] <- max(I_M[t] + beta_M * I_H[t] * S_M[t] - mu_M * I_M[t] - C_M * I_M[t], 0)
  
  if((D_H[t] > thres) & (jk1 < CD))
  {
    C_M = C_int
    
    jk1 = jk1 + 1
  }
  
  else 
  {
    C_M = 0
  }
  
  A1 <- c(A1, C_M)
}

C_M_sita_case2 <- A1

S_H_sita_case2 <- S_H

I_H_sita_case2 <- I_H

R_H_sita_case2 <- R_H

D_H_sita_case2 <- D_H

S_M_sita_case2 <- S_M

I_M_sita_case2 <- I_M


con_eff_sita_case2 <- (round(sum(D_H_NoControl)) - round(sum(D_H_sita_case2)))/sum(C_M_sita_case2)


######case 3 for theta change

C_int <- 1.5

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


C_M <- 0

A1 <- c()

A1[1] <- C_M ##keep track of all control over the whole year

b_H = 2.4657534e-05*2000/700 #human birth rate
beta_H = 0.00005 #human transmission rate
mu_H = 2.35616e-05*2000/700#natural mortality rate of human
r = 0.037 #recovery rate of human from infection 
delta_H = 3*5.468913e-05 ##rate how infected human to die
eta_M = 5 #egg laying rate - determine mosquito birth 
#K_M = 3500 #mosquito carrying capacity
beta_M = 0.0003#mosquito transmission rate
mu_M = 1/13#natural mortality of mosquito

jk1 = 0

for(t in 1: Timesteps_new)
{
  S_H[t+1] <- max(S_H[t] + b_H * (S_H[t] + I_H[t] + R_H[t]) - beta_H * I_M[t] * S_H[t] - mu_H *  S_H[t], 0)
  I_H[t+1] <- max(I_H[t] + beta_H * I_M[t] * S_H[t] - r * I_H[t] - mu_H * I_H[t] - delta_H * I_H[t], 0)
  R_H[t+1] <- max(R_H[t] + r * I_H[t] - mu_H * R_H[t], 0)
  D_H[t+1] <- max(D_H[t] + delta_H * I_H[t], 0)
  S_M[t+1] <- max(S_M[t] + eta_M * (S_M[t] + I_M[t]) * (1 -  eta_M * (S_M[t] + I_M[t])/ K_M) - beta_M *  I_H[t] * S_M[t] - mu_M * S_M[t] - C_M * S_M[t], 0)
  I_M[t+1] <- max(I_M[t] + beta_M * I_H[t] * S_M[t] - mu_M * I_M[t] - C_M * I_M[t], 0)
  
  if((D_H[t] > thres) & (jk1 < CD))
  {
    C_M = C_int
    
    jk1 = jk1 + 1
  }
  
  else 
  {
    C_M = 0
  }
  
  A1 <- c(A1, C_M)
}

C_M_sita_case3 <- A1

S_H_sita_case3 <- S_H

I_H_sita_case3 <- I_H

R_H_sita_case3 <- R_H

D_H_sita_case3 <- D_H

S_M_sita_case3 <- S_M

I_M_sita_case3 <- I_M


con_eff_sita_case3 <- (round(sum(D_H_NoControl)) - round(sum(D_H_sita_case3)))/sum(C_M_sita_case3)


tiff(paste("Figure S1-Sita-OneControl_NoDecay_thres =",thres, "CD =", CD, "K =", K_M, ".tiff"), width = 10, height = 10, units='in',res=600)

par(mfrow = c(2, 2))
par(mar = c(2,4.5,2,2))

plot(0:(Timesteps_new), C_M_sita_case3, type = "l", lwd = 2, ylab = "Control", col = "blue", lty = 1,  ylim = c(min(C_M_sita_case1, C_M_sita_case2, C_M_sita_case3), max(C_M_sita_case1, C_M_sita_case2, C_M_sita_case3)), cex.lab = 2)
points(0:(Timesteps_new), C_M_sita_case2, type = "l", lwd = 2, col = "purple", cex.lab = 2)
points(0:(Timesteps_new), C_M_sita_case1, type = "l", lty = 5, lwd = 2, col = "red", cex.lab = 2)

plot(0:(Timesteps_new), S_H_sita_case1, type = "l", lwd = 2, ylab = "Susceptible Human", col = "red",lty = 5,  ylim = c(min(S_H_sita_case1, S_H_sita_case2, S_H_sita_case3), max(S_H_sita_case1, S_H_sita_case2, S_H_sita_case3)), cex.lab = 2)
points(0:(Timesteps_new), S_H_sita_case2, type = "l", lwd = 2, col = "purple", cex.lab = 2)
points(0:(Timesteps_new), S_H_sita_case3, type = "l", lwd = 2, col = "blue", cex.lab = 2)
legend("top", legend = c("Initial Terror", expression(paste(~theta, ~"= 0")), expression(paste(~theta, ~"= 1")), expression(paste(~theta, ~"= 1.5"))), lty = c(1, 5, 1 , 1), lwd = c(2, 2, 2, 2), col = c("white", "red", "purple", "blue"), bty = "n", cex = 2)

plot(0:(Timesteps_new), I_H_sita_case1, type = "l", lwd = 2, col = "red",lty = 5, ylim = c(min(I_H_sita_case1, I_H_sita_case2, I_H_sita_case3) - 80, max(I_H_sita_case1, I_H_sita_case2, I_H_sita_case3)), ylab = "Infected Human", cex.lab = 2)
points(0:(Timesteps_new), I_H_sita_case2, type = "l", lwd = 2, col = "purple", cex.lab = 2)
points(0:(Timesteps_new), I_H_sita_case3, type = "l", lwd = 2, col = "blue", cex.lab = 2)
legend("topright", legend = c("log(Tot_Infec)", paste("=",round(log(sum(I_H_sita_case1)), 2)), paste("=", round(log(sum(I_H_sita_case2)), 5)), paste("=", round(log(sum(I_H_sita_case3)), 2))), lty = c(1, 5, 1 , 1), lwd = c(2, 2, 2, 2), col = c("white", "red", "purple", "blue"), bty = "n", cex = 2)

plot(0:(Timesteps_new), D_H_sita_case1, type = "l", lwd = 2, ylab = "Death Cases", col = "red", lty = 5, ylim = c(min(D_H_sita_case1, D_H_sita_case2, D_H_sita_case3),max(D_H_sita_case1, D_H_sita_case2, D_H_sita_case3)), cex.lab = 2)
points(0:(Timesteps_new), D_H_sita_case2, type = "l", lwd = 2, col = "purple3", cex.lab = 2)
points(0:(Timesteps_new), D_H_sita_case3, type = "l", lwd = 2, col = "blue", cex.lab = 2)

dev.off()




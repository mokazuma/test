
##### SEIR model

pacman::p_load(magrittr, ggplot2)

rkode <- function(tmax, init, dt) {
  
  S = init$S; E = init$E; I = init$I; R = init$R
  k0 = k1 = k2 = k3 = vector()
  d = data.frame() %>% dplyr::tbl_df()
  
  for(t in seq(0, tmax-dt, by=dt)) {
    k0[1] = dt * dS(t,S,E,I,R)
    k0[2] = dt * dE(t,S,E,I,R)
    k0[3] = dt * dI(t,S,E,I,R)
    k0[4] = dt * dR(t,S,E,I,R)
    
    k1[1] = dt * dS(t+dt/2, S+k0[1]/2, E+k0[2]/2, I+k0[3]/2, R+k0[4]/2)
    k1[2] = dt * dE(t+dt/2, S+k0[1]/2, E+k0[2]/2, I+k0[3]/2, R+k0[4]/2)
    k1[3] = dt * dI(t+dt/2, S+k0[1]/2, E+k0[2]/2, I+k0[3]/2, R+k0[4]/2)
    k1[4] = dt * dR(t+dt/2, S+k0[1]/2, E+k0[2]/2, I+k0[3]/2, R+k0[4]/2)
    
    k2[1] = dt * dS(t+dt/2, S+k1[1]/2, E+k1[2]/2, I+k1[3]/2, R+k1[4]/2)
    k2[2] = dt * dE(t+dt/2, S+k1[1]/2, E+k1[2]/2, I+k1[3]/2, R+k1[4]/2)
    k2[3] = dt * dI(t+dt/2, S+k1[1]/2, E+k1[2]/2, I+k1[3]/2, R+k1[4]/2)
    k2[4] = dt * dR(t+dt/2, S+k1[1]/2, E+k1[2]/2, I+k1[3]/2, R+k1[4]/2)
    
    k3[1] = dt * dS(t+dt, S+k2[1], E+k2[2], I+k2[3], R+k2[4])
    k3[2] = dt * dE(t+dt, S+k2[1], E+k2[2], I+k2[3], R+k2[4])
    k3[3] = dt * dI(t+dt, S+k2[1], E+k2[2], I+k2[3], R+k2[4])
    k3[4] = dt * dR(t+dt, S+k2[1], E+k2[2], I+k2[3], R+k2[4])
    
    S = S + (k0[1] + 2*k1[1] + 2*k2[1] + k3[1])/6
    E = E + (k0[2] + 2*k1[2] + 2*k2[2] + k3[2])/6
    I = I + (k0[3] + 2*k1[3] + 2*k2[3] + k3[3])/6
    R = R + (k0[4] + 2*k1[4] + 2*k2[4] + k3[4])/6
    
    d = d %>% dplyr::bind_rows(data.frame(Time = t, S, E, I, R))
  }
  
  return(d)
}

##### Parameters
Ro = 2.5          # basic reproduction number: https://twitter.com/ClusterJapan/status/1247463049662889985
latent = 5        # latent period: https://twitter.com/ClusterJapan/status/1252845333366796288
infectious = 10   # infectious period: about...

epsilon = 1 / latent    # rate at which an exposed person becomes infective
gamma = 1 / infectious  # recovery rate
beta = Ro * gamma       # infection rate
init = data.frame(S = 1-1e-5*1, E = 0, I = 1e-5, R = 0.0)

##### Ordinary Differential Equations, Runge-Kutta Method
dS <- function(t, S, E, I, R) r <- -beta * S * I
dE <- function(t, S, E, I, R) r <-  beta * S * I - epsilon * E
dI <- function(t, S, E, I, R) r <-  epsilon * E - gamma * I
dR <- function(t, S, E, I, R) r <-  gamma * I
out = rkode(tmax=250, init=init, dt=0.1)

##### Plot time series
gdat = out %>% dplyr::rename(Susceptible=S, Exposed=E, Infected=I, Recovered=R) %>% 
  tidyr::gather(Index, dv, -Time) %>% 
  transform(Index = factor(Index, levels=c('Susceptible', 'Exposed', 'Infected', 'Recovered'))) %>% 
  dplyr::mutate(dv = 100*dv)
ggplot(gdat, aes(x=Time, y=dv, color=Index)) + geom_line() + 
  theme_bw(base_size = 18) + theme(legend.position = c(0.25, 0.5)) + 
  xlab('Time (day)') + ylab('Proportion (%)') + ggtitle('SEIR Model')

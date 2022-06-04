#Install deSolve
R=2
r=1.5
gamma2=2/3 #gamma squared
omegae=1 #omega_dot_e
g=1
l=5
#Defining the periodic function lambda(t)
t1=l/(r*omegae)
t2=t1+pi*(R-r)/(r*omegae)
lambda = function(t){
  s = t%%t2
  L = (1/(R-r))*(s>t1 & s<=t2)
  return(L)
}
########################################

#Defining the periodic function D(t)
Height = function(t){
  s = t%%(2*t2)
  D = (l*s/t1)*(s<=t1)
  D = D+(l+(R-r)*sin(r*omegae*(s-t1)/(R-r)))*(s>t1 & s<=t2)
  D = D+(l-(l/t1)*(s-t2))*(s>t2 & s<=t1+t2)
  D = D-(R-r)*sin(r*omegae*(s-(t1+t2))/(R-r))*(s>t1+t2)
  return(D)
}
######################################
c1 =-(gamma2/(1+gamma2))*(r^2)*omegae
c2 = omegae
c3 = -g/(1+gamma2)
parameters = c(c1,c2,c3)

Rolling = function(t,y,parameters){
  c1 = parameters[1]
  c2 = parameters[2]
  c3 = parameters[3]
  with(as.list(y),{
    dtime_var     = 1
    dh            = sigma
    dsigma        = c1*lambda(time_var)*w+c3
    dw            = c2*lambda(time_var)*sigma
    list(c(dtime_var,dh,dsigma,dw))
  })
}

yinit = c(time_var=0,h=0,sigma=1,w=1)
times = seq(from=0,to=500,by=0.01)
out   = ode(y=yinit, times=times, func=Rolling, parms=parameters)

plot(out[,"time_var"],out[,"h"],type="l",xlab="time",ylab="height")
grid()

h1=out[,"h"]

#x=seq(from=0,to=30,by=0.01)
#plot(x,Height(x),type='l')
#plot(x,lambda(x),type='l')

plot(h1,Height(times),type="l",xlab="longitude",ylab="height")
grid()



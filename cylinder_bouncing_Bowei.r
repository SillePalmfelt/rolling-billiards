#Cross-product of two vectors
cross_prod = function(u,v){
  w = c(u[2]*v[3]-u[3]*v[2],u[3]*v[1]-u[1]*v[3],u[1]*v[2]-u[2]*v[1])
  return(w)
}

####################################
#PVT Function takes in last a, last post-collision linear and angular velocities u and omega
#force phi, gamma, R, r, m
#Gives back new position a, new postcollision linear and angular velocities,
#and intercollision time.
PVT = function(a,u,omega,phi,gamma,R,r,m){
  s             = 2*gamma/(1+gamma^2)
  c             = (1-gamma^2)/(1+gamma^2)
  a[1:2]        = (R-r)*a[1:2]/sqrt(a[1]^2+a[2]^2)
  nu            = -c(a[1],a[2],0)/(R-r)
  t             = 2*(R-r)*sum(nu*u)/sum(u[1:2]*u[1:2])
  a_new         = a+t*u+(t^2/2)*(phi/m)*c(0,0,1) #new position
  u_n           = u+t*(phi/m)*c(0,0,1) #new pre-coll linear velocity
  omega_n       = omega #new pre-coll angular velocity
  nu_n          = -c(a_new[1],a_new[2],0)/(R-r) #new unit normal vector
  u_new         = c*u_n-(s/gamma)*sum(u_n*nu_n)*nu_n+(s*r*gamma)*cross_prod(omega_n,nu_n)
  omega_new     = (s/(r*gamma))*cross_prod(nu_n,u_n)-c*omega_n+(s/gamma)*sum(omega_n*nu_n)*nu_n
  return(list(a_new,u_new,omega_new,t))
}
#############################
#Experiment
#############################
N          = 200
Position   = matrix(0,nrow=3,ncol=N)
Linear_v   = matrix(0,nrow=3,ncol=N)
Angular_v  = matrix(0,nrow=3,ncol=N)
Coll_times = matrix(0,nrow=1,ncol=N)
#############################
#Parameters
gamma = sqrt(2/5)
R     = 10
r     = 1
m     = 1
phi   = -0.5
############################
#Initial conditions: start position: (R-r)*c(1,0,0)
#                    initial linear velocity:  speed*c(-epsilon,sqrt(1-epsilon^2),v_3)
#                    initial angular velocity: c(omega1,omega2,omega3)                   
Position[,1]   = c((R-r)*cos(0),(R-r)*sin(0),0)
epsilon        = 0.8
Linear_v[,1]   = c(-3*epsilon,3*sqrt(1-epsilon^2),0)
ss=1.7
Angular_v[,1]  = c(0,0,-ss*3*sqrt(1-epsilon^2)/r)
for (i in 2:N){
  a              = Position[,i-1]
  u              = Linear_v[,i-1]
  omega          = Angular_v[,i-1]
  PVT_new        = PVT(a,u,omega,phi,gamma,R,r,m)
  Position[,i]   = PVT_new[[1]]
  Linear_v[,i]   = PVT_new[[2]]
  Angular_v[,i]  = PVT_new[[3]]
  Coll_times[,i] = PVT_new[[4]]+Coll_times[,i-1]
}
#Kinetic energy
KE = (m/2)*(apply(Linear_v*Linear_v,2,sum)+(r*gamma)^2*apply(Angular_v*Angular_v,2,sum))-Position[3,]*m*phi

Pos_X = Position[1,]
Pos_Y = Position[2,]
Pos_Z = Position[3,]

################
#Collision time limit
N_limit = 1102
theta_view=pi/20
X_image=Pos_Y[1:N_limit]
Y_image=-sin(theta_view)*Pos_X[1:N_limit]+cos(theta_view)*Pos_Z[1:N_limit]
plot(X_image,Y_image,type='l',asp=1,xlim=c(-(R-r),(R-r)),ylim=c(min(Y_image)-1,max(Y_image)+3))
X_line1=c(-(R-r),-(R-r))
Y_line1=c(0,cos(theta_view)*min(Pos_Z[1:N_limit]))
lines(X_line1,Y_line1)
X_line2=c((R-r),(R-r))
Y_line2=c(0,cos(theta_view)*min(Pos_Z[1:N_limit]))
lines(X_line2,Y_line2)
angle  = seq(from=0,to=2*pi,by=2*pi/100)
l_angle = length(angle)

X_circ1=(R-r)*sin(angle)
zz=(0*c(1:l_angle)+0)
Y_circ1=-sin(theta_view)*(R-r)*cos(angle)+cos(theta_view)*zz
lines(X_circ1,Y_circ1)
Y_circ2=-sin(theta_view)*(R-r)*cos(angle)+cos(theta_view)*(zz+min(Pos_Z[1:N_limit]))
lines(X_circ1,Y_circ2)
################



#Determine step number for given collision time. Use ss=1.15
#########################
#aa=which(Coll_times>250)
#number=min(aa)
#number
#Coll_times[number]
#########################
#Time = 36.8,  step number = 222
#Time = 73.4,  step number = 442
#Time = 110.0, step number = 662
#Time = 146.6, step number = 882
#Time = 183.2, step number = 1102


#Cross-section path
#########################
angle  = seq(from=0,to=2*pi,by=2*pi/100)
X_circ = (R-r)*cos(angle)
Y_circ = (R-r)*sin(angle)
plot(X_circ,Y_circ,type='l',asp=1)
lines(Pos_X,Pos_Y)
#########################
plot(Coll_times,Pos_Z,type='l')
grid()
points(c(36.8,73.4,110.0,146.6,183.2),c(0,-40,-50,-80,-100))
max(KE)-min(KE)

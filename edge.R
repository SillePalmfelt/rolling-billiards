r=1 #radius of ball
theta=pi/10 #angle between edges
x= 25 #distance from the corner
N = 1000 #count of iterations
c = cos(theta)
s = sin(theta)
tan = s/c
eta= 0.5
E1 = c(0,1) #normal vector of the horizontol edge
E2 = c(s,-c) #normal vector of the slant edge  
e1 = -1*c(1,0) #unit vector for the horizontal edge 
e2 = (-1)*c(c,s) #unit vector for the slant edge
J=matrix(0,nrow=2,ncol=2)  #J is rotation by pi/2
J[1,2]=-1
J[2,1]=1



Roll = function(sheet,edge,a,u,sp){
  sigma= (-1)^(sheet==1) #different sign depending on which side the ball is on
  t = (a[2]-a[1]*tan)/(u[1]*tan-u[2])*(edge==0)-a[2]/u[2]*(edge==1) #time to reach the other end
  a_n = a+t*u #new position after reaching the other edge, but before the displacement on the edge
  edge_new = 1*(edge==0)
  E = E1*(edge_new==0)+E2*(edge_new==1) #normal vector of the edge where rolling occurs
  u0 = u-sum(u[1:2]*E[1:2])*E #tangential component of linear velocity
  W=sp*J%*%E #convert spin to W (as in the paper)
  u0_new = cos(pi*eta)*u0+sin(pi*eta)*W  #new tangential component
  u_n = u0_new+sum(u[1:2]*E[1:2])*E #new linear velocity
  e_n=e1*(edge_new==0)+e2*(edge_new==1)
  u_new=2*sum(u_n[1:2]*e_n[1:2])/sum(u_n[1:2]*u_n[1:2])*e_n-u_n
  W_new = sin(pi*eta)*u0-cos(pi*eta)*W #new W (as in the paper)
  sheet_new =  1*(sheet==0)
  mu = sum(u[1:2]*E[1:2]) 
  tm = pi*r/abs(mu) #time to roll around
  omega = eta*mu/r #omega is as in the paper
  X0=a0-sigma*W0/omega
  X_new=(cos(omega*eta*tm)-1)*W0/omega+sigma*sin(omega*eta*tm)*u1/omega
  dis = sum(X_new[1:2]*E[1:2])  #These three lines are to compute displacement on the edge according to example 9
  JE = J%*%E
  sp_new = sum(W_new[1:2]*JE[1:2]) #new spin
  a_new = a_n-dis*e2*(edge_new==1)*sigma-dis*e1*(edge_new==0)*sigma #new position after the displacement
  t_new = t+tm
  KE=sum(u_new[1:2]*u_new[1:2]+W_new[1:2]*W_new[1:2])
  return(list(a_n,a_new,u_new,t_new,sp_new,sheet_new,edge_new,KE))
}

Position   = matrix(0,nrow=2,ncol=2*N)
Linear_v   = matrix(0,nrow=2,ncol=N)
Spin_v = matrix(0,nrow=1,ncol=N)
Coll_times = matrix(0,nrow=1,ncol=N)
Sheet_side = matrix(0,nrow=1,ncol=N) #keep track of the sheet side
Edge_side = matrix(0,nrow=1,ncol=N) #keep track of the edge
Kinetic_Energy=matrix(0,nrow=1,ncol=N)

Position[,1]   = c(x,0)
epsilon        = 0.03
Linear_v[,1]   = c(0,1)  
Spin_v[,1]  = 1


#####initial conditions, since the displacement depends on the initial conditions

a0            = Position[,1]
u1             = Linear_v[,1] #u0 is taken so I used u1 for initial linear velocity
sp0 = Spin_v[,1]
sheet0         = Sheet_side[,1]
edge0        = Edge_side[,1]
mu0 = sum(u1[1:2]*E1[1:2])*(edge0==0)+sum(u1[1:2]*E2[1:2])*(edge0==1) 
W0=sp0*J%*%E1*(edge0==1)+sp0*J%*%E2*(edge0==0)

for (i in 2:N){
  a              = Position[,2*i-3] #there is an intermediate position a_n, so two slots are needed to record position in each iteration
  u              = Linear_v[,i-1]
  sheet         = Sheet_side[,i-1]
  edge          = Edge_side[,i-1]
  sp = Spin_v[,i-1]
  Roll_new      = Roll(sheet,edge,a,u,sp)
  Position[,2*i-2] = Roll_new[[1]]
  Position[,2*i-1]   = Roll_new[[2]]
  Linear_v[,i]   = Roll_new[[3]]
  Coll_times[,i] = Roll_new[[4]]+Coll_times[,i-1]
  Spin_v[,i] = Roll_new[[5]]
  Sheet_side[,i] = Roll_new[[6]]
  Edge_side[,i] = Roll_new[[7]]
  Kinetic_Energy[,i]=Roll_new[[8]]
}


Pos_X = Position[1,]
Pos_Y = Position[2,] 


X_image=Pos_X[1:(2*N)]
Y_image=Pos_Y[1:(2*N)]
plot(X_image,Y_image,type='l',asp=1,xlim=c(0,50),ylim=c(0,50*tan),xlab="x",ylab="y")
segments(x0 = 0, y0 = 0, x1 = 50, y1 = 50*tan, col = "red")  
segments(x0 = 0, y0 = 0, x1 = 50, y1 = 0, col = "red")  
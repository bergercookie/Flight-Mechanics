% Close old plot so that new plot appears on top
close
%
% Boeing 747, p 165 Etkin
%
% State vector is x = (delta_u  w  q  delta_theta)
%

% US units -----------------------------------------------------

% Geometry

cbar=27.31

% Mass properties

g=32.2
W=636636
m=W/g
Ix=0.183*10^8
Iy=0.331*10^8

% Aerodynamics

Xu=-135.8
Xw=275.8
Zu=-1778
Zw=-6188
Zq=-101700
Zwdot=130.8
Mu=3581
Mw=-35150
Mq=-1.122*10^7
Mwdot=-3826

Xde=-3.717
Zde=-3.551*10^5
Mde=-3.839*10^7

Xdp=0.3*m*g

% Flight condition

u0=774
theta0=0

% The system matrix

A(1,1)=Xu/m
A(1,2)=Xw/m
A(1,3)=0
A(1,4)=-g*cos(theta0)
A(2,1)=Zu/(m-Zwdot)
A(2,2)=Zw/(m-Zwdot)
A(2,3)=(Zq+m*u0)/(m-Zwdot)
A(2,4)=-m*g*sin(theta0)/(m-Zwdot)
A(3,1)=(Mu+Mwdot*Zu/(m-Zwdot))/Iy
A(3,2)=(Mw+Mwdot*Zw/(m-Zwdot))/Iy
A(3,3)=(Mq+Mwdot*(Zq+m*u0)/(m-Zwdot))/Iy
A(3,4)=-Mwdot*m*g*sin(theta0)/(Iy*(m-Zwdot))
A(4,1)=0
A(4,2)=0
A(4,3)=1
A(4,4)=0

% (7.6.4)
B(1,1)=Xde/m
B(2,1)=Zde/(m-Zwdot)
B(3,1)=Mde/Iy + Mwdot*Zde/(Iy*(m-Zwdot))
B(4,1)=0

B(1,2)=0.3*g
B(2,2)=0
B(3,2)=0
B(4,2)=0

% State vector is x = (delta_u  w  q  delta_theta)
% Output vector y = [ delta_u  alfa  gamma q ]

C =[1   0     0  0 
    0  1/u0   0  0
    0 -1/u0   0  1
    0   0     1  0];
D =zeros(4,2);


%-----------------------------------------------------------------
% Sample uses of Matlab commands on the state space model

initial(A,B,C,D,[1 0 0 0],100) 
temp=input('<cr> for next');

t = 0:0.01:50;
u = sin(t);
lsim(A,B(:,1),C,D(:,1),u,t)
temp=input('<cr> for next');

istate=1
bode(A,B,C(istate,:),D(istate,:),1)
temp=input('<cr> for next');

% This sequence gives Figures 7.19 and 7.20

T=0:.01:10;
step(A,B,(pi/180)*C(1,:),D(1,:),1,T);title('Figure 7.19 (a)')
temp=input('<cr> for next');

T=0:.01:60;
step(A,B,(pi/180)*C(2,:),D(2,:),1,T);title('Figure 7.19 (b)')
temp=input('<cr> for next');

T=0:.01:10;
step(A,B,(pi/180)*C(3,:),D(3,:),1,T);title('Figure 7.19 (c)')
temp=input('<cr> for next');

T=0:1:600;
step(A,B,(pi/180)*C(1,:),D(1,:),1,T);title('Figure 7.20 (a)')
temp=input('<cr> for next');

T=0:1:600;
step(A,B,(pi/180)*C(2,:),D(2,:),1,T);title('Figure 7.20 (b)')
temp=input('<cr> for next');

T=0:1:600;
step(A,B,(pi/180)*C(3,:),D(3,:),1,T);title('Figure 7.20 (c)')
temp=input('<cr> for next');

% Response to throttle, Figure 7.21

T=0:1:600;
step(A,B,(1/6)*C(1,:),D(1,:),2,T) ; title('Figure 7.21 (a)')
temp=input('<cr> for next');

step(A,B,(1/6)*C(2,:),D(2,:),2,T) ; title('Figure 7.21 (b)')
temp=input('<cr> for next');

step(A,B,(1/6)*C(3,:),D(3,:),2,T) ; title('Figure 7.21 (c)')
temp=input('<cr> for next');


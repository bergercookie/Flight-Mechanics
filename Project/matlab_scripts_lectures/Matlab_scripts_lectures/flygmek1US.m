%
% Boeing 747, p 165 Etkin
%
% State vector is x = (delta_u  w  q  delta_theta)
%

% US units -----------------------------------------------------

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

cbar=27.31

% Mass properties

g=32.2
W=636636
m=W/g
Ix=0.183*10^8
Iy=0.331*10^8

% Flight condition

u0=774
theta0=0

% Data in SI-units

% Aerodynamics

XuSI=-1982
XwSI=4025
ZuSI=-25950
ZwSI=-90300
ZqSI=-452400
ZwdotSI=1909
MuSI=15930
MwSI=-156300
MqSI=-15210000
MwdotSI=-17020

cbarSI=8.324

% Mass properties

gSI=9.81
WSI=2831760
mSI=WSI/gSI
IxSI=0.247*10^8
IySI=0.449*10^8

% Flight condition

u0SI=235.9
theta0SI=0

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

A_Etkin=[-0.006868    0.01395     0.0     -32.2
         -0.09055    -0.3151    773.98      0.0
          0.0001187  -0.001026   -0.4285    0.0
          0.0         0.0         1.0       0.0]

ASI(1,1)=XuSI/mSI
ASI(1,2)=XwSI/mSI
ASI(1,3)=0
ASI(1,4)=-gSI*cos(theta0SI)
ASI(2,1)=ZuSI/(mSI-ZwdotSI)
ASI(2,2)=ZwSI/(mSI-ZwdotSI)
ASI(2,3)=(ZqSI+mSI*u0SI)/(mSI-ZwdotSI)
ASI(2,4)=-mSI*gSI*sin(theta0SI)/(mSI-ZwdotSI)
ASI(3,1)=(MuSI+MwdotSI*ZuSI/(mSI-ZwdotSI))/IySI
ASI(3,2)=(MwSI+MwdotSI*ZwSI/(mSI-ZwdotSI))/IySI
ASI(3,3)=(MqSI+MwdotSI*(ZqSI+mSI*u0SI)/(mSI-ZwdotSI))/IySI
ASI(3,4)=-MwdotSI*mSI*gSI*sin(theta0SI)/(IySI*(mSI-ZwdotSI))
ASI(4,1)=0
ASI(4,2)=0
ASI(4,3)=1
ASI(4,4)=0

% Dependent state variables x_2 = (delta_xdot_E delta_zdot_E ) = B*x

B(1,1)=cos(theta0)
B(1,2)=sin(theta0)
B(1,3)=0
B(1,4)=-u0*sin(theta0)
B(2,1)=-sin(theta0)
B(2,2)=cos(theta0)
B(2,3)=0
B(2,4)=-u0*cos(theta0)

BSI(1,1)=cos(theta0SI)
BSI(1,2)=sin(theta0SI)
BSI(1,3)=0
BSI(1,4)=-u0SI*sin(theta0SI)
BSI(2,1)=-sin(theta0SI)
BSI(2,2)=cos(theta0SI)
BSI(2,3)=0
BSI(2,4)=-u0SI*cos(theta0SI)

% Solve eigenvalue problem

[V,D]=eig(A)

% Non-dimensional eigenvectors

scale=[1/u0 1/u0 cbar/(2*u0) 1]'
scale=[scale scale scale scale]
Vnond=V.*scale
scale2=abs(Vnond(4,:))
scale2=[scale2
        scale2
        scale2
        scale2]
Vnond=Vnond./scale2
Mag=abs(Vnond)
Phase=(180/pi)*angle(Vnond)
scale3=Phase(4,:)
scale3=[scale3
        scale3
        scale3
        scale3]
Phase=Phase-scale3

eig(A)
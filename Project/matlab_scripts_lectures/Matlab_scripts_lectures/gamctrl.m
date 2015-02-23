% Boeing 747 Longitudinal example from Etkin

close all
theta0 = 0.0;
u0     = 774.0;
A = [ -0.00686854   0.01394951    0.0           -32.2   0.0
      -0.09052721  -0.31506319  773.97653828      0.0   0.0
       0.00011865  -0.00102552   -0.42843608      0.0   0.0
       0.0          0.0           1.0             0.0   0.0
       0.0          -1/u0         0.0             1.0   0.0 ];
Bde = [  -0.000187
      -17.85
       -1.158
        0.0       
        0.0        ];

% Throttle input

Bdp= [9.66 0 0 0 0]'


[V,D]=eig(A(1:4,1:4));


T=2*pi./imag(diag(D));
thalf=log(2)./abs(real(diag(D)));
nhalf=(log(2)/(2*pi))*abs(imag(diag(D)))./abs(real(diag(D)));

% Set up linear time invariant system (LTI)

k=[-0.5 -0.5 -0.5];
k=[20 1 300];

Chat=[0 -1/u0 0 1 0
      0 0 0 0 1];
Chat(3,:)=-A(2,1:5)/u0+[0 0 1 0 0];

Ahat=A-Bdp*k*Chat

Bhat=k(1)*Bdp+[0 0 0 0 -1]';

sysdp=ss(Ahat,(pi/180)*Bhat,(180/pi)*[0 -1/u0 0 1 0],[0]);

step(sysdp,[0:300])

break

% Root locus plot

k=[0 0 0];
Ahat=A-B*k*Chat;
ev0(:,1)=eig(Ahat);
for j=1:51
Ahat=A-B*k*Chat;
ev1(:,j)=eig(Ahat);
k(1)=k(1)-0.01;
end
k=[-0.5 0 0];
for j=1:151
Ahat=A-B*k*Chat;
ev2(:,j)=eig(Ahat);
k(2)=k(2)-0.01;
end
k=[-0.5 -1.5 0];
for j=1:51
Ahat=A-B*k*Chat;
ev3(:,j)=eig(Ahat);
k(3)=k(3)-0.01;
end
close all
plot(real(ev1),imag(ev1),'k.',real(ev0),imag(ev0),'ko')
axis([-0.5 1 0 1.5])
xlabel('real(\lambda)')
ylabel('imag(\lambda)')
title('Root locus plot of attitude controller')
hold on
plot(real(ev2),imag(ev2),'b.')
plot(real(ev3),imag(ev3),'r.')

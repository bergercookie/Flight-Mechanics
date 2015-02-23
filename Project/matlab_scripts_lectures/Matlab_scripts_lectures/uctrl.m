% Boeing 747 Longitudinal example from Etkin

close all

theta0 = 0.0;
u0     = 774.0;
Ahat = [ -0.00686854   0.01394951    0.0           -32.2   0.0
         -0.09052721  -0.31506319  773.97653828      0.0   0.0
          0.00011865  -0.00102552   -0.42843608      0.0   0.0
          0.0          0.0           1.0             0.0   0.0
          1.0          0.0           0.0             0.0   0.0 ];
Bhat = [  -0.000187   9.66
         -17.85       0.0
          -1.158      0.0
           0.0        0.0
           0.0        0.0  ];

[V,D]=eig(Ahat(1:4,1:4));

T=2*pi./imag(diag(D));
thalf=log(2)./abs(real(diag(D)));
nhalf=(log(2)/(2*pi))*abs(imag(diag(D)))./abs(real(diag(D)));

% Set up linear time invariant system (LTI)

% PID control for following commanded delta_u = 0

k=[0.005 0.0001 0.015]

%k=[0.0108 .0045 0.018]
% Ziegler Nichols
% k0=0.018 T0=8
%k=[0.0081  0.0027 0]

Chat=[1 0 0 0 0
      0 0 0 0 1
      Ahat(1,:) ];

Aloop=Ahat-Bhat(:,1)*k*Chat

Bloop=(-1/6)*Bhat(:,2);

C = [1    0   0  0  0        % u
     0 -1/u0  0  1  0        % -alpha+theta=gamma
        -k*Chat              % delta_e
     0  1/u0  0  0  0   ];   % alpha

D=zeros(1,1);

sys=ss(Aloop,Bloop,C,D);

step(sys,50)
%step(sys,100)
%figure
%bode(sys,{0.01 10})
%figure
%impulse(sys)


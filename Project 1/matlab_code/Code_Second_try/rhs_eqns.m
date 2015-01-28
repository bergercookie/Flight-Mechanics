function [y_out] = rhs_eqns(t, y)

%% Variable Definitions
global sref tepsr xtpref xtpfpl rmfpl rmfull
global rmfuel iebk
global g

%% "Proper" Variable Naming
v = y(1); h = y(2); m = y(3);
h_km = h / 1000; %altitude in km

%% Aerodynamic Properties
[rho,aspeed,temp,press] = stdatm(h_km); %given altitude compute these properties
mach = v / aspeed;
[T,b]=rm6cal(mach,h_km,iebk); %compute thrust and fuel burn

qdyn=0.5*rho*v^2;  % Dynamic pressure

% Use the aerodynamic data base to extract the coefficients

ca0=ca0cal(mach);
cd0=cd0cal(mach);
cla=clacal(mach);
clarad=cla*180/pi;
[eta,dxeta]=etacal(mach);
xtp=xtpcal(rmfuel);
eta=eta+(xtp-xtpref)*dxeta;

%% Gamma Function
% gamma = gamma_func(t); %todo
gamma = 0;

%% Solving for alfa
alfa = find_a(T, m, gamma); alfadeg = to_degrees(alfa);

%% Drag
cl=cla*(alfadeg-ca0);
D  = ( cd0 + sref*eta*(cl)^2 ) * qdyn;

%% Equations Configuration
y_calc = zeros(3,1);
y_calc(1) = (T * cos(alfa + tepsr) - D - m*g*sin(gamma)) / m;
y_calc(2) = v * sin(gamma);
y_calc(3) = - b;

%% Return the derivative vector
y_out = y_calc;

end
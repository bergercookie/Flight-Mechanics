% Demo script for J35 Draken model

global sref tepsr xtpref xtpfpl rmfpl rmfull qmax

% Initialize data structures

initfpl35
initrm6

g=9.81;

iebk=1;            % Engine at full afterburner
rmfuel=0.6*rmfull; % Fuel tank is 60% full
gamma=0;           % Level flight
mass=rmfpl+rmfuel  % Total mass of aircraft

mach=0.5;  % Set the Mach number or compute it from some other airspeed
altkm=5.0; % Set the altitude

alt=altkm*1000;  % Maintain also the altitude in meters
[rho,aspeed,temp,press] = stdatm(altkm); % ISA 0 data

qdyn=0.5*rho*(aspeed*mach)^2;  % Dynamic pressure

% Use the aerodynamic data base to extract the coefficients

ca0=ca0cal(mach);
cd0=cd0cal(mach);
cla=clacal(mach);
clarad=cla*180/pi;
[eta,dxeta]=etacal(mach);
xtp=xtpcal(rmfuel);
eta=eta+(xtp-xtpref)*dxeta;

% Use engine model to get the thrust (Newtons) and fuelburn (kg/s)

[thrust,fuelb]=rm6cal(mach,altkm,iebk);

% Use fixed angle of attack or insert code to find equilibrium value

alfa=0.0873;            % Radians
alfadeg=alfa*180/pi;  % Degrees

% Lift and drag

cl=cla*(alfadeg-ca0);
rlift = cl*qdyn*sref;
drag  = ( cd0 + sref*eta*(cl)^2 ) * qdyn;


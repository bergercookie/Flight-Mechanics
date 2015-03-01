function [Tex] = texcess35(altkm, mach, mass, gamma)
% function used to calculate the Excess thrust for the given circumstance.
% TEXCESS35 takes as arguements the altitude [km], the mach number, the mass [kg], the flight
% math angle [rad] and returns the excess thrust [N].

global sref tepsr xtpref xtpfpl rmfpl rmfull qmax
global e_t_rad rmfuel alfa alfadeg iebk
global g

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

% Lift and drag

cl=cla*(alfadeg-ca0);
% rlift = cl*qdyn*sref;
rlift = mass * g; % load factor = 1
drag  = ( cd0 + sref*eta*(cl)^2 ) * qdyn;

Tex = thrust * cos(alfa + e_t_rad) - drag - mass * g * sin(gamma);
end

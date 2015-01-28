function [Tex, Sep, alfa] = TexSep(altkm, mach, mass, gamma)
% function used to calculate the Excess thrust for the given circumstance.
% TexSep takes as arguements the altitude [km], the mach number, the mass [kg], the flight
% math angle [rad] and returns the excess thrust [N].


%% variable definitions
global sref tepsr xtpref xtpfpl rmfpl rmfull qmax
global rmfuel iebk
global g
global qdyn_max qdyn_index qdyn_mach qdyn_alt qdyn_tex


%% aerodynamic properties
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

%% Lift and drag

% angle of attack is not constant
alfa = find_a(thrust, mass, gamma); alfadeg = to_degrees(alfa);

cl=cla*(alfadeg-ca0);
rlift = cl*qdyn*sref;

drag  = ( cd0 + sref*eta*(cl)^2 ) * qdyn;
    
vel = mach * aspeed;

%% output variables
Tex = thrust * cos(alfa + tepsr) - drag - mass * g * sin(gamma);
Sep = Tex * vel / (mass * g);
alfa_out = alfa;

%% dynamic pressure limit
if qdyn > qdyn_max
    qdyn_alt(qdyn_index) = altkm;
    qdyn_mach(qdyn_index) = mach;
    qdyn_tex(qdyn_index) = Tex;
    qdyn_index = qdyn_index + 1;
end

end

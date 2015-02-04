function cd0=cd0cal(m)
% Zero lift drag coefficient vs Mach
% Compute CD_0*Sref depending on Mach number
% NOTE: The coefficient given is dimensional (m^2)!

global cd0dat
mtem=max(0,min(2,m));
cd0=interp1(cd0dat(:,1),cd0dat(:,2),mtem);


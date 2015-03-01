function cla=clacal(m)

% Compute CL_alfa depending on Mach number
% NOTE: alfa in degrees!

global cladat
mtem=max(0,min(2,m));
cla=interp1(cladat(:,1),cladat(:,2),mtem);


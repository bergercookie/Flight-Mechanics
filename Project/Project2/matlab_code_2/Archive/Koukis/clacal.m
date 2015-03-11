function cla=clacal(m)

% Compute CL_alfa depending on Mach number
% NOTE: alfa in degrees!

global cladat
cladat=[
0.00  0.037 
0.50  0.037 
1.00  0.035 
1.20  0.031 
1.50  0.027 
2.00  0.022];
mtem=max(0,min(2,m));
cla=interp1(cladat(:,1),cladat(:,2),mtem);


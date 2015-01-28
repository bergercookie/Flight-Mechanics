function ca0=ca0cal(m)
% Zero lift angle of attack as a function of the Mach number
global ca0dat
mtem=max(0,min(2,m));
ca0=interp1(ca0dat(:,1),ca0dat(:,2),mtem);


function ca0=ca0cal(m)
% Zero lift angle of attack as a function of the Mach number
global ca0dat
ca0dat=[
0.0D0  1.70
0.50   1.10
1.00   0.60
1.20   0.40
1.50   0.20
2.00   0.10];
mtem=max(0,min(2,m));
ca0=interp1(ca0dat(:,1),ca0dat(:,2),mtem);


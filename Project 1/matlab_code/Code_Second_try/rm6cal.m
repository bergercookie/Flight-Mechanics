function [thrust,fuelb]=rm6cal(mach,altkm,iebk)

global rm6 mrm6 arm6

% iebk = 0 is full thrust no afterburner ISA 0
% iebk = 1 is full thrust with afterburner ISA 0
% iebk = 2 is full thrust with afterburner ISA -15

% Fix input so that it is in the feasible range

mtem=min(2.1,max(0,mach));
atem=min(16.0,max(0,altkm));
item=min(2,max(0,iebk));

% Compute addresses depending on iebk input value
switch item
case 0
it1=1; ib1=52;
case 1
it1=18; ib1=69;
case 2
it1=35; ib1=86;
end

% Linear interpolation

thrust=interp2(mrm6,arm6,rm6(it1:it1+16,1:22),mtem,atem);
fuelb=interp2(mrm6,arm6,rm6(ib1:ib1+16,1:22),mtem,atem);

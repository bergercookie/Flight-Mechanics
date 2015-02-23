% Estimate torsional stiffness of the flexible sting
%
L=0.81;         % Length in meters
d=0.01;         % Diameter in meters
E=206*1000^3;   % Modulus of elasticity for this steel type
nu=0.3;         % Poisson ratio
G=E/(2*(1+nu)); % Shear modulus
K=pi/2*(d/2)^4; % 
k=(G*K/L);      % Estimated stiffness
kexp=98.2951;   % From an experiment

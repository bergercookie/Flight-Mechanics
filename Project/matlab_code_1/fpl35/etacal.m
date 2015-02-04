function [eta,dxeta]=etacal(m)

% Compute coefficient for induced drag depending on Mach number
% Use data as: eta+(xtp-xtpref)*dxeta

global etadat
mtem=max(0,min(2,m));
eta=interp1(etadat(:,1),etadat(:,3),mtem);

dxeta=(interp1(etadat(:,1),etadat(:,4),mtem) - ...
       interp1(etadat(:,1),etadat(:,2),mtem) )/0.1;


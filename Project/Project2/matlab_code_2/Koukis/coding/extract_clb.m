function clb = extract_clb(y2, clp, aoa)
% EXTRACT_CLB computes the Clb coefficient given the quadratic coefficient 
% of the v-wm curve 
% calculated throught the least squares method. The formula for the
% Wn^2 is given in the report of the Project Work

% global variable definition
global sref span Ix
% calculate rho given the airspeed
altkm = 0; % zero altitude - sea level
[rho, ~, ~, ~] = stdatm(0);

clb = -2*Ix*y2/(aoa*rho*span*sref) - clp^2*sref*span^3*Ix/(32*Ix^2*aoa*span);
end
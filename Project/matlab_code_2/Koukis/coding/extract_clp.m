function clp = extract_clp(slope)
% EXTRACT_CLP computes the Clp coefficient given the slope of the v-n curve 
% calculated throught the least squares method. The formula for the damping 
% coefficient is given in the report of the Project Work

% global variable definition
global sref span Ix

% calculate rho given the airspeed
altkm = 0; % zero altitude - sea level
[rho, ~, ~, ~] = stdatm(0);

clp = 8*Ix*slope/(rho*sref*span^2);
end
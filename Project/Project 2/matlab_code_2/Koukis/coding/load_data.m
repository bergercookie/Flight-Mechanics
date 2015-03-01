function [time_out, vibration_out] = load_data(group, fname)
% LOAD_DATA reads the data from the filename specified (in the specified
% group) computes the necessary properties (see below) and returns the
% sampling times and the vibration data
% 
% Input variables:
%   group -> 'A', 'B', ... (no 'group')
%   a00_u05, a00_15... (no suffix)
% Local variables:
%   offset -> offset from the zero reference line
% Output Variables:
%   time_out -> meaningfull time range of data (start -end of oscillation)
%   vibration_out -> actual data corresponding to time_out
%% variable initialization
global path_to_experimental suffix

%% actual computation

% extract the data
full_path = strcat(path_to_experimental, 'group', group, '/', fname, suffix);
half_path = strcat('group', group, '/', fname, suffix);

results=load(full_path);
time=results(:,1);
vibration=results(:,2);

% define the offset as the mean value of the oscillation
offset = mean(vibration);
vibration = vibration - offset;

% extract meaningfull range
[amplitude,debut]=max(vibration);
for j =1:length(vibration)
    if vibration(j) >= amplitude/8
        fin = j-1;
    end
end

time = time(debut:fin);   
vibration = vibration(debut:fin);

% return function outputs
time_out = time;
vibration_out = vibration;

% display message to make sure
disp(sprintf('File: %s loaded successfully', half_path));
end
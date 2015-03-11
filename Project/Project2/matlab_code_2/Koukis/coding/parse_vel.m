function vel = parse_vel(str)
% PARSE_VEL parses the velocity data out of the string given.
% Essentially takes the integer which is between u and .log
% or whatever is after u if the .log suffix is not given

str1 = 'u';
str2 = '.log';

ind1 = find(str == str1);
wo_suffix = length(str) - length(str2);

if strfind(str, str2)
    % .log suffix given
    vel = str(ind1+1:wo_suffix);
else
    % clear name given
    vel = str(ind1+1:end);
end
vel = strrep(vel, '_', '.');
vel = str2num(vel);
end
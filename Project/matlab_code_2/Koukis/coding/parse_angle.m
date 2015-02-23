function a = parse_angle(str)
% PARSE_ANGLE parses the angle data out of the string given.
% Essentially takes the value which is between 'a' and 'u'.


str1 = 'a';
str2 = 'u';

ind1 = find(str == str1);
ind2 = find(str == str2);

angle = str(ind1+1:ind2-1);
% replace _ with .
angle = strrep(angle, '_', '.');

a = str2num(angle);
end
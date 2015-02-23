function isthere = strfind2(str, pattern1, pattern2)
% STRFIND2 finds whether a0 or a00 exists in a given string
% the need for this function is because matlab takes []|1 as 1 (?!) check
% it!
%todo
case1 = strfind(str, pattern1);
case2 = strfind(str, pattern2);

if case1 == 1 | case2 == 1
    isthere = 1;
else
    isthere = 0;
end

end
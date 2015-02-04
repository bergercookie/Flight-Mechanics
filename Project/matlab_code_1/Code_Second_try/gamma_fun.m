function [gamma] = gamma_fun(t)
%Given the time of the simulation GAMMA_FUN computes the needed gamma value
global gammavec
gamma = interp1(gammavec(:, 1), gammavec(:, 2), t);

end
function [gamma] = gamma_fun(t)
%Given the time of the simulation GAMMA_FUN computes the needed gamma value
gammavec = [0. 0;
    100, 0;
    200, 0.2;
    250, 0.2;
    300, 0;
    10000, 0];

gamma = interp1(gammavec(:, 1), gammavec(:, 2), t);
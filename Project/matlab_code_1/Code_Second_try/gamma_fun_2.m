function [gamma] = gamma_fun_2(t)
%Given the time of the simulation GAMMA_FUN computes the needed gamma value
global gamma_f t_s t_f gamma_s
gamma = interp1([t_s; t_f], [gamma_s; gamma_f], t);

end
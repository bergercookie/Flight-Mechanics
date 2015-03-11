function y = quadratic_model(x, td)
% QUADRATIC_MODEL provides a function which is used for approximating the
% of velocity damping pairs calculated
% func = a - b*t
% func = a*t^2 + b*t + c;
% x(1) = a
% x(2) = b
% x(3) = c

y =  x(1).*td.^2 + x(2).*td + x(3);
end

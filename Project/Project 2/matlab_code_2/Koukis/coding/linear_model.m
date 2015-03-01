function y = linear_model(x, td)
% LINEAR_MODEL provides a function which is used for approximating the
% of velocity damping pairs calculated
% func = a - b*t
% x(1) = y(x = 0)
% x(2) = slope

% y = x(4) * exp(x(2)*td) .* sin(x(1)*td + x(3));
y = x(1) - x(2) .* td;
end

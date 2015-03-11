function y = curve_model(x, td)
% CURVE_MODEL provides a function which is used for approximating the set
% of data time - vibration given.
% func = amp*exp(lamda*t+phase) but no complex numbers inserted
%   x(1)= omega
%   x(2)= damping
%   x(3)= phase
%   x(4)=amplitude
y = x(4) * exp(x(2)*td) .* sin(x(1)*td + x(3));

end

function y = myfun(x, td)
    y = x(4) * exp(x(2)*td) .* sin(x(1)*td + x(3));
end

% func = amp*exp(lamda*t+phase) but no complex numbers inserted
% x(1)= omega
% x(2)= damping
% x(3)= phase
% x(4)=amplitude
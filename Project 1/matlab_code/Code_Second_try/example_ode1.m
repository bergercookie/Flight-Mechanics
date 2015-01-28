function [y_dot] = rhs(t, y)
y_dot = [+ 8 * y(1) + 15*y(2); ...
    -3*y(1)  - 2 * y(2)];
% function [y_out] = rhs(t, y)
% y_out = -0.2*y;
% end
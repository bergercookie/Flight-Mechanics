function [J_out] = linearize_mdl(x, step)
% the purpose of LINEARIZE_MDL function is to linearize the set of
% non-linear algebraic equations given.

for j=1:length(x)
xh=x;
xmh=x;
xh(j)=xh(j)+step;
xmh(j)=xmh(j)-step;
[xdoth]=fplmod(0,xh);
[xdotmh]=fplmod(0,xmh);
J(:,j)=(xdoth'-xdotmh')/(2*step);
end

% return the final Jacobian (return it once)
J_out = J;
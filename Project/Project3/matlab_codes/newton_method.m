function [xtrim_new] = newton_method(xtrim_old, funs, J, ifun, mach)
% NEWTON_METHOD is an implementation of the Newton-Raphson numerical
% approximation method *for the trim problem*.

% INPUTS:
%    xtrim_old: set of variables for which I solve for
%    funs: set of (nonlinear) funs (xdot in the code given)
%    J: Jacobian calculated through the linearization procedure
%    ifun: boolean array indicating the subset of the funs which I have to
%    set to zero
%    mach: Current Mach number

% OUTPUT:
%    xtrim_new: trimmed state at which the newton method converges to

ftrim=funs(ifun)';
ftrim(5)=ftrim(5)-mach;
Jtrim=J(ifun,ivar);
xtrim_new=xtrim_old-Jtrim\ftrim;
fprintf('|f| %e\n',norm(ftrim));
end
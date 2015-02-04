function xtp=xtpcal(fuelmass)

% Compute actual c.g position for this fuel level

global xbmdat rmfpl xtpfpl

fmtem=max(0,min(2323,fuelmass));
xbm=interp1(xbmdat(:,1),xbmdat(:,2),fmtem);

xtp=(rmfpl*xtpfpl+xbm)/(rmfpl+fmtem);

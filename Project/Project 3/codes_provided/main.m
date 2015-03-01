% Main program
% Remember that this model is for low speed but high alfa

% Initialize global data structures

initfpl35
initrm6

% Call model for trimmed point

x=[
       138.126782239696    % u (m/s)
       14.5298136566714    % w (m/s)
       0.0                 % q (rad/s)
       0.104806425361396   % theta (rad)
       0.0                 % Distance (m)
       1000                % alt (m)
       500                 % fuel (kg)
       -0.0544410758220454 % de (rad)
       0.288038080583802   % dp
]

[xdot]=fplmod(0,x)

% Example uses illustrated below

% Linearization

step=10^(-5)
for j=1:length(x)
xh=x;
xmh=x;
xh(j)=xh(j)+step;
xmh(j)=xmh(j)-step;
[xdoth]=fplmod(0,xh);
[xdotmh]=fplmod(0,xmh);
J(:,j)=(xdoth'-xdotmh')/(2*step);
end

% Trim loop

ivar=[1 2 4 8 9];   % Selects variables to change during trim iteration
ifun=[1 2 3 6 11];  % Selects functions that should be zero in trim

xtrim=[100 10 0 0 0]'; % Initial guess

machset=0.30; % In this case Mach number is set

% Iterate until convergence

for iter=1:5

x(ivar)=xtrim;
[xdot]=fplmod(0,x);

step=10^(-5);
for j=1:length(x)
xh=x;
xmh=x;
xh(j)=xh(j)+step;
xmh(j)=xmh(j)-step;
[xdoth]=fplmod(0,xh);
[xdotmh]=fplmod(0,xmh);
J(:,j)=(xdoth'-xdotmh')/(2*step);
end

ftrim=xdot(ifun)';
ftrim(5)=ftrim(5)-machset;
Jtrim=J(ifun,ivar);
xtrim=xtrim-Jtrim\ftrim;
fprintf('|f| %e\n',norm(ftrim));

end % End of Newton iteration




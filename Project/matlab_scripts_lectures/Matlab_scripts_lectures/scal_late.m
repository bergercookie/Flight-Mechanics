% B747
load B747late
u0=774
b=195.7
theta0=0;
vscal=[1/u0 b/(2*u0)  b/(2*u0) 1];
[V D]=eig(ALATE)
for j=1:4
v(:,j)=V(:,j).*vscal'
v(:,j)=v(:,j)/norm(v(:,j))
end

% Extended system

A=[ALATE zeros(4,1)
0 0 sec(theta0) 0 0]

[V D]=eig(A)

vscal=[1/u0 b/(2*u0)  b/(2*u0) 1 1];
for j=1:5
ve(:,j)=V(:,j).*vscal'
ve(:,j)=ve(:,j)/norm(ve(:,j))
end

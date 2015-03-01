clear all
load vibdata.mat

t=t0(4000:6000)-4;   
y=x0(4000:6000);

tsample=t(2)-t(1);
data=iddata(y,[],tsample);
model=n4sid(data);
cmodel=d2c(model);
A=cmodel.A; 

plot(t,y)
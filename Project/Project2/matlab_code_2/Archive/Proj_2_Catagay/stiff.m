clear all, close all, clc
global ca0dat cladat
global sref tepsr xtpref xtpfpl rmfpl rmfull qmax gamvec1 gamvec2
% Estimate torsional stiffness of the flexible sting
%
L=0.81;         % Length in meters
d=0.01;         % Diameter in meters
E=206*1000^3;   % Modulus of elasticity for this steel type
nu=0.3;         % Poisson ratio
G=E/(2*(1+nu)); % Shear modulus
K=pi/2*(d/2)^4; % 
span= 9.4;
k=(G*K/L);      % Estimated stiffness
kexp=98.2951;   % From an experiment

%% Ground Vibration Test

% extract the data
results=load('../groupD/a00_u0.log');
offset=2.5;
time=results(:,1);
vibration=results(:,2)-offset;

% extract meaningfull range
[amplitude,debut]=max(vibration);
for j =1:length(vibration)
    if vibration(j) >= amplitude/8
        fin = j-1;
    end
end

time = time(debut:fin);   
vibration = vibration(debut:fin);

figure(1)
plot(time,vibration)

te= time(2)-time(1);
fe=1/te;
L=length(vibration);
NFFT = 2^nextpow2(L);
Y = fft(vibration,NFFT)/L;
f = fe/2*linspace(0,1,NFFT/2+1);

[~,index_f]=max(2*abs(Y(1:NFFT/2+1)));
freq=f(index_f);
omega=2*pi*freq;

coeff0= [omega; 0; 0; amplitude];
coeff= lsqcurvefit(@myfun,coeff0,time,vibration);
omega=coeff(1);
damping=coeff(2);
phase=coeff(3);
amplitude=coeff(4);
y_fit = amplitude.* exp(damping*time).* sin(omega*time + phase);

figure(2)
plot(time,vibration,'b',time,y_fit,'r')
legend('exp data','fit curve')
xlabel('time (s)')

% Mass moment of inertia
Ix=k/(omega^2+damping^2);

%% Clp and Clb


[rho,aspeed,temp,press] = stdatm(0);
aspeed = round(aspeed);
alpha = [zeros(8,1)' 10*ones(6,1)' 15*ones(6,1)' 20*ones(6,1)' 2*ones(6,1)' 4*ones(6,1)' 6*ones(6,1)' 8*ones(6,1)'];
v = [10:5:40 5 10:5:30 5 10:5:30 5 10:5:30 5 10:5:30 5 10:5:30 5 10:5:30 5 10:5:30 5];
mach = v/aspeed;

the_path = '../groupD';
contents = dir(the_path); % get its contents

C_lp = zeros(1, length(contents) - 3);
% for all the logfiles in the log folder
for n = 3 : length(contents) - 1

    name = [the_path, '/', contents(n).name];
    results = load(name);
    
    time=results(:,1);
    vibration=results(:,2)-offset;
    
    % Extraction of the effective measurment
    [amplitude,debut]=max(vibration);
    for j =1:length(vibration)
        if vibration(j) >= amplitude/8
            fin = j;
        end
    end
    time = time(debut:fin);   
    vibration = vibration(debut:fin);
    
    % Fourier Transform on the signal
    te= time(2)-time(1); % sampling time
    fe=1/te; % sampling frequency
    L=length(vibration);
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    Y = fft(vibration,NFFT)/L;
    f = fe/2*linspace(0,1,NFFT/2+1);
    
    coeff0= [omega; 0; 0; amplitude];
    coeff= lsqcurvefit(@myfun,coeff0,time,vibration);
    Omega=coeff(1);
    Damping=coeff(2);
    Phase=coeff(3);
    Amplitude=coeff(4);
    
    Mach = mach(n-2);
    Alpha = alpha(n-2);
    
    ca0=ca0cal(Mach);
    cla=clacal(Mach);
    cl=cla*(Alpha-ca0);
    
    c_lp= cl/span;
%     C_betha(n-3,1)=
    C_lp(n-2) = c_lp;
    
end

figure(5)
plot(alpha(find(v==5)),C_lp(find(v==5)),'r*')
hold on
plot(alpha(find(v==10)),C_lp(find(v==10)),'b*')
hold on
plot(alpha(find(v==15)),C_lp(find(v==15)),'c*')
hold on
plot(alpha(find(v==20)),C_lp(find(v==20)),'m*')
hold on
plot(alpha(find(v==25)),C_lp(find(v==25)),'k*')
hold on
plot(alpha(find(v==30)),C_lp(find(v==30)),'y*')
hold on
plot(alpha(find(v==35)),C_lp(find(v==35)),'g*')
hold on
plot(alpha(find(v==40)),C_lp(find(v==40)),'+')
legend('5 m/s','10 m/s','15 m/s','20 m/s','25 m/s','30 m/s','35 m/s','40 m/s')
xlabel('Angle of attack'); ylabel('Clp');

% C_lp=horzcat(C_lp(1:8),[C_lp(9:14);0],[C_lp(15:20);0],[C_lp(21:26);0],[C_lp(27:32);0],[C_lp(33:38);0],[C_lp(39:44);0],[C_lp(45:50);0]);
% C_betha=horzcat(C_betha(1:8),[C_betha(9:14);0],[C_betha(15:20);0],[C_betha(21:26);0],[C_betha(27:32);0],[C_betha(33:38);0],[C_betha(39:44);0],[C_betha(45:50);0]);
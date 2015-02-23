clear all
close all
clc

% Estimate torsional stiffness of the flexible sting
%
L=0.81;         % Length in meters
d=0.01;         % Diameter in meters
E=206*1000^3;   % Modulus of elasticity for this steel type
nu=0.3;         % Poisson ratio
G=E/(2*(1+nu)); % Shear modulus
K=pi/2*(d/2)^4; % 
k=(G*K/L);      % Estimated stiffness
kexp=98.2951;   % From an experiment


% Experimental data
% load vibdata
% time=t20;
% vibration=x20;

figure(1)
results=load('a0_u0.log');
offset=2.5;
time=results(:,1);
vibration=results(:,2)-offset;

plot(time,vibration)
xlabel('time (s)')
ylabel('signal')

% Extraction of the effective measurment
[amplitude,debut]=max(vibration);
for j =1:length(vibration)
    if vibration(j) >= amplitude/10
        fin = j;
    end
end
time = time(debut:fin);   
vibration = vibration(debut:fin);
figure
plot(time,vibration)

% Fourier Transform on the signal
te= time(2)-time(1); % sampling time
fe=1/te; % sampling frequency
L=length(vibration);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(vibration,NFFT)/L;
f = fe/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
figure
plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

% Omega
[~,index_f]=max(2*abs(Y(1:NFFT/2+1)));
freq=f(index_f)
omega=2*pi*freq

% Damping and Phase Shift
coeff0= [omega; 0; 0; amplitude];
coeff= lsqcurvefit(@myfun,coeff0,time,vibration);
omega=coeff(1)
damping=coeff(2)
phase=coeff(3)
amplitude=coeff(4)
y_fit = amplitude.* exp(damping*time).* sin(omega*time + phase);

% Comparaison between exp and fit data
figure
plot(time,vibration,'b',time,y_fit,'r')
legend('exp data','fit curve')
xlabel('time (s)')

% Mass moment of inertia
Ix=k/(omega^2+damping^2)

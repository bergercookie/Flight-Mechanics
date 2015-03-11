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


%% GVT

results=load('groupD/a00_u0.log');
offset=2.5;
time=results(:,1);
vibration=results(:,2)-offset;

% figure(1)
% plot(time,vibration)
% xlabel('time (s)')
% ylabel('signal')

% Extraction of the effective measurment
[amplitude,debut]=max(vibration);
for j =1:length(vibration)
    if vibration(j) >= amplitude/8
        fin = j;
    end
end
time = time(debut:fin);   
vibration = vibration(debut:fin);

% figure(2)
% plot(time,vibration)

% Fourier Transform on the signal
te= time(2)-time(1); % sampling time
fe=1/te; % sampling frequency
L=length(vibration);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(vibration,NFFT)/L;
f = fe/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
% figure
% plot(f,2*abs(Y(1:NFFT/2+1))) 
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')

% Omega
[~,index_f]=max(2*abs(Y(1:NFFT/2+1)));
freq=f(index_f)
omega=2*pi*freq

% Damping and Phase Shift
coeff0= [omega; 0; 0; amplitude];
coeff= lsqcurvefit(@myfun,coeff0,time,vibration);
omega=coeff(1);
damping=coeff(2);
phase=coeff(3);
amplitude=coeff(4);
y_fit = amplitude.* exp(damping*time).* sin(omega*time + phase);

% Comparaison between exp and fit data
% figure
% plot(time,vibration,'b',time,y_fit,'r')
% legend('exp data','fit curve')
% xlabel('time (s)')

% Mass moment of inertia
Ix=k/(omega^2+damping^2)



%% C_lp and C_betha

group = 'groupD';
names = ls(group);

C_lp=zeros(length(names)-3,1);
C_betha=zeros(length(names)-3,1);
    
for n = 4 : size(names,1)

    name = [group,'/',names(n,:)];
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
    
%     C_lp(n-3,1)=
%     C_betha(n-3,1)=
    
end

C_lp=horzcat(C_lp(1:8),[C_lp(9:14);0],[C_lp(15:20);0],[C_lp(21:26);0],[C_lp(27:32);0],[C_lp(33:38);0],[C_lp(39:44);0],[C_lp(45:50);0]);
C_betha=horzcat(C_betha(1:8),[C_betha(9:14);0],[C_betha(15:20);0],[C_betha(21:26);0],[C_betha(27:32);0],[C_betha(33:38);0],[C_betha(39:44);0],[C_betha(45:50);0]);








  

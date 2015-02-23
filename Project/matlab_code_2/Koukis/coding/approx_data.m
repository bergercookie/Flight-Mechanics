function [omega, damping, phase, amplitude] = approx_data(time, vibration)
% APPROX_DATA approximates the decaying sinusoidal function
% corresponding to the experimental data given, then returns the
% coefficients.
% 
[amplitude,debut]=max(vibration);
te= time(2)-time(1);
fe=1/te;
L=length(vibration);
NFFT = 2^nextpow2(L);
Y = fft(vibration,NFFT)/L;
f = fe/2*linspace(0,1,NFFT/2+1);

[~,index_f]=max(2*abs(Y(1:NFFT/2+1)));
freq=f(index_f);
omega=2*pi*freq;

% initial estimation of coefficients
coeff0= [omega; 0; 0; amplitude];

%% Implementation of the least squares method
options = optimset('Display','off'); % suppress lsqcurvefit output
coeff= lsqcurvefit(@curve_model,coeff0,time,vibration,[], [], options);
omega=coeff(1);
damping=coeff(2);
phase=coeff(3);
amplitude=coeff(4);
end
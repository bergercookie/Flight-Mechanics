function plot_spectrum(t0, x0)

fs = 1 / (t0(2) - t0(1)); %assume sampling time constant
f_vals = linspace(0, fs, length(x0));
X0 = fft(x0);

% magnitude of fft
X0_mag = abs(X0);

figure(1); clf()
subplot(211)
hold on
plot(X0_mag);
xlabel('bins')
ylabel('Magnitude X_0');
grid on

xlim([0, length(X0_mag)]);

subplot(212);
title('Frequency Spectrum')
hold on;
plot(f_vals, X0_mag);
xlabel('frequency [Hz]')
ylabel('Magnitude X_0');
grid on

end
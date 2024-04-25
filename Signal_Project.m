clc
close all

%% Define sampling frequency and time vector
fs = 1000;
t = 0:1/fs:1;

%% Define modulation parameters
fm = 5; % frequency of modulating signal
Am = 4.5; % amplitude of modulating signal
fc = 100; % frequency of carrier signal
Ac = 20; % amplitude of carrier signal
Kf = 15; % frequency deviation constant

%% Task 1:
% Generate message signal
m = 0.5*cos(2*pi*fm*t) + cos(2*pi*fm*t) + 3*cos(2*pi*fm*t);

%% Task 2:
% Generate carrier signal
c = Ac*cos(2*pi*fc*t);


%% Task 3:
% Generate FM signal
beta = ((Kf*Am)/(fm)); % calculate peak frequency deviation from modulation parameters
fm_signal = Ac*cos(2*pi*fc*t + beta*sin(2*pi*fm*t));

%% Plot the signals
figure;
subplot(3,1,1)
plot(t, m);
xlabel('Time (s)');
ylabel('Amplitude');
title('Sinusoidal Signal');

subplot(3,1,2)
plot(t, c);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Carrier Signal');

subplot(3,1,3)
plot(t, fm_signal);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('FM Signal');


%% Task 4:
% Calculate and plot FFT of FM signal for different amplitudes
N = length(t);
f = (-N/2:N/2-1)*fs/N; % frequency vector
Am_vals = [0, 1/12, 1/3, 3]; % different modulation indices to test
figure;
for i = 1:length(Am_vals)
    beta = (Kf*Am_vals(i))/fm;
    fm_signal_i = Ac*cos(2*pi*fc*t + beta*sin(2*pi*fm*t));
    Signal_i = abs(fftshift(fft(fm_signal_i)))/N;
    subplot(2,2,i)
    plot(f,Signal_i)
    xlim([-200,200])
    title(sprintf('Beta =%.2f', beta))
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
end


%% Task 5:
% Demodulate FM signal and plot original and demodulated signals
demod_signal = demod(fm_signal,fc,fs,"fm");
figure;
subplot(2,1,1);
plot(t,m);
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');
subplot(2,1,2); 
plot(t,demod_signal);
title('Demodulated Signal');
xlabel('Time (s)');
ylabel('Amplitude');







%% Signal generation parameters
fs = 1000;         % Sample rate
t = 0:1/fs:1;      % Time vector
fm = 5;            % Message frequency
fc = 100;          % Carrier frequency
Ac = 20;           % Carrier amplitude

%% Message signal generation
m = 0.5*cos(2*pi*fm*t) + cos(2*pi*fm*t) + 3*cos(2*pi*fm*t);

%% Carrier signal generation
c = Ac*cos(2*pi*fc*t);

%% DSB-LC AM modulation
modulated = m .* c;

%% DSB-SC AM modulation
modulated_sc = ammod(m, fc, fs);

%% SSB-SC modulation
h = hilbert(m);
usb = h .* exp(1j*2*pi*fc*t);
lsb = m .* exp(-1j*2*pi*fc*t);

%% Time domain plots
figure;
subplot(5,1,1);
plot(t, m);
title('Message Signal (m)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(5,1,2);
plot(t, c);
title('Carrier Signal (c)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(5,1,3);
plot(t, modulated);
title('DSB-LC Modulated Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(5,1,4);
plot(t, modulated_sc);
title('DSB-SC Modulated Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(5,1,5);
plot(t, usb);
hold on;
plot(t, lsb);
title('SSB Modulated Signal');
xlabel('Time (s)');
ylabel('Amplitude');
legend('USB', 'LSB');

%% Frequency domain plots
f = (-fs/2:fs/length(t):fs/2-fs/length(t));
m_fft = fftshift(fft(m));
c_fft = fftshift(fft(c));
modulated_fft = fftshift(fft(modulated));
modulated_sc_fft = fftshift(fft(modulated_sc));
usb_fft = fftshift(fft(usb));
lsb_fft = fftshift(fft(lsb));

figure;
subplot(5,1,1);
plot(f, abs(m_fft));
title('Message Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot(5,1,2);
plot(f, abs(c_fft));
title('Carrier Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot(5,1,3);
plot(f, abs(modulated_fft));
title('DSB-LC Modulated Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot(5,1,4);
plot(f, abs(modulated_sc_fft));
title('DSB-SC Modulated Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot(5,1,5);
plot(f, abs(usb_fft));
hold on;
plot(f, abs(lsb_fft));
title('SSB Modulated Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
legend('USB', 'LSB');
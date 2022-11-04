%% Creation of Signal
close all; clc; clear workspace;
%Let the carrier frequency be 10 KHz
Fc = 10000;
Fs = Fc * 16;
N_bits = 1024;
baseband_dataRate = 1000;

t = 0: 1/Fs : N_bits/baseband_dataRate;
signalLen = Fs* N_bits/baseband_dataRate;
carrier_sig = 1 .* cos(2*pi*Fc*t);
% figure; plot(t, carrier_sig); ylim([-1.25 1.25]); xlim([0 0.0001]);title("Carrier Signal 1 cycle");
% figure; plot(t, carrier_sig); ylim([-1.25 1.25]); xlim([0 0.001]);title("Carrier Signal 1 bit");

%% Seeing the signal at freq domain
Y = fft(carrier_sig);
f = 1:1000:20000;
P2 = abs(Y/signalLen);
P1 = P2(1:1000:20000);

% plot(f,P1) 
% title("Single-Sided Amplitude Spectrum of carrier_sig")
% xlabel("f (Hz)")
% ylabel("Y(f)")

%% Calculating the power of the Signal

SignalPower1 = (norm(carrier_sig)^2)/signalLen;
SignalPower2 = rms(carrier_sig)^2;
SignalPower3 = bandpower(carrier_sig);

fprintf('SignalPower1: %.5f\n', SignalPower1);
fprintf('SignalPower2: %.5f\n', SignalPower2);
fprintf('SignalPower3: %.5f\n', SignalPower3);

%% Power Spectral Density

Nfft = 1024;
[Pxx,f] = pwelch(carrier_sig,gausswin(Nfft),Nfft/2,Nfft,Fs);

figure; plot(f, Pxx);
grid on
title("PSD")
xlabel("Frequency (Hz)")
ylabel("Power/Frequency (dB/Hz)")


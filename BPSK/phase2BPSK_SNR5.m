%% Phase 2: Modulation for communication
clc; close all; clear;

%Generated baseband data
N_bits = 1024;

%Let the carrier frequency be 10 KHz
Fc = 10000;

%carrier signal is 16 times oversampled
Fs = Fc * 16;

%consider the baseband data rate as 1 kbps
baseband_dataRate = 1000;
SamplePerBit = Fs / baseband_dataRate; % sampling period OR for each bit, sample 160 times.

%Modulate the data samples with carrier signal (cos(2pft))
Amp = 2;
t = 0: 1/Fs : N_bits/baseband_dataRate;
carrier_sig = Amp .* cos(2*pi*Fc*t);
% figure; plot(carrier_sig); ylim([-1.25 1.25]); xlim([0 160]); title("Carrier Sig");



signalLen = Fs* N_bits /baseband_dataRate + 1;
SNR = 5;

rng(0);
Data = randi([0 1], 1 , N_bits);
figure; stairs(Data); ylim([-5 5]); xlim([1 9]); title("Binary data from baseband");

%% Bandpass Modulation

Spower2Npower = (10.^(SNR/10));

DataStream = zeros(1, signalLen);
for i = 1: signalLen - 1
    DataStream(i) = Data(ceil(i*baseband_dataRate/Fs));
end
DataStream(signalLen) = DataStream(signalLen - 1);
DataStream = DataStream .* 2 - 1;
figure; stairs(DataStream); ylim([-5 5]); xlim([1 1440]); title("Extended to Fs");

% Modulated
Signal = carrier_sig .* DataStream;
figure; 
% stairs(DataStream); hold on; 
plot(Signal); ylim([-5 5]); xlim([1 1440]); title("Mod Signal");
% legend("Data", "Mod Signal");

%% Simulating channel and channel noise
% 
% SignalPower = bandpower(Signal);
% NoisePower_variance = SignalPower ./ Spower2Npower;
% rng(0);
% Noise = sqrt(NoisePower_variance/2) .*randn(1,signalLen);
% % figure; plot(OOK_Noise); ylim([-5 5]); xlim([1 1440]); title("Noise");
% 
% Signal_Received = Signal + Noise;
% figure; 
% % plot(OOK_Signal); hold on; 
% plot(Signal_Received); ylim([-5 5]); xlim([1 1440]); title("Signal Received");
% % legend("With noise", "Without noise");

%% At the receiver
Signal_Received = Signal;
Squared = Signal_Received .^ 2; 
figure; plot(Squared); ylim([-5 5]); xlim([1 1440]); title("Mul with 2");

Nfft = 1024;
[Pxx,f] = pwelch(Squared,gausswin(Nfft),Nfft/2,Nfft,Fs);

figure; plot(f, Pxx);
grid on
title("PSD0")
xlabel("Frequency (Hz)")
ylabel("Power/Frequency (dB/Hz)")

%% 1
[b_low, a_low] = butter(6, 0.2, 'high');
Filtered = filtfilt(b_low, a_low, Squared);
figure; plot(Filtered); ylim([-3.5 3.5]); xlim([1 1440]); title("Filtered");

Nfft = 1024;
[Pxx,f] = pwelch(Filtered,gausswin(Nfft),Nfft/2,Nfft,Fs);

figure; plot(f, Pxx);
grid on
title("PSD1")
xlabel("Frequency (Hz)")
ylabel("Power/Frequency (dB/Hz)")

%% 2
% Use the decision threshold logic for decoding of received signals
Sampled = sample(Filtered, SamplePerBit, N_bits);
figure; stairs(Sampled, "-o"); ylim([-2 2]); xlim([1 9]); title("Sampled");

%% 567
Result = zeros(1, N_bits);
for x = 1:N_bits
    if (Sampled(x) > 0)
        Result(x) = 1;
    else
        Result(x) = 0;
    end
end

figure; 
stairs(Sampled, "-o"); hold on;
yline(0); hold on;
stairs(Result); ylim([-1.25 1.25]); xlim([1 9]); title("Binary output");

figure;
stairs(Data); hold on;
stairs(Result); ylim([-5 5]); xlim([1 9]); title("Binary output");
legend("input", "Output");

function sampled = sample(x,sampling_period,num_bit)
    sampled = zeros(1, num_bit);
    for n = 1: num_bit
        sampled(n) = x((2 * n - 1) * sampling_period / 2);
    end
end

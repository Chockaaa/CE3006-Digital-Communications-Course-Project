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
Amp = 1;
t = 0: 1/Fs : N_bits/baseband_dataRate;
carrier_sig = Amp .* cos(2*pi*Fc*t);
carrier_sig2 = Amp .* cos(2*pi*(10*Fc)*t);
figure; plot(carrier_sig); title("1st freq");ylim([-1.25 1.25]); xlim([1 100]);
figure; plot(carrier_sig2); title("2nd freq");ylim([-1.25 1.25]); xlim([1 100]);

%define 6th order LP butterworth filter with 0.2 normalized cutoff frequency
[b_low,a_low] = butter(6, 0.2);
%define 6th order HP butterworth filter with 0.2 normalized cutoff frequency
[b_high,a_high] = butter(6, 0.2, 'high');

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
% DataStream = DataStream .* 2 - 1;
figure; stairs(DataStream); ylim([-5 5]); xlim([1 1440]); title("Extended to Fs");

% Modulated
Signal_BFSK1 = DataStream .* carrier_sig2;
Signal_BFSK2 = (1 - DataStream) .* carrier_sig;
figure; plot(Signal_BFSK1); title("1");ylim([-1.25 1.25]); xlim([1 1440]);
figure; plot(Signal_BFSK2); title("0");ylim([-1.25 1.25]); xlim([1 1440]);

Signal = Signal_BFSK1 + Signal_BFSK2;
figure; plot(Signal); title("10");ylim([-1.25 1.25]); xlim([1 1440]);

%% Simulating channel and channel noise

SignalPower = bandpower(Signal);
NoisePower_variance = SignalPower ./ Spower2Npower;
rng(0);
Noise = sqrt(NoisePower_variance/2) .*randn(1,signalLen);
% figure; plot(OOK_Noise); ylim([-5 5]); xlim([1 1440]); title("Noise");

Signal_Received = Signal + Noise;
figure; 
% plot(OOK_Signal); hold on; 
plot(Signal_Received); ylim([-5 5]); xlim([1 1440]); title("Signal Received");
% legend("With noise", "Without noise");

%% At the receiver

%non coherent detection -- bandpass filters
ReceiveFSK_LOW = filtfilt(b_low,a_low,Signal_Received);
figure; plot(ReceiveFSK_LOW); ylim([-3.5 3.5]); xlim([1 1440]); title("filtered low");
ReceiveFSK_HIGH = filtfilt(b_high,a_high,Signal_Received);
figure; plot(ReceiveFSK_HIGH); ylim([-3.5 3.5]); xlim([1 1440]); title("filtered high");
%Square Law
SquaredFSK_LOW = ReceiveFSK_LOW.^2;
figure; plot(SquaredFSK_LOW); ylim([-3.5 3.5]); xlim([1 1440]); title("squared low");
SquaredFSK_HIGH = ReceiveFSK_HIGH.^2;
figure; plot(SquaredFSK_HIGH); ylim([-3.5 3.5]); xlim([1 1440]); title("squared high");
SquaredFSK = SquaredFSK_HIGH - SquaredFSK_LOW;
figure; plot(SquaredFSK); ylim([-3.5 3.5]); xlim([1 1440]); title("total squared");
%sample and decision device
Sampled = sample(SquaredFSK, SamplePerBit, N_bits);
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

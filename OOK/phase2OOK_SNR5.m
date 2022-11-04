%% Phase 2: Modulation for communication
clc; close all; clear workspace;

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
% figure; plot(carrier_sig); ylim([-1.25 1.25]); xlim([0 160]); title("Carrier Sig");

%Assume a 6th order LPF with cut-off frequency 0.2 in the function
[b_low, a_low] = butter(6, 0.2);

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
% figure; stairs(DataStream); ylim([-5 5]); xlim([1 1440]); title("Extended to Fs");

% Modulated
OOK_Signal = carrier_sig .* DataStream;
figure; 
% stairs(DataStream); hold on; 
plot(OOK_Signal); ylim([-5 5]); xlim([1 1440]); title("Mod Signal");
% legend("Data", "Mod Signal");

%% Simulating channel and channel noise

OOK_SignalPower = bandpower(OOK_Signal);
OOK_NoisePower_variance = OOK_SignalPower ./ Spower2Npower;
rng(0);
OOK_Noise = sqrt(OOK_NoisePower_variance/2) .*randn(1,signalLen);
% figure; plot(OOK_Noise); ylim([-5 5]); xlim([1 1440]); title("Noise");

OOK_Signal_Received = OOK_Signal + OOK_Noise;
figure; 
% plot(OOK_Signal); hold on; 
plot(OOK_Signal_Received); ylim([-5 5]); xlim([1 1440]); title("Signal + Noise");
% legend("With noise", "Without noise");

%% At the receiver

OOK_Squared = OOK_Signal_Received.^2; %square law device (detection)
figure; plot(OOK_Squared); ylim([-0.25 3.5]); xlim([1 1440]); title("Squared");

% filtering of the demodulated signal
OOK_Filtered = filtfilt(b_low, a_low, OOK_Squared);
figure; plot(OOK_Filtered); ylim([-0.25 3.5]); xlim([1 1440]); title("Filtered");

% Use the decision threshold logic for decoding of received signals
OOK_Sampled = sample(OOK_Filtered, SamplePerBit, N_bits);
figure; stairs(OOK_Sampled, "-o"); ylim([-0.25 2]); xlim([1 9]); title("Sampled");

OOK_Result = zeros(1, N_bits);
for x = 1:N_bits
    if (OOK_Sampled(x) > ((Amp*Amp)/2))
        OOK_Result(x) = 1;
    else
        OOK_Result(x) = 0;
    end
end

figure; 
stairs(OOK_Sampled, "-o"); hold on;
yline((Amp*Amp)/2); hold on;
stairs(OOK_Result); ylim([-0.25 1.25]); xlim([1 9]); title("Binary output");

figure;
stairs(Data); hold on;
stairs(OOK_Result); ylim([-5 5]); xlim([1 9]); title("Binary output");
legend("input", "Output");

function sampled = sample(x,sampling_period,num_bit)
    sampled = zeros(1, num_bit);
    for n = 1: num_bit
        sampled(n) = x((2 * n - 1) * sampling_period / 2);
    end
end

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
Ts = Fs / baseband_dataRate; % sampling period OR for each bit, sample 160 times.

%Modulate the data samples with carrier signal (cos(2pft))
A = 1; % multiplying twice the carrier signal
t = 0: 1/Fs : N_bits/baseband_dataRate;
carrier_sig = A .* cos(2*pi*(Fc-1000)*t);
carrier_sig2 = A .* cos(2*pi*(Fc+1000)*t);

%%

%gen LPF
%Assume a 6th order filter with cut-off frequency 0.2 in the function
[b_low, a_low] = butter(6, 0.2);

signalLen = Fs* N_bits /baseband_dataRate + 1;
SNR_db_Values_Array = 0:5:50; 
ER_OOK = zeros(1, length(SNR_db_Values_Array));

rng(0);
Data = randi([0 1], 1 , N_bits);

%fill the data stream
DataStream = zeros(1, signalLen);
for i = 1: signalLen - 1
    DataStream(i) = Data(ceil(i*baseband_dataRate/Fs));
end
DataStream(signalLen) = DataStream(signalLen - 1);

DataStream_BPSK = DataStream .* 2 - 1;
Signal_BPSK = carrier_sig .* DataStream_BPSK;


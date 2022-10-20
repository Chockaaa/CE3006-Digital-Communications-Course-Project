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
carrier_sig = A .* cos(2*pi*(Fc-5000)*t);
carrier_sig2 = A .* cos(2*pi*(Fc+5000)*t);
figure; plot(carrier_sig); title("1st freq");ylim([-1.25 1.25]); xlim([1 100]);
figure; plot(carrier_sig2); title("2nd freq");ylim([-1.25 1.25]); xlim([1 100]);

%%
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
figure; stairs(DataStream); ylim([-5 5]); xlim([1 1440]);

Signal_BFSK1 = DataStream .* carrier_sig;
Signal_BFSK2 = (1 - DataStream) .* carrier_sig2;
figure; plot(Signal_BFSK1); title("1");ylim([-1.25 1.25]); xlim([1 1440]);
figure; plot(Signal_BFSK2); title("0");ylim([-1.25 1.25]); xlim([1 1440]);

Signal = Signal_BFSK1 + Signal_BFSK2;
figure; plot(Signal); title("10");ylim([-1.25 1.25]); xlim([1 1440]);




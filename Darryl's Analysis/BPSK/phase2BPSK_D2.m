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
carrier_sig = A .* cos(2*pi*Fc*t);
% figure; plot(carrier_sig); ylim([-1.25 1.25]); xlim([0 100]);
% fprintf('1\n');
% pause;

%gen LPF
%Assume a 6th order filter with cut-off frequency 0.2 in the function
[b_low, a_low] = butter(6, 0.2);

signalLen = Fs* N_bits /baseband_dataRate + 1;
SNR_db_Values_Array = 0:5:50; 
ER_BPSK = zeros(1, length(SNR_db_Values_Array));

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
figure; plot(Signal_BPSK); xlim([0 1440]); ylim([-5 5]);
fprintf('1\n');
pause;

for k = 1:length(SNR_db_Values_Array)
    SNR = (10.^(SNR_db_Values_Array(k)/10));
    
    %generate noise
    Signal_Power_BPSK = (norm(Signal_BPSK)^2)/signalLen;
    
    Noise_Power_BPSK = Signal_Power_BPSK ./SNR;
    NoiseBPSK = sqrt(Noise_Power_BPSK/2) .*randn(1,signalLen);

    %transmission
    ReceiveBPSK = Signal_BPSK+NoiseBPSK;
    %non-coherent detection -- square law
    SquaredBPSK = ReceiveBPSK .* carrier_sig;
    
    OutputBPSK = filtfilt(b_low, a_low, SquaredBPSK);

    %sample and decision device
    sampledBPSK = sample(OutputBPSK, Ts, N_bits);  

    BPSK_Result = zeros(1, N_bits);
    for x = 1:N_bits
        if (sampledBPSK(x) > 0)
            BPSK_Result(x) = 1;
        else
            BPSK_Result(x) = 0;
        end
    end

% Calculate the bit error rate performance
    BPSK_Error = 0;

    for i = 1: N_bits - 1
        if(BPSK_Result(i) ~= Data(i))
            BPSK_Error = BPSK_Error + 1;
        end
    end


    if(SNR_db_Values_Array(k) == 5)
        plot_signal = Data;
        plot_mod_BPSK = Signal_BPSK;
        plot_receive_BPSK = ReceiveBPSK;
        plot_demod_BPSK = OutputBPSK;
        plot_decoded_BPSK = BPSK_Result;
    end

    ER_BPSK(k) = (BPSK_Error/N_bits)+eps;
end

% plot the result using semilogy function
figure(1);
semilogy (SNR_db_Values_Array,ER_BPSK,'-o');
title('Error rate performance for BPSK');
ylabel('Pe');
% ylim([0.001 1]);
xlabel('Eb/No');

% plot the signals at different stages (data waveform, modulated
% signal, received signal, demodulated signal and decoded signal) 
% for a selected SNR value
figure(2);
% stairs(plot_signal);hold on;stairs(plot_decoded_OOK);title('Generated Data');ylim([-0.25 1.25]);legend('ori', 'out');
subplot(511);stairs(plot_signal);title('Generated Data');ylim([-0.25 1.25]);
subplot(512);plot(plot_mod_BPSK,'k');title('Modulated OOK');ylim([-0.25 1.25]);
subplot(513);plot(plot_receive_BPSK, 'k');title('Received Signal OOK');ylim([-0.25 1.25]);
subplot(514);plot(plot_demod_BPSK, 'k');title('Demodulated OOK');ylim([-0.25 1.25]);
subplot(515);stairs(plot_decoded_BPSK);title('Decoded Data');ylim([-0.25 1.25]);

function sampled = sample(x,sampling_period,num_bit)
    sampled = zeros(1, num_bit);
    for n = 1: num_bit
        sampled(n) = x((2 * n - 1) * sampling_period / 2);
    end
end

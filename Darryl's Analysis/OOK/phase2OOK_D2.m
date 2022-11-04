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
figure; plot(carrier_sig); ylim([-1.25 1.25]); xlim([0 100]);
fprintf('1\n');
pause;
%%
%gen LPF
%Assume a 6th order filter with cut-off frequency 0.2 in the function
[b_low, a_low] = butter(6, 0.2);

signalLen = Fs* N_bits /baseband_dataRate + 1;
SNR_db_Values_Array = 0:5:50; 
ER_OOK = zeros(1, length(SNR_db_Values_Array));

% rng(0);
% Original_Binary = randi([0 1], 1 , N_bits);

for k = 1:length(SNR_db_Values_Array)
    SNR = (10.^(SNR_db_Values_Array(k)/10));
      
    rng(0);
    Data = randi([0 1], 1 , N_bits);
    figure; plot(linspace(1, N_bits, N_bits), Data); ylim([-5 5]); xlim([1 9]);
    figure; stairs(Data); ylim([-5 5]); xlim([1 9]);
    fprintf('SNR: %d\n', SNR_db_Values_Array(k));
    fprintf('SNR2: %d\n', SNR);
    pause;

    %fill the data stream
    DataStream = zeros(1, signalLen);
    for i = 1: signalLen - 1
        DataStream(i) = Data(ceil(i*baseband_dataRate/Fs));
    end
    DataStream(signalLen) = DataStream(signalLen - 1);
    figure; plot(linspace(1, signalLen, signalLen), DataStream); ylim([-5 5]); xlim([1 1440]);
    figure; stairs(DataStream); ylim([-5 5]); xlim([1 1440]);
    fprintf('2\n');
    pause;


    % OOK
    OOK_Signal = carrier_sig .* DataStream;

    % generate noise
    OOK_SignalPower = (norm(OOK_Signal)^2)/signalLen;
    OOK_NoisePower_variance = OOK_SignalPower ./ SNR;
    OOK_Noise = sqrt(OOK_NoisePower_variance/2) .*randn(1,signalLen);

    OOK_Signal_Received = OOK_Signal + OOK_Noise;

    OOK_Squared = OOK_Signal_Received.^2; %square law device (detection)

    % filtering of the demodulated signal
    OOK_Filtered = filtfilt(b_low, a_low, OOK_Squared);

    % Use the decision threshold logic for decoding of received signals
    OOK_Sampled = sample(OOK_Filtered, Ts, N_bits);

    OOK_Result = zeros(1, N_bits);
    for x = 1:N_bits
        if (OOK_Sampled(x) > A/2)
            OOK_Result(x) = 1;
        else
            OOK_Result(x) = 0;
        end
    end

% Calculate the bit error rate performance
    OOK_Error = 0;

    for i = 1: N_bits - 1
        if(OOK_Result(i) ~= Data(i))
            OOK_Error = OOK_Error + 1;
        end
    end


    if(SNR_db_Values_Array(k) == 5)
        plot_signal = Data;
        plot_mod_OOK = OOK_Signal;
        plot_receive_OOK = OOK_Signal_Received;
        plot_demod_OOK = OOK_Filtered;
        plot_decoded_OOK = OOK_Result;
    end

    ER_OOK(k) = (OOK_Error/N_bits)+eps;
end

% plot the result using semilogy function
figure(1);
semilogy (SNR_db_Values_Array,ER_OOK,'-o');
title('Error rate performance for OOK');
ylabel('Pe');
ylim([0.001 1]);
xlabel('Eb/No');

% plot the signals at different stages (data waveform, modulated
% signal, received signal, demodulated signal and decoded signal) 
% for a selected SNR value
figure(2);
% stairs(plot_signal);hold on;stairs(plot_decoded_OOK);title('Generated Data');ylim([-0.25 1.25]);legend('ori', 'out');
subplot(511);stairs(plot_signal);title('Generated Data');ylim([-0.25 1.25]);
subplot(512);plot(plot_mod_OOK,'k');title('Modulated OOK');ylim([-0.25 1.25]);
subplot(513);plot(plot_receive_OOK, 'k');title('Received Signal OOK');ylim([-0.25 1.25]);
subplot(514);plot(plot_demod_OOK, 'k');title('Demodulated OOK');ylim([-0.25 1.25]);
subplot(515);stairs(plot_decoded_OOK);title('Decoded Data');ylim([-0.25 1.25]);

function sampled = sample(x,sampling_period,num_bit)
    sampled = zeros(1, num_bit);
    for n = 1: num_bit
        sampled(n) = x((2 * n - 1) * sampling_period / 2);
    end
end

clc; close all; clear workspace;

%Generated baseband data
N_bits = 1024 ;
Raw_Data = randi([0 1], 1 , N_bits);
Signal =  2 .* (Raw_Data - 0.5);

%Let the carrier frequency be 10 KHz
Fc = 10000;

%carrier signal is 16 times oversampled
Fs = Fc * 16;

%consider the baseband data rate as 1 kbps
baseband_dataRate = 1000;
Ts = Fs / baseband_dataRate; % sampling period

%Modulate the data samples with carrier signal (cos(2pft))
A = 2; %multiplying twice the carrier signal
t = 0: 1/Fs : N_bits/baseband_dataRate;
carrier_sig = A .* cos(2*pi*Fc*t);

%gen LPF
%Assume a 6th order filter with cut-off frequency 0.2 in the function
[b_low, a_low] = butter(6, 0.2);
[b_high, a_high] = butter(6, 0.2, "high");

signalLen = Fs*N_bits /baseband_dataRate + 1;

SNR_db_Values_Array = -20:1:20; %0:5:50;

ER_OOK = zeros(length(SNR_db_Values_Array));

for k = 1:length(SNR_db_Values_Array)
    SNR = (10.^(SNR_db_Values_Array(k)/10));

    avg_OOK_error =0;

    % generate data
    for j = 1 : 10      % each SNR avg the error over x times
       
        Data = round(rand(1,N_bits));
        
        %fill the data stream
        DataStream = zeros(1, signalLen);
        for i = 1: signalLen - 1
            DataStream(i) = Data(ceil(i*baseband_dataRate/Fs));
        end
        DataStream(signalLen) = DataStream(signalLen - 1);

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
        OOK_Result = decision_device(OOK_Sampled, N_bits, A/2);

   % Calculate the bit error rate performance
        OOK_Error = 0;

        for i = 1: N_bits - 1
            if(OOK_Result(i) ~= Data(i))
                OOK_Error = OOK_Error + 1;
            end
        end
        avg_OOK_error = OOK_Error + avg_OOK_error;

    end

    if(SNR_db_Values_Array(k) == 5)
        plot_signal = Data;
        plot_mod_OOK = OOK_Signal;
        plot_receive_OOK = OOK_Signal_Received;
        plot_demod_OOK = OOK_Filtered;
        plot_decoded_OOK = OOK_Result;
    end

    ER_OOK(k) = (avg_OOK_error / 10)/N_bits;
end

% plot the result using  semilogy’ function
figure(1);
semilogy (SNR_db_Values_Array,ER_OOK,'k-*');
title('Error rate performance for OOK');
ylabel('Pe');
xlabel('Eb/No');

% plot the signals at different stages (data waveform, modulated
% signal, received signal, demodulated signal and decoded signal) 
% for a selected SNR value
figure(2);
subplot(511);title('Generated Data');plot(plot_signal);
subplot(512);title('Modulated OOK');plot(plot_mod_OOK,'k');
subplot(513);title('Received Signal OOK');plot(plot_receive_OOK, 'k')
subplot(514);title('Demodulated OOK');plot(plot_demod_OOK, 'k');
subplot(515);title('Decoded Data');plot(plot_decoded_OOK);


function sampled = sample(x, samplingPeriod, numBit)
    sampled = zeros(1, numBit);
    for i = 1:numBit
        sampled(i) = x((2 * i - 1) * samplingPeriod / 2);
    end
end

function bin_out = decision_device(sampled, numBit, threshold)
    bin_out = zeros(1, numBit);
    for i = 1:numBit
        if(sampled(i) > threshold)
            bin_out(i) = 1;
        else 
            bin_out(i) = 0;
        end
    end
end
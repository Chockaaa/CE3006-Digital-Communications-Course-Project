%% Phase 2: Modulation for communication
clc; close all; clear;

% Generated baseband data
N_bits = 1024;

% Let the carrier frequency be 10 KHz
Fc = 10000;

% Carrier signal is 16 times oversampled
Fs = Fc * 16;

% onsider the baseband data rate as 1 kbps
baseband_dataRate = 1000;
Ts = Fs / baseband_dataRate; % sampling period

% Modulate the data samples with carrier signal (cos(2pft))
A = 1; % multiplying twice the carrier signal
t = 0: 1/Fs : N_bits/baseband_dataRate;
carrier_sig = A .* cos(2*pi*Fc*t);
No_runs = 100;

% Gen LPF
% Assume a 6th order filter with cut-off frequency 0.2 in the function
[b_low, a_low] = butter(6, 0.2);

signalLen = Fs* N_bits /baseband_dataRate + 1;
SNR_db_Values_Array = -50:5:50; %0:5:50;
ER_OOK = zeros(1, length(SNR_db_Values_Array));

for k = 1:length(SNR_db_Values_Array)
    SNR = (10.^(SNR_db_Values_Array(k)/10));   

    avg_OOK_error = 0;

    % Generate data
    for j = 1 : No_runs     % Each SNR avg the error over 10 times
        Data = randi([0 1], 1 , N_bits);
        

        % Fill the data stream
        DataStream = zeros(1, signalLen);
        
        for i = 1: signalLen - 1
            DataStream(i) = Data(ceil(i*baseband_dataRate/Fs));
            
        end
        
        DataStream(signalLen) = DataStream(signalLen - 1);

        % OOK
        OOK_Signal = carrier_sig .* DataStream;

        % Generate noise
        OOK_SignalPower = (norm(OOK_Signal)^2)/signalLen;
        
        OOK_NoisePower_variance = OOK_SignalPower ./ SNR;
        OOK_Noise = sqrt(OOK_NoisePower_variance/2) .*randn(1,signalLen);

        OOK_Signal_Received = OOK_Signal + OOK_Noise;

        OOK_Squared = OOK_Signal_Received.^2; %square law device (detection)

        % Filtering of the demodulated signal
        OOK_Filtered = filtfilt(b_low, a_low, OOK_Squared);

        % Use the decision threshold logic for decoding of received signals
        OOK_Sampled = sample(OOK_Filtered, Ts, N_bits);
        
        OOK_Result = decision_logic(OOK_Sampled,N_bits,(A*A)/2);
        
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

    ER_OOK(k) = (avg_OOK_error / No_runs)/N_bits + eps;
end

% plot the result using  semilogy’ function
figure(1);
semilogy (SNR_db_Values_Array,ER_OOK,'k-*');
title('Error rate performance for OOK');
ylabel('Pe');
xlabel('Eb/No');

% Plot the signals at different stages (data waveform, modulated
% Signal, received signal, demodulated signal and decoded signal) 
% for a selected SNR value
% figure(2);
% subplot(511);plot(plot_signal);title('Generated Data');
% subplot(512);plot(plot_mod_OOK,'k');title('Modulated OOK');
% subplot(513);plot(plot_receive_OOK, 'k');title('Received Signal OOK');
% subplot(514);plot(plot_demod_OOK, 'k');title('Demodulated OOK');
% subplot(515);plot(plot_decoded_OOK);title('Decoded Data');

function Result_Out = decision_logic(sampled,N_bits,threshold)
    Result_Out = zeros(1, N_bits);
    for x = 1:N_bits
        if (sampled(x) > threshold)
            Result_Out(x) = 1;
        else
            Result_Out(x) = 0;
        end
    end
end


function sampled = sample(x, samplingPeriod, numBit)
    sampled = zeros(1, numBit);
    for i = 1:numBit
        sampled(i) = x((2 * i - 1) * samplingPeriod / 2);
    end
end

%% Phase 2: Modulation for communication
clc; clear;
% close all;

% Generated baseband data
N_bits = 1024;

% Let the carrier frequency be 10 KHz
Fc = 10000;

% Carrier signal is 16 times oversampled
Fs = Fc * 16;

% onsider the baseband data rate as 1 kbps
baseband_dataRate = 1000;
SamplesPerBit = Fs / baseband_dataRate; % sampling period

% Modulate the data samples with carrier signal (cos(2pft))
A = 1; 
t = 0: 1/Fs : N_bits/baseband_dataRate;
carrier_sig = A .* cos(2*pi*Fc*t);
No_runs = 20;

% Gen LPF
% Assume a 6th order filter with cut-off frequency 0.2 in the function
[b_low, a_low] = butter(6, 0.2);

signalLen = Fs* N_bits /baseband_dataRate + 1;
SNR_db_Values_Array = -50:5:50; %0:5:50;
Bit_Error_Rate = zeros(1, length(SNR_db_Values_Array));

for k = 1:length(SNR_db_Values_Array)
    Spower_2_Npower = (10.^(SNR_db_Values_Array(k)/10));   

    avg_error = 0;

    % Generate data
    for j = 1 : No_runs     % Each SNR avg the error over 100 times
        Data = randi([0 1], 1 , N_bits);
        
        % Fill the data stream
        DataStream = zeros(1, signalLen);
        
        for i = 1: signalLen - 1
            DataStream(i) = Data(ceil(i*baseband_dataRate/Fs));
            
        end
        
        DataStream(signalLen) = DataStream(signalLen - 1);

        Signal = carrier_sig .* DataStream;

        % Generate noise
        SignalPower = (norm(Signal)^2)/signalLen;
        
        NoisePower_variance = SignalPower ./ Spower_2_Npower;
        Noise = sqrt(NoisePower_variance/2) .*randn(1,signalLen);

        Signal_Received = Signal + Noise;

        Squared = Signal_Received.^2; %square law device (detection)

        % Filtering of the demodulated signal
        Filtered = filtfilt(b_low, a_low, Squared);

        % Use the decision threshold logic for decoding of received signals
        Sampled = sample(Filtered, SamplesPerBit, N_bits);
        
        Result = decision_logic(Sampled,N_bits,(A*A)/2);
        
        % Calculate the bit error rate performance
        No_Bit_Error = 0;    
    
        for i = 1: N_bits - 1
            if(Result(i) ~= Data(i))
                No_Bit_Error = No_Bit_Error + 1;
            end
        end     
        avg_error = No_Bit_Error + avg_error;
    
    end

    if(SNR_db_Values_Array(k) == 5)
        plot_DS = DataStream;
        plot_signal = Data;
        plot_mod = Signal;
        plot_receive = Signal_Received;
        plot_demod = Filtered;
        plot_decoded = Result;
    end

    Bit_Error_Rate(k) = (avg_error / No_runs)/N_bits;
end

% plot the result using  semilogy’ function
figure;
semilogy (SNR_db_Values_Array,Bit_Error_Rate,'k-*');
title('Error rate performance for OOK');
ylabel('Pe');
xlabel('Eb/No');

% Plot the signals at different stages (data waveform, modulated
% Signal, received signal, demodulated signal and decoded signal) 
% for a selected SNR value

bits_y_range = [-0.25 1.25];
cont_y_range = [-2 2];

% figure(2);
% subplot(511); stairs(plot_signal);title('Generated Data');ylim(bits_y_range); xlim([1 9]);
% subplot(512); plot(plot_mod);title('Modulated');ylim(cont_y_range); xlim([1 1280]);xticks(0:160:1280);
% subplot(513); plot(plot_receive);title('Received Signal');ylim(cont_y_range); xlim([1 1280]);xticks(0:160:1280);
% subplot(514); plot(plot_demod);title('Demodulated');ylim(bits_y_range); xlim([1 1280]);xticks(0:160:1280);
% subplot(515); stairs(plot_decoded);title('Decoded Data');ylim(bits_y_range); xlim([1 9]);

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

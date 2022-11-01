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
SamplesPerBit = Fs / baseband_dataRate; % sampling period

% Modulate the data samples with carrier signal (cos(2pft))
Amp = 1; % multiplying twice the carrier signal
t = 0: 1/Fs : N_bits/baseband_dataRate;
carrier_sig = Amp .* cos(2*pi*Fc*t);
carrier_sig2 = Amp .* cos(2*pi*(10*Fc)*t);
No_runs = 100;

% Assume a 6th order filter with cut-off frequency 0.2 in the function
[b_low, a_low] = butter(6, 0.2);
%define 6th order HP butterworth filter with 0.2 normalized cutoff frequency
[b_high,a_high] = butter(6, 0.2, 'high');

signalLen = Fs* N_bits /baseband_dataRate + 1;
SNR_db_Values_Array = -50:5:50; %0:5:50;
Bit_Error_Rate = zeros(1, length(SNR_db_Values_Array));

for k = 1:length(SNR_db_Values_Array)
    Spower_2_Npower = (10.^(SNR_db_Values_Array(k)/10));   

    avg_error = 0;

    % Generate data
    for j = 1 : No_runs     % Each SNR avg the error over 10 times
        Data = randi([0 1], 1 , N_bits);
        
        % Fill the data stream
        DataStream = zeros(1, signalLen);
        
        for i = 1: signalLen - 1
            DataStream(i) = Data(ceil(i*baseband_dataRate/Fs));
            
        end
        
        DataStream(signalLen) = DataStream(signalLen - 1);
        
        Signal_1 = DataStream .* carrier_sig2;
        Signal_0 = (1 - DataStream) .* carrier_sig;
        
        Signal = Signal_0 + Signal_1;

        % Generate noise
        SignalPower = bandpower(Signal);
        
        NoisePower_variance = SignalPower ./ Spower_2_Npower;
        Noise = sqrt(NoisePower_variance/2) .*randn(1,signalLen);

        Signal_Received = Signal + Noise;
        
        % Fiter the signal received first.
        Filtered_0 = filtfilt(b_low,a_low,Signal_Received);
        Filtered_1 = filtfilt(b_high,a_high,Signal_Received);

        % The squaring of output;
        Squared_0 = Filtered_0.^2;
        Squared_1 = Filtered_1.^2;
        
        % Adding both;
        Summation = Squared_1 - Squared_0;
        
        Sampled = sample(Summation, SamplesPerBit, N_bits);
        
        Result = decision_logic(Sampled,N_bits,0);
        
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
        plot_signal = Data;
        plot_mod = Signal;
        plot_receive = Signal_Received;
        plot_demod = Summation;
        plot_decoded = Result;
    end

    Bit_Error_Rate(k) = (avg_error / No_runs)/N_bits + eps;
end

% plot the result using  semilogy’ function
figure(1);
semilogy (SNR_db_Values_Array,Bit_Error_Rate,'k-*');
title('Error rate performance for BFSK');
ylabel('Pe');
xlabel('Eb/No');

% Plot the signals at different stages (data waveform, modulated
% Signal, received signal, demodulated signal and decoded signal) 
% for a selected SNR value
figure(2);
subplot(511);plot(plot_signal);title('Generated Data');
subplot(512);plot(plot_mod,'k');title('Modulated');
subplot(513);plot(plot_receive, 'k');title('Received Signal');
subplot(514);plot(plot_demod, 'k');title('Demodulated');
subplot(515);plot(plot_decoded);title('Decoded Data');

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

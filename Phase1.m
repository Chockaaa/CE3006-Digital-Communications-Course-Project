%% Phase 1: Data Generation
clc; close all; clear workspace;

% Assume the number of bits for transmission is 1024 (means, N= 1024).
% Generate random binary digits (0 or 1)
% Convert binary digits to ±1 (means 1 to +1 and 0 to -1). This is your data for transmission.

N_bits = 1024 ;
Raw_Data = randi([0 1], 1 , N_bits);
Signal =  2 .* (Raw_Data - 0.5);

% Generate equal number of noise samples.
% The generated noise should have normal distribution with zero mean and unit variance (use the function randn in MATLAB).

Noise = randn(1, N_bits);

% Change the noise variance with respect to SNR (signal to noise ratio) value.
% For that, fix the SNR value, (For example, let SNR = 10 dB)
% Use SNR value to generate noise variance.

SNR_db_Values_Array = 0:5:50;

Result = zeros([1 11]);

for k = 1:length(SNR_db_Values_Array)
    % SNR (in dB) = 10log10 (S/N) where S is the Signal power (or variance) and N is the Noise power (or variance)
    % Assume signal (the input data) has unit power. 
    % That is, S=1.Obtain the noise variance (=N) from the previous relation and use it together with the noise samples to generate the required noise.
    SNR = SNR_db_Values_Array(k);
    SignalPower = 1;
    NoisePower_variance = SignalPower ./ (10 .^ (SNR/10));

    Noise = sqrt(NoisePower_variance) .* Noise;

    %  Add noise samples with transmitted data. This is assumed as your received signal.
    Signal_Received = Signal + Noise;

    % Consider a threshold logic at the receiver.
    % Fix the threshold value as 0 (the transmitted data is +1 and -1, and 0 is the mid value)
    Threshold = 0;
    
    % If the received signal is above or equal to the threshold level, take it as 1
    % If the received signal is below the threshold value, take it as 0.
    
    Output_Signal = zeros(1,N_bits);
    Error_Count = 0;


    for i = 1:N_bits
        if (Signal_Received(i) > Threshold)
            Output_Signal(i) = 1;
        else
            Output_Signal(i) = 0;
        end

        % Compute the bit error rate during transmission.
        % Compare the output values from the threshold logice with the input binary digits (1 or 0 format)
        % Compute the bit error rate using the relation Bit error rate = (number of errors during transmission)/(total number of bits for transmission)
        if Output_Signal(i) ~= Raw_Data(i) 
            Error_Count = Error_Count + 1;
        end
    end
    
    Bit_Error_Rate = Error_Count ./ N_bits;
    Result(k) = Bit_Error_Rate;

% Repeat the steps for different SNRs
% Repeat the steps from noise addition to bit error rate computation
% Do the simulation for the same data set with different SNRs
% Consider different SNR values from 0 dB to 50 dB (in multiples of 5 dB).
end


% Plot your result in a graph.
% X axis should be SNR values and Y axis should be bit error rate.

figure(1)
plot(SNR_db_Values_Array,Result)
xlabel('SNR Values (dB)');
ylabel('Bit Error Rate (BER)');
title("BER vs SNR");


%% Phase 2: Modulation for Communication
clc; close all; clear workspace;




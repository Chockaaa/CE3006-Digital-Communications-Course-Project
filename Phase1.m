%% Phase 1: Data Generation
clc; close all; clear workspace;

% Assume the number of bits for transmission is 1024 (means, N= 1024).
% Generate random binary digits (0 or 1)
% Convert binary digits to ±1 (means 1 to +1 and 0 to -1). This is your data for transmission.

rng(0);
N_bits = 1024;
Original_Binary = randi([0 1], 1 , N_bits);
Original_Signal =  2 .* (Original_Binary - 0.5);
% figure; stairs(Original_Binary, '-o'); hold on; stairs(Original_Signal, '-o'); 
% ylim([-2 2]); xlim([1 9]); legend('Binary', 'Signal');


% fprintf('Program paused. Press enter to continue.\n');
% pause;

% Generate equal number of noise samples.
% The generated noise should have normal distribution with zero mean and unit variance (use the function randn in MATLAB).
rng(0);
Noise = randn(1, N_bits);
% figure; plot(linspace(1, N_bits, N_bits), Noise); ylim([-5 5]); xlim([1 9]);
% figure; stairs(Noise); ylim([-5 5]); xlim([1 9]);

% Change the noise variance with respect to SNR (signal to noise ratio) value.
% For that, fix the SNR value, (For example, let SNR = 10 dB)
% Use SNR value to generate noise variance.

SNR_db_Values_Array = 0:5:50;

Result = zeros([1 11]);

% fprintf('Program paused. Press enter to continue.\n');
% pause;

Noise_bef = Noise;
for k = 1:length(SNR_db_Values_Array)
    % SNR (in dB) = 10log10 (S/N) where S is the Signal power (or variance) and N is the Noise power (or variance)
    % Assume signal (the input data) has unit power. 
    % That is, S=1.Obtain the noise variance (=N) from the previous relation and use it together with the noise samples to generate the required noise.
    SNR = SNR_db_Values_Array(k);
    
    fprintf('SNR: %d\n', SNR);
    
    SignalPower = 1;
    NoisePower = SignalPower ./ (10 .^ (SNR/10));
    
    Noise = sqrt(NoisePower) .* Noise;
    
%     figure; stairs(Noise_bef, '-o'); hold on; stairs(Noise, '-o'); 
%     xlim([1 9]); legend('Noise_bef', 'Noise');

    %  Add noise samples with transmitted data. This is assumed as your received signal.
    Signal_Received = Original_Signal + Noise;

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
        if Output_Signal(i) ~= Original_Binary(i) %returns True when output and raw data are not equal.
            Error_Count = Error_Count + 1;
        end
    end
    
    Bit_Error_Rate = Error_Count ./ N_bits;
    Result(k) = Bit_Error_Rate;
    
%     fprintf('Program paused. Press enter to continue.\n');
%     pause;
    
%     figure;
%     stairs(Raw_Data, '-o');
%     hold on
% %     stairs(Signal_Received, '-o');
% %     hold on
%     stairs(Output_Signal, '-o');
%     ylim([-10 10]); 
% %     xlim([1 9]);
%     title("SNR: " + SNR);
%     legend('Original', 'Output');

% Repeat the steps for different SNRs
% Repeat the steps from noise addition to bit error rate computation
% Do the simulation for the same data set with different SNRs
% Consider different SNR values from 0 dB to 50 dB (in multiples of 5 dB).
end


% Plot your result in a graph.
% X axis should be SNR values and Y axis should be bit error rate.

figure;
plot(SNR_db_Values_Array,Result)
ylim([-0.1 0.2]); xlim([-1 51]);
xlabel('SNR Values (dB)');
ylabel('Bit Error Rate (BER)');
title("BER vs SNR");

% ASSUMPTIONS:
fprintf('Program paused. Press enter to continue.\n');
pause;

%% Clear everything
clc; close all; clear workspace;




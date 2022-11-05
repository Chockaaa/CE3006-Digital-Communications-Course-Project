%--Admin stuff--%
clear; close all; clc;

%define carrier frequency
fc = 10000; %10kHz
%16 times oversampled -> sample freq = 16 fc
fs = 16 * fc;

%define data rate of 1kbps
dataRate = 1000;
%define number of data bits
nBits = 1024;
CC_nBits = nBits/4*7;   %we doing (7,4) code
%define sampling rate
samplingPeriod = fs / dataRate;

%define Amplitude
Amp = 1;
%define time steps
t = 0: 1/fs : CC_nBits/dataRate;
t_pure = 0:1/fs : nBits/dataRate;    %pure -- no encode version

%generate carrier frequency
Carrier = Amp .* cos(2*pi*fc*t);
%enconder carrier frequency
Carrier_pure = Amp.* cos(2*pi*fc*t_pure);

%calculate signal length
SignalLength = fs*CC_nBits/dataRate + 1;
SignalLength_pure = fs*nBits/dataRate +1;

%SNR_dB = 10 log (Signal_Power/Noise_Power)                 
SNR_dB = -50:5:50;
%==> SNR = Signal_Power/Noise_Power = 10^(SNR_dB/10)
SNR = (10.^(SNR_dB/10));

%set run times
Total_Run = 20;

%define placeholder for error calculation
% Error_Rate_Hamming = zeros(length(SNR));
Error_Rate_NoEncode = zeros(length(SNR));

%for each SNR value
for i = 1 : length(SNR)
% 	Avg_Error_Hamming = 0;
    Avg_Error_NoEncode = 0;
    
    %for each SNR value, average the error over %Total_Run times
	for j = 1 : Total_Run
        
        %-----Data generation-----%
        Data = round(rand(1,nBits));
%         EncodeHamming = encode(Data, 7, 4, 'hamming/fmt'); 

        %fill the data stream
%         DataStream_Hamming = zeros(1, SignalLength);

%         for k = 1: SignalLength - 1
%             DataStream_Hamming(k) = EncodeHamming(ceil(k*dataRate/fs));
%         end
%         DataStream_Hamming(SignalLength) = DataStream_Hamming(SignalLength - 1);

        DataStream_NoEncode = zeros(1, SignalLength_pure);
        for k = 1:SignalLength_pure -1
            DataStream_NoEncode(k)= Data(ceil(k*dataRate/fs));
        end
        DataStream_NoEncode(SignalLength_pure) = DataStream_NoEncode(SignalLength_pure-1);
        
        
        %----- OOK -----%
%         resultOOK_Hamming = OOK_transmission(DataStream_Hamming,SNR(i),Carrier,SignalLength,samplingPeriod,CC_nBits,Amp);
        resultOOK_NoEncode = OOK_transmission(DataStream_NoEncode,SNR(i),Carrier_pure,SignalLength_pure,samplingPeriod,nBits,Amp);
        
%         decodedHamming = decode(resultOOK_Hamming,7,4,'hamming/fmt');
        
        %--Calculate Error--%
%         ErrorHamming = 0;
        ErrorNoEncode = 0;
        for k = 1: nBits
%             if(decodedHamming(k) ~= Data(k))
%                 ErrorHamming = ErrorHamming + 1;
%             end
            if (resultOOK_NoEncode(k) ~= Data(k))
                ErrorNoEncode = ErrorNoEncode + 1;
            end
        end
%         Avg_Error_Hamming = ErrorHamming + Avg_Error_Hamming;
        Avg_Error_NoEncode = ErrorNoEncode + Avg_Error_NoEncode;
    end
   
%     Error_Rate_Hamming(i) = Avg_Error_Hamming/Total_Run/nBits;
    Error_Rate_NoEncode(i) = Avg_Error_NoEncode/Total_Run/nBits;
end


%Error plot
figure;
% s1 = semilogy (SNR_dB, Error_Rate_Hamming,'r-*');
% hold on
s2 = semilogy(SNR_dB, Error_Rate_NoEncode, 'k-*');
hold on
title('Error rate');
% legend([s1(1),s2(1)],'hamming','None');
ylabel('Pe');
xlabel('Eb/No')




%%--HELPER FUNCTION--%%
function sampled = sample(x,sampling_period,num_bit)
    sampled = zeros(1, num_bit);
    for n = 1: num_bit
        sampled(n) = x((2 * n - 1) * sampling_period / 2);
    end
end


%This function simulates the decision device
function binary_out = decision_logic(sampled,num_bit,threshold)
    binary_out = zeros(1,num_bit);
    for n = 1:num_bit
        if(sampled(n) > threshold)
            binary_out(n) = 1;
        else 
            binary_out(n) = 0;
        end
    end
end

%This function is a wrapper for Phase 2 OOK 
function Result = OOK_transmission(DataStream,Spower_2_Npower,carrier_sig,signalLen,SamplesPerBit,N_bits,A)
        [b_low, a_low] = butter(6, 0.2);        
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

end
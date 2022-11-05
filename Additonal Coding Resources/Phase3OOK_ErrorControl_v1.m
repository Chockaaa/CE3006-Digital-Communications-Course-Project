%% Phase 3: Basic Error Control Coding to improve the performacne
%% Cyclic Block Code
clc; close all; clear workspace;

n = 7; % Codeword length
k = 4; % Message length

% Create a generator polynomial for a cyclic code. 
% Create a parity-check matrix by using the generator polynomial.
% Create a syndrome decoding table by using the parity-check matrix.
genpoly = cyclpoly(n,k);
parmat = cyclgen(n,genpoly);
trt = syndtable(parmat);


% Generated baseband data
N_bits = 1024;
Enc_N_bits = N_bits/k*n; %need to find values for cycle block code


% Let the carrier frequency be 10 KHz
Fc = 10000;
% Carrier signal is 16 times oversampled
Fs = Fc * 16;
% Consider the baseband data rate as 1 kbps
baseband_dataRate = 1000;
Ts = Fs / baseband_dataRate; % sampling period

% Modulate the data samples with carrier signal (cos(2pft))
A = 1;
t = 0: 1/Fs : N_bits/baseband_dataRate;
Enc_t = 0: 1/Fs : Enc_N_bits/baseband_dataRate;

carrier_sig = A .* cos(2*pi*Fc*t);
Enc_carrier_sig = A .* cos(2*pi*Fc*Enc_t);

% Gen LPF
% Assume a 6th order filter with cut-off frequency 0.2 in the function
[b_low, a_low] = butter(6, 0.2);
SNR_db_Values_Array = -50:5:50; %0:5:50;

signalLen = Fs* N_bits /baseband_dataRate + 1;
Enc_signalLen = Fs* Enc_N_bits /baseband_dataRate + 1;



ER_OOK = zeros(length(SNR_db_Values_Array));
CBC_ER_OOK = zeros(length(SNR_db_Values_Array));
Hamming_ER_OOK = zeros(length(SNR_db_Values_Array));

for x = 1:length(SNR_db_Values_Array)
    SNR = (10.^(SNR_db_Values_Array(x)/10));   

    avg_OOK_error = 0;
    CBC_avg_OOK_error = 0;
    Hamming_avg_OOK_error = 0;
    % Generate data
    for j = 1 : 20     % Each SNR avg the error over 10 times
        Data = randi([0 1], 1 , N_bits);
        

        % OOK: Data stream
        DataStream = zeros(1, signalLen);      
        
        for i = 1: signalLen - 1
            DataStream(i) = Data(ceil(i*baseband_dataRate/Fs));            
        end          
        DataStream(signalLen) = DataStream(signalLen - 1);
        OOK_Signal = carrier_sig .* DataStream;

        % OOK Cyclic Block Code: Data stream 
        CBCSignal = encode(Data,n,k,'cyclic/binary',genpoly);
        CBC_DataStream = zeros(1,Enc_signalLen);

        for i = 1:Enc_signalLen -1
            CBC_DataStream(i) = CBCSignal(ceil(i*baseband_dataRate/Fs));
        end
        CBC_DataStream(Enc_signalLen) = CBC_DataStream(Enc_signalLen -1);
        CBC_OOK_Signal = Enc_carrier_sig .* CBC_DataStream;



        % OOK Hamming: Data stream         
        HammingSignal = encode(Data,n,k,'hamming/binary');
        Hamming_DataStream = zeros(1,Enc_signalLen);
        for i = 1:Enc_signalLen -1
            Hamming_DataStream(i) = HammingSignal(ceil(i*baseband_dataRate/Fs));
        end
        Hamming_DataStream(Enc_signalLen) = Hamming_DataStream(Enc_signalLen -1);
        Hamming_OOK_Signal = Enc_carrier_sig .* Hamming_DataStream;

       
        % Generate noise OOK
        OOK_SignalPower = (norm(OOK_Signal)^2)/signalLen;        
        OOK_NoisePower_variance = OOK_SignalPower ./ SNR;
        OOK_Noise = sqrt(OOK_NoisePower_variance/2) .*randn(1,signalLen);
        % Generate noise Cyclic Block Code OOK
        CBC_OOK_SignalPower = (norm(CBC_OOK_Signal)^2)/Enc_signalLen;        
        CBC_OOK_NoisePower_variance = CBC_OOK_SignalPower ./ SNR;
        CBC_OOK_Noise = sqrt(CBC_OOK_NoisePower_variance/2) .*randn(1,Enc_signalLen);
        % Generate noise Hamming OOK
        Hamming_OOK_SignalPower = (norm(Hamming_OOK_Signal)^2)/Enc_signalLen;        
        Hamming_OOK_NoisePower_variance = Hamming_OOK_SignalPower ./ SNR;
        Hamming_OOK_Noise = sqrt(Hamming_OOK_NoisePower_variance/2) .*randn(1,Enc_signalLen);

        % Transmit Signal OOK
        OOK_Signal_Received = OOK_Signal + OOK_Noise;     
        % Transmit Signal Cyclic Block Code: OOK
        CBC_OOK_Signal_Received = CBC_OOK_Signal + CBC_OOK_Noise;
        % Transmit Signal Hamming: OOK
        Hamming_OOK_Signal_Received = Hamming_OOK_Signal + Hamming_OOK_Noise;
        
        
        %square law device (detection) OOK
        OOK_Squared = OOK_Signal_Received.^2; %square law device (detection)
        %square law device (detection) Cyclic Block Code: OOK
        CBC_OOK_Squared = CBC_OOK_Signal_Received.^2; 
        %square law device (detection) Hamming: OOK
        Hamming_OOK_Squared = Hamming_OOK_Signal_Received.^2; 


        % Filtering of the demodulated signal OOK
        OOK_Filtered = filtfilt(b_low, a_low, OOK_Squared);
        % Filtering of the demodulated signal Cyclic Block Code: OOK
        CBC_OOK_Filtered = filtfilt(b_low, a_low, CBC_OOK_Squared);
        % Filtering of the demodulated signal Cyclic Block Code: OOK
        Hamming_OOK_Filtered = filtfilt(b_low, a_low, Hamming_OOK_Squared);


        % Use the decision threshold logic for decoding of received signals
        % OOK
        OOK_Sampled = sample(OOK_Filtered, Ts, N_bits);        
        OOK_Result = decision_logic(OOK_Sampled,N_bits,A*A/2);

        % Use the decision threshold logic for decoding of received signals
        % Cyclic Block Code:
        CBC_OOK_Sampled = sample(CBC_OOK_Filtered, Ts, Enc_N_bits);        
        CBC_OOK_Result = decision_logic(CBC_OOK_Sampled,Enc_N_bits,(A*A)/2);
        CBC_OOK_DecodedResult = decode(CBC_OOK_Result,n,k,'cyclic/binary',genpoly,trt);

        % Use the decision threshold logic for decoding of received signals
        % Cyclic Block Code:
        Hamming_OOK_Sampled = sample(Hamming_OOK_Filtered, Ts, Enc_N_bits);        
        Hamming_OOK_Result = decision_logic(Hamming_OOK_Sampled,Enc_N_bits,(A*A)/2);
        Hamming_OOK_DecodedResult = decode(Hamming_OOK_Result,n,k,'hamming/binary');
        
        % Calculate the bit error rate performance OOK
        OOK_Error = 0;        
        for i = 1: N_bits - 1
            if(OOK_Result(i) ~= Data(i))
                OOK_Error = OOK_Error + 1;
            end
        end     
        avg_OOK_error = OOK_Error + avg_OOK_error;

        % Calculate the bit error rate performance Cyclic Block Code OOK
        CBC_OOK_Error = 0;        
        for i = 1: N_bits - 1
            if(CBC_OOK_DecodedResult(i) ~= Data(i))
                CBC_OOK_Error = CBC_OOK_Error + 1;
            end
        end     
        CBC_avg_OOK_error = CBC_OOK_Error + CBC_avg_OOK_error;

        % Calculate the bit error rate performance Hamming OOK
        Hamming_OOK_Error = 0;        
        for i = 1: N_bits - 1
            if(Hamming_OOK_DecodedResult(i) ~= Data(i))
                Hamming_OOK_Error = Hamming_OOK_Error + 1;
            end
        end     
        Hamming_avg_OOK_error = Hamming_OOK_Error + Hamming_avg_OOK_error;
    
    end

    if(SNR_db_Values_Array(x) == 5)
        % OOK
        plot_signal = Data;
        plot_mod_OOK = OOK_Signal;
        plot_receive_OOK = OOK_Signal_Received;
        plot_demod_OOK = OOK_Filtered;
        plot_decoded_OOK = OOK_Result;

        % OOK Cyclic Block Coding
        CBC_plot_signal = Data;
        CBC_Encoded_signal = CBCSignal;
        CBC_plot_mod_OOK = CBC_OOK_Signal;
        CBC_plot_receive_OOK = CBC_OOK_Signal_Received;
        CBC_plot_demod_OOK = CBC_OOK_Filtered;
        CBC_plot_decoded_OOK = CBC_OOK_DecodedResult;

        % OOK Hamming
        Hamming_plot_signal = Data;
        Hamming_Encoded_signal = HammingSignal;
        Hamming_plot_mod_OOK = Hamming_OOK_Signal;
        Hamming_plot_receive_OOK = Hamming_OOK_Signal_Received;
        Hamming_plot_demod_OOK = Hamming_OOK_Filtered;
        Hamming_plot_decoded_OOK = Hamming_OOK_DecodedResult;

    end

    ER_OOK(x) = (avg_OOK_error / 20)/N_bits;
    CBC_ER_OOK(x) = (CBC_avg_OOK_error / 20)/N_bits;
    Hamming_ER_OOK(x) = (Hamming_avg_OOK_error / 20)/N_bits;
end

% plot the result using  semilogy’ function
figure(1);
s1 = semilogy (SNR_db_Values_Array,ER_OOK,'k-*');
hold on;
s2 = semilogy (SNR_db_Values_Array,CBC_ER_OOK,'r-o');
hold on;
s3 = semilogy (SNR_db_Values_Array,Hamming_ER_OOK,'b-+');
hold off
hold on

title('Error rate performance for OOK, Cyclic Block Code, Hamming Code');
ylabel('Pe');
xlabel('Eb/No');
legend([s1(1),s2(1),s3(1)],'OOK','OOK-CBC','OOK-Hamming','location','northeast');
hold off

% Plot the signals at different stages (data waveform, modulated
% Signal, received signal, demodulated signal and decoded signal) 
% for a selected SNR value
figure(2);
subplot(511);plot(plot_signal);title('OOK: Generated Data');
subplot(512);plot(plot_mod_OOK,'k');title('OOK: Modulated OOK');
subplot(513);plot(plot_receive_OOK, 'k');title('OOK: Received Signal OOK');
subplot(514);plot(plot_demod_OOK, 'k');title('OOK: Demodulated OOK');
subplot(515);plot(plot_decoded_OOK);title('OOK: Decoded Data');

figure(3);
subplot(611);plot(CBC_plot_signal);title('Cyclic Block Coding OOK: Generated Data');
subplot(612);plot(CBC_Encoded_signal,'k');title('Cyclic Block Coding OOK: Encoded Data');
subplot(613);plot(CBC_plot_mod_OOK,'k');title('Cyclic Block Coding OOK: Modulated OOK');
subplot(614);plot(CBC_plot_receive_OOK, 'k');title('Cyclic Block Coding OOK: Received Signal OOK');
subplot(615);plot(CBC_plot_demod_OOK, 'k');title('Cyclic Block Coding OOK: Demodulated OOK');
subplot(616);plot(CBC_OOK_DecodedResult);title('Cyclic Block Coding OOK: Decoded Data');

figure(4);
subplot(611);plot(Hamming_plot_signal);title('Hamming Code OOK: Generated Data');
subplot(612);plot(Hamming_Encoded_signal,'k');title('Hamming Code OOK: Encoded Data');
subplot(613);plot(Hamming_plot_mod_OOK,'k');title('Hamming Code OOK: Modulated OOK');
subplot(614);plot(Hamming_plot_receive_OOK, 'k');title('Hamming Code OOK: Received Signal OOK');
subplot(615);plot(Hamming_plot_demod_OOK, 'k');title('Hamming Code OOK: Demodulated OOK');
subplot(616);plot(Hamming_OOK_DecodedResult);title('Hamming Code OOK: Decoded Data');



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

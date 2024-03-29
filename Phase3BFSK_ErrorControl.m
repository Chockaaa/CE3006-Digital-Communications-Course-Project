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

gen_mat = gen2par(parmat);

% Generated baseband data
N_bits = 1024;
Enc_N_bits = N_bits/k * n; %need to find values for cycle block code


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
carrier_sig2 = A .* cos(2*pi*(10*Fc)*t);
Enc_carrier_sig = A .* cos(2*pi*Fc*Enc_t);
Enc_carrier_sig2 = A .* cos(2*pi*(10*Fc)*Enc_t);

% Gen LPF
% Assume a 6th order filter with cut-off frequency 0.2 in the function
[b_low, a_low] = butter(6, 0.2);
%define 6th order HP butterworth filter with 0.2 normalized cutoff frequency
[b_high,a_high] = butter(6, 0.2, 'high');

SNR_db_Values_Array = -50:5:50; %0:5:50;

signalLen = Fs* N_bits /baseband_dataRate + 1;
Enc_signalLen = Fs* Enc_N_bits /baseband_dataRate + 1;



ER_BFSK = zeros(1,length(SNR_db_Values_Array));
CBC_ER_BFSK = zeros(1,length(SNR_db_Values_Array));
LBC_ER_BFSK = zeros(1,length(SNR_db_Values_Array));
Hamming_ER_BFSK = zeros(1,length(SNR_db_Values_Array));

num_of_runs = 10;
for x = 1:length(SNR_db_Values_Array)
    SNR = (10.^(SNR_db_Values_Array(x)/10));   

    avg_BFSK_error = 0;
    CBC_avg_BFSK_error = 0;
    LBC_avg_BFSK_error = 0;
    Hamming_avg_BFSK_error = 0;

    % Generate data
    for j = 1 : num_of_runs     % Each SNR avg the error over 10 times
        Data = randi([0 1], 1 , N_bits);
        

        % BFSK: Data stream
        DataStream = zeros(1, signalLen);      
        
        for i = 1: signalLen - 1
            DataStream(i) = Data(ceil(i*baseband_dataRate/Fs));            
        end          
        DataStream(signalLen) = DataStream(signalLen - 1);
        
        BFSK_Signal_1 = DataStream .* carrier_sig2;
        BFSK_Signal_0 = (1 - DataStream) .* carrier_sig;
        
        BFSK_Signal = BFSK_Signal_0 + BFSK_Signal_1;

        % BFSK Cyclic Block Code: Data stream 
        CBCSignal = encode(Data,n,k,'cyclic/binary',genpoly);
        CBC_DataStream = zeros(1,Enc_signalLen);

        for i = 1:Enc_signalLen -1
            CBC_DataStream(i) = CBCSignal(ceil(i*baseband_dataRate/Fs));
        end
        CBC_DataStream(Enc_signalLen) = CBC_DataStream(Enc_signalLen -1);
        
        CBC_BFSK_Signal_1 = CBC_DataStream .* Enc_carrier_sig2;
        CBC_BFSK_Signal_0 = (1 - CBC_DataStream) .* Enc_carrier_sig;
        
        CBC_BFSK_Signal = CBC_BFSK_Signal_0 + CBC_BFSK_Signal_1;

        % BFSK Linear Block Code: Data stream 
        LBCSignal = encode(Data,n,k,'linear/binary',gen_mat);
        LBC_DataStream = zeros(1,Enc_signalLen);

        for i = 1:Enc_signalLen -1
            LBC_DataStream(i) = LBCSignal(ceil(i*baseband_dataRate/Fs));
        end
        LBC_DataStream(Enc_signalLen) = LBC_DataStream(Enc_signalLen -1);
        
        LBC_BFSK_Signal_1 = LBC_DataStream .* Enc_carrier_sig2;
        LBC_BFSK_Signal_0 = (1 - LBC_DataStream) .* Enc_carrier_sig;
        
        LBC_BFSK_Signal = LBC_BFSK_Signal_0 + LBC_BFSK_Signal_1;

        % BFSK Hamming: Data stream         
        HammingSignal = encode(Data,n,k,'hamming/binary');
        Hamming_DataStream = zeros(1,Enc_signalLen);
        for i = 1:Enc_signalLen -1
            Hamming_DataStream(i) = HammingSignal(ceil(i*baseband_dataRate/Fs));
        end
        Hamming_DataStream(Enc_signalLen) = Hamming_DataStream(Enc_signalLen -1);
        
        Hamming_BFSK_Signal_1 = Hamming_DataStream .* Enc_carrier_sig2;
        Hamming_BFSK_Signal_0 = (1 - Hamming_DataStream) .* Enc_carrier_sig;
        
        Hamming_BFSK_Signal = Hamming_BFSK_Signal_0 + Hamming_BFSK_Signal_1;

       
        % Generate noise BFSK
        BFSK_SignalPower = (norm(BFSK_Signal).^2)/signalLen;        
        BFSK_NoisePower_variance = BFSK_SignalPower ./ SNR;
        BFSK_Noise = sqrt(BFSK_NoisePower_variance/2) .*randn(1,signalLen);
        % Generate noise Cyclic Block Code BFSK
        CBC_BFSK_SignalPower = (norm(CBC_BFSK_Signal).^2)/Enc_signalLen;        
        CBC_BFSK_NoisePower_variance = CBC_BFSK_SignalPower ./ SNR;
        CBC_BFSK_Noise = sqrt(CBC_BFSK_NoisePower_variance/2) .*randn(1,Enc_signalLen);
        % Generate noise Linear Block Code BFSK
        LBC_BFSK_SignalPower = (norm(LBC_BFSK_Signal).^2)/Enc_signalLen;        
        LBC_BFSK_NoisePower_variance = LBC_BFSK_SignalPower ./ SNR;
        LBC_BFSK_Noise = sqrt(LBC_BFSK_NoisePower_variance/2) .*randn(1,Enc_signalLen);
        % Generate noise Hamming BFSK
        Hamming_BFSK_SignalPower = (norm(Hamming_BFSK_Signal).^2)/Enc_signalLen;        
        Hamming_BFSK_NoisePower_variance = Hamming_BFSK_SignalPower ./ SNR;
        Hamming_BFSK_Noise = sqrt(Hamming_BFSK_NoisePower_variance/2) .*randn(1,Enc_signalLen);

        % Transmit Signal BFSK
        BFSK_Signal_Received = BFSK_Signal + BFSK_Noise;     
        % Transmit Signal Cyclic Block Code: BFSK
        CBC_BFSK_Signal_Received = CBC_BFSK_Signal + CBC_BFSK_Noise;
        % Transmit Signal Linear Block Code: BFSK
        LBC_BFSK_Signal_Received = LBC_BFSK_Signal + LBC_BFSK_Noise;
        % Transmit Signal Hamming: BFSK
        Hamming_BFSK_Signal_Received = Hamming_BFSK_Signal + Hamming_BFSK_Noise;
        
        %{
        %Coherent detection BFSK
        BFSK_demod = BFSK_Signal_Received.*(2.*carrier_sig); %square law device (detection)
        %Coherent detection Cyclic Block Code: BFSK
        CBC_BFSK_demod = CBC_BFSK_Signal_Received.*(2.*Enc_carrier_sig);
        %Coherent detection Linear Block Code: BFSK
        LBC_BFSK_demod = LBC_BFSK_Signal_Received.*(2.*Enc_carrier_sig);
        %Coherent detection Hamming: BFSK
        Hamming_BFSK_demod = Hamming_BFSK_Signal_Received.*(2.*Enc_carrier_sig); 
        %}

        % Filtering of the demodulated signal BFSK
        BFSK_Filtered_0 = filtfilt(b_low, a_low, BFSK_Signal_Received);
        BFSK_Filtered_1 = filtfilt(b_high, a_high, BFSK_Signal_Received);
        % Filtering of the demodulated signal Cyclic Block Code: BFSK
        CBC_BFSK_Filtered_0 = filtfilt(b_low, a_low, CBC_BFSK_Signal_Received);
        CBC_BFSK_Filtered_1 = filtfilt(b_high, a_high, CBC_BFSK_Signal_Received);
        % Filtering of the demodulated signal Cyclic Block Code: BFSK
        LBC_BFSK_Filtered_0 = filtfilt(b_low, a_low, LBC_BFSK_Signal_Received);
        LBC_BFSK_Filtered_1 = filtfilt(b_high, a_high, LBC_BFSK_Signal_Received);
        % Filtering of the demodulated signal Cyclic Block Code: BFSK
        Hamming_BFSK_Filtered_0 = filtfilt(b_low, a_low, Hamming_BFSK_Signal_Received);
        Hamming_BFSK_Filtered_1 = filtfilt(b_high, a_high, Hamming_BFSK_Signal_Received);


        % Use the decision threshold logic for decoding of received signals
        % BFSK
        BFSK_Squared_0 = BFSK_Filtered_0.^2;
        BFSK_Squared_1 = BFSK_Filtered_1.^2;
        BFSK_Summation = BFSK_Squared_1 - BFSK_Squared_0;
        BFSK_Sampled = sample(BFSK_Summation,Ts,N_bits);      
        BFSK_Result = decision_logic(BFSK_Sampled,N_bits,0);

        % Use the decision threshold logic for decoding of received signals
        % Cyclic Block Code:
        CBC_BFSK_Squared_0 = CBC_BFSK_Filtered_0.^2;
        CBC_BFSK_Squared_1 = CBC_BFSK_Filtered_1.^2;
        CBC_BFSK_Summation = CBC_BFSK_Squared_1 - CBC_BFSK_Squared_0;
        CBC_BFSK_Sampled = sample(CBC_BFSK_Summation,Ts,Enc_N_bits);  
        CBC_BFSK_Result = decision_logic(CBC_BFSK_Sampled,Enc_N_bits,0);
        CBC_BFSK_DecodedResult = decode(CBC_BFSK_Result,n,k,'cyclic/binary',genpoly,trt);

        % Linear Block Code:
        LBC_BFSK_Squared_0 = LBC_BFSK_Filtered_0.^2;
        LBC_BFSK_Squared_1 = LBC_BFSK_Filtered_1.^2;
        LBC_BFSK_Summation = LBC_BFSK_Squared_1 - LBC_BFSK_Squared_0;
        LBC_BFSK_Sampled = sample(LBC_BFSK_Summation,Ts,Enc_N_bits);             
        LBC_BFSK_Result = decision_logic(LBC_BFSK_Sampled,Enc_N_bits,0);
        LBC_BFSK_DecodedResult = decode(LBC_BFSK_Result,n,k,'linear/binary',gen_mat,trt);

        % Hamming Block Code:
        Hamming_BFSK_Squared_0 = Hamming_BFSK_Filtered_0.^2;
        Hamming_BFSK_Squared_1 = Hamming_BFSK_Filtered_1.^2;
        Hamming_BFSK_Summation = Hamming_BFSK_Squared_1 - Hamming_BFSK_Squared_0;
        Hamming_BFSK_Sampled = sample(Hamming_BFSK_Summation,Ts,Enc_N_bits);    
        Hamming_BFSK_Result = decision_logic(Hamming_BFSK_Sampled,Enc_N_bits,0);
        Hamming_BFSK_DecodedResult = decode(Hamming_BFSK_Result,n,k,'hamming/binary');
        
        % Calculate the bit error rate performance BFSK
        BFSK_Error = 0;        
        for i = 1: N_bits - 1
            if(BFSK_Result(i) ~= Data(i))
                BFSK_Error = BFSK_Error + 1;
            end
        end     
        avg_BFSK_error = BFSK_Error/N_bits + avg_BFSK_error;

        % Calculate the bit error rate performance Cyclic Block Code BFSK
        CBC_BFSK_Error = 0;        
        for i = 1: N_bits - 1
            if(CBC_BFSK_DecodedResult(i) ~= Data(i))
                CBC_BFSK_Error = CBC_BFSK_Error + 1;
            end
        end     
        CBC_avg_BFSK_error = CBC_BFSK_Error/N_bits + CBC_avg_BFSK_error;

        % Calculate the bit error rate performance Linear Block Code BFSK
        LBC_BFSK_Error = 0;        
        for i = 1: N_bits - 1
            if(LBC_BFSK_DecodedResult(i) ~= Data(i))
                LBC_BFSK_Error = LBC_BFSK_Error + 1;
            end
        end     
        LBC_avg_BFSK_error = LBC_BFSK_Error/N_bits + LBC_avg_BFSK_error;

        % Calculate the bit error rate performance Hamming BFSK
        Hamming_BFSK_Error = 0;        
        for i = 1: N_bits - 1
            if(Hamming_BFSK_DecodedResult(i) ~= Data(i))
                Hamming_BFSK_Error = Hamming_BFSK_Error + 1;
            end
        end     
        Hamming_avg_BFSK_error = Hamming_BFSK_Error/N_bits + Hamming_avg_BFSK_error;
    
    end


    ER_BFSK(x) = (avg_BFSK_error / num_of_runs);
    CBC_ER_BFSK(x) = (CBC_avg_BFSK_error / num_of_runs);
    LBC_ER_BFSK(x) = (LBC_avg_BFSK_error / num_of_runs);
    Hamming_ER_BFSK(x) = (Hamming_avg_BFSK_error / num_of_runs);

end

% plot the result using  semilogy’ function
figure(1);
s1 = semilogy (SNR_db_Values_Array,ER_BFSK,'k-*');
hold on;
s2 = semilogy (SNR_db_Values_Array,CBC_ER_BFSK,'r-o');
hold on;
s3 = semilogy (SNR_db_Values_Array,LBC_ER_BFSK,'g-+');
hold on;
s4 = semilogy (SNR_db_Values_Array,Hamming_ER_BFSK,'b-x');
hold off
hold on

title('Error rate performance for BFSK, Cyclic Block Code,Linear Block Code, Hamming Code');
ylabel('Pe');
xlabel('Eb/No');
legend([s1(1),s2(1),s3(1),s4(1)],'BFSK','BFSK-CBC','BFSK-LBC','BFSK-Hamming','location','northeast');
hold off

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

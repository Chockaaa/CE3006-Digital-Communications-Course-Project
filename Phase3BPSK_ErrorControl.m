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
Enc_carrier_sig = A .* cos(2*pi*Fc*Enc_t);

% Gen LPF
% Assume a 6th order filter with cut-off frequency 0.2 in the function
[b_low, a_low] = butter(6, 0.2);
SNR_db_Values_Array = -50:5:50; %0:5:50;

signalLen = Fs* N_bits /baseband_dataRate + 1;
Enc_signalLen = Fs* Enc_N_bits /baseband_dataRate + 1;



ER_BPSK = zeros(1,length(SNR_db_Values_Array));
CBC_ER_BPSK = zeros(1,length(SNR_db_Values_Array));
LBC_ER_BPSK = zeros(1,length(SNR_db_Values_Array));
Hamming_ER_BPSK = zeros(1,length(SNR_db_Values_Array));

num_of_runs = 20;
for x = 1:length(SNR_db_Values_Array)
    SNR = (10.^(SNR_db_Values_Array(x)/10));   

    avg_BPSK_error = 0;
    CBC_avg_BPSK_error = 0;
    LBC_avg_BPSK_error = 0;
    Hamming_avg_BPSK_error = 0;

    % Generate data
    for j = 1 : num_of_runs     % Each SNR avg the error over 10 times
        Data = randi([0 1], 1 , N_bits);
        

        % BPSK: Data stream
        DataStream = zeros(1, signalLen);      
        
        for i = 1: signalLen - 1
            DataStream(i) = Data(ceil(i*baseband_dataRate/Fs));            
        end          
        DataStream(signalLen) = DataStream(signalLen - 1);
        BPSK_TransmitSignal = DataStream .* 2 - 1;
        BPSK_Signal = carrier_sig .* BPSK_TransmitSignal;

        % BPSK Cyclic Block Code: Data stream 
        CBCSignal = encode(Data,n,k,'cyclic/binary',genpoly);
        CBC_DataStream = zeros(1,Enc_signalLen);

        for i = 1:Enc_signalLen -1
            CBC_DataStream(i) = CBCSignal(ceil(i*baseband_dataRate/Fs));
        end
        CBC_DataStream(Enc_signalLen) = CBC_DataStream(Enc_signalLen -1);
        CBC_BPSK_TransmitSignal = CBC_DataStream .* 2 - 1;
        CBC_BPSK_Signal = Enc_carrier_sig .* CBC_BPSK_TransmitSignal;

        % BPSK Linear Block Code: Data stream 
        LBCSignal = encode(Data,n,k,'linear/binary',gen_mat);
        LBC_DataStream = zeros(1,Enc_signalLen);

        for i = 1:Enc_signalLen -1
            LBC_DataStream(i) = LBCSignal(ceil(i*baseband_dataRate/Fs));
        end
        LBC_DataStream(Enc_signalLen) = LBC_DataStream(Enc_signalLen -1);
        LBC_BPSK_TransmitSignal = LBC_DataStream .* 2 - 1;
        LBC_BPSK_Signal = Enc_carrier_sig .* LBC_BPSK_TransmitSignal;


        % BPSK Hamming: Data stream         
        HammingSignal = encode(Data,n,k,'hamming/binary');
        Hamming_DataStream = zeros(1,Enc_signalLen);
        for i = 1:Enc_signalLen -1
            Hamming_DataStream(i) = HammingSignal(ceil(i*baseband_dataRate/Fs));
        end
        Hamming_DataStream(Enc_signalLen) = Hamming_DataStream(Enc_signalLen -1);
        Hamming_BPSK_TransmitSignal = Hamming_DataStream .* 2 - 1;
        Hamming_BPSK_Signal = Enc_carrier_sig .* Hamming_BPSK_TransmitSignal;

       
        % Generate noise BPSK
        BPSK_SignalPower = (norm(BPSK_Signal).^2)/signalLen;        
        BPSK_NoisePower_variance = BPSK_SignalPower ./ SNR;
        BPSK_Noise = sqrt(BPSK_NoisePower_variance/2) .*randn(1,signalLen);
        % Generate noise Cyclic Block Code BPSK
        CBC_BPSK_SignalPower = (norm(CBC_BPSK_Signal).^2)/Enc_signalLen;        
        CBC_BPSK_NoisePower_variance = CBC_BPSK_SignalPower ./ SNR;
        CBC_BPSK_Noise = sqrt(CBC_BPSK_NoisePower_variance/2) .*randn(1,Enc_signalLen);
        % Generate noise Linear Block Code BPSK
        LBC_BPSK_SignalPower = (norm(LBC_BPSK_Signal).^2)/Enc_signalLen;        
        LBC_BPSK_NoisePower_variance = LBC_BPSK_SignalPower ./ SNR;
        LBC_BPSK_Noise = sqrt(LBC_BPSK_NoisePower_variance/2) .*randn(1,Enc_signalLen);
        % Generate noise Hamming BPSK
        Hamming_BPSK_SignalPower = (norm(Hamming_BPSK_Signal).^2)/Enc_signalLen;        
        Hamming_BPSK_NoisePower_variance = Hamming_BPSK_SignalPower ./ SNR;
        Hamming_BPSK_Noise = sqrt(Hamming_BPSK_NoisePower_variance/2) .*randn(1,Enc_signalLen);

        % Transmit Signal BPSK
        BPSK_Signal_Received = BPSK_Signal + BPSK_Noise;     
        % Transmit Signal Cyclic Block Code: BPSK
        CBC_BPSK_Signal_Received = CBC_BPSK_Signal + CBC_BPSK_Noise;
        % Transmit Signal Linear Block Code: BPSK
        LBC_BPSK_Signal_Received = LBC_BPSK_Signal + LBC_BPSK_Noise;
        % Transmit Signal Hamming: BPSK
        Hamming_BPSK_Signal_Received = Hamming_BPSK_Signal + Hamming_BPSK_Noise;
        
        
        %Coherent detection BPSK
        BPSK_demod = BPSK_Signal_Received.*(2.*carrier_sig); %square law device (detection)
        %Coherent detection Cyclic Block Code: BPSK
        CBC_BPSK_demod = CBC_BPSK_Signal_Received.*(2 .*Enc_carrier_sig);
        %Coherent detection Linear Block Code: BPSK
        LBC_BPSK_demod = LBC_BPSK_Signal_Received.*(2 .*Enc_carrier_sig);
        %Coherent detection Hamming: BPSK
        Hamming_BPSK_demod = Hamming_BPSK_Signal_Received.*(2 .*Enc_carrier_sig); 
        

        % Filtering of the demodulated signal BPSK
        BPSK_Filtered = filtfilt(b_low, a_low, BPSK_demod);
        % Filtering of the demodulated signal Cyclic Block Code: BPSK
        CBC_BPSK_Filtered = filtfilt(b_low, a_low, CBC_BPSK_demod);
        % Filtering of the demodulated signal Cyclic Block Code: BPSK
        LBC_BPSK_Filtered = filtfilt(b_low, a_low, LBC_BPSK_demod);
        % Filtering of the demodulated signal Cyclic Block Code: BPSK
        Hamming_BPSK_Filtered = filtfilt(b_low, a_low, Hamming_BPSK_demod);


        % Use the decision threshold logic for decoding of received signals
        % BPSK
        BPSK_Sampled = sample(BPSK_Filtered, Ts, N_bits);        
        BPSK_Result = decision_logic(BPSK_Sampled,N_bits,0);

        % Use the decision threshold logic for decoding of received signals
        % Cyclic Block Code:
        CBC_BPSK_Sampled = sample(CBC_BPSK_Filtered, Ts, Enc_N_bits);        
        CBC_BPSK_Result = decision_logic(CBC_BPSK_Sampled,Enc_N_bits,0);
        CBC_BPSK_DecodedResult = decode(CBC_BPSK_Result,n,k,'cyclic/binary',genpoly,trt);

        % Linear Block Code:
        LBC_BPSK_Sampled = sample(LBC_BPSK_Filtered, Ts, Enc_N_bits);        
        LBC_BPSK_Result = decision_logic(LBC_BPSK_Sampled,Enc_N_bits,0);
        LBC_BPSK_DecodedResult = decode(LBC_BPSK_Result,n,k,'linear/binary',gen_mat);

        % Hamming Block Code:
        Hamming_BPSK_Sampled = sample(Hamming_BPSK_Filtered, Ts, Enc_N_bits);        
        Hamming_BPSK_Result = decision_logic(Hamming_BPSK_Sampled,Enc_N_bits,0);
        Hamming_BPSK_DecodedResult = decode(Hamming_BPSK_Result,n,k,'hamming/binary');
        
        % Calculate the bit error rate performance BPSK
        BPSK_Error = 0;        
        for i = 1: N_bits - 1
            if(BPSK_Result(i) ~= Data(i))
                BPSK_Error = BPSK_Error + 1;
            end
        end     
        avg_BPSK_error = BPSK_Error/N_bits + avg_BPSK_error;

        % Calculate the bit error rate performance Cyclic Block Code BPSK
        CBC_BPSK_Error = 0;        
        for i = 1: N_bits - 1
            if(CBC_BPSK_DecodedResult(i) ~= Data(i))
                CBC_BPSK_Error = CBC_BPSK_Error + 1;
            end
        end     
        CBC_avg_BPSK_error = CBC_BPSK_Error/N_bits + CBC_avg_BPSK_error;

        % Calculate the bit error rate performance Linear Block Code BPSK
        LBC_BPSK_Error = 0;        
        for i = 1: N_bits - 1
            if(LBC_BPSK_DecodedResult(i) ~= Data(i))
                LBC_BPSK_Error = LBC_BPSK_Error + 1;
            end
        end     
        LBC_avg_BPSK_error = LBC_BPSK_Error/N_bits + LBC_avg_BPSK_error;

        % Calculate the bit error rate performance Hamming BPSK
        Hamming_BPSK_Error = 0;        
        for i = 1: N_bits - 1
            if(Hamming_BPSK_DecodedResult(i) ~= Data(i))
                Hamming_BPSK_Error = Hamming_BPSK_Error + 1;
            end
        end     
        Hamming_avg_BPSK_error = Hamming_BPSK_Error/N_bits + Hamming_avg_BPSK_error;
    
    end


    ER_BPSK(x) = (avg_BPSK_error / num_of_runs);
    CBC_ER_BPSK(x) = (CBC_avg_BPSK_error / num_of_runs);
    LBC_ER_BPSK(x) = (LBC_avg_BPSK_error / num_of_runs);
    Hamming_ER_BPSK(x) = (Hamming_avg_BPSK_error / num_of_runs);

end

% plot the result using  semilogy’ function
figure(1);
s1 = semilogy (SNR_db_Values_Array,ER_BPSK,'k-*');
hold on;
s2 = semilogy (SNR_db_Values_Array,CBC_ER_BPSK,'r-o');
hold on;
s3 = semilogy (SNR_db_Values_Array,LBC_ER_BPSK,'g-+');
hold on;
s4 = semilogy (SNR_db_Values_Array,Hamming_ER_BPSK,'b-x');
hold off
hold on

title('Error rate performance for BPSK, Cyclic Block Code,Linear Block Code, Hamming Code');
ylabel('Pe');
xlabel('Eb/No');
legend([s1(1),s2(1),s3(1),s4(1)],'BPSK','BPSK-CBC','BPSK-LBC','BPSK-Hamming','location','northeast');
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

SNR = 0:1:24;
N = 1e5; % # of symbols
Nt = 4; % # of Tx Antennas
Nr = 1; % # of Rx Antennas
M_EVGSM = 4; 
M_GSM = 8;
M_SM = 8;

BER_SM = SM_BERvsSNR_Simulation(SNR, N, Nt, Nr, M_SM);
BER_GSM = GSM_BERvsSNR_Simulation(SNR, N, Nt, Nr, M_GSM);
BER_EVGSM = EVGSM_BERvsSNR_Simulation(SNR, N, Nt, Nr, M_EVGSM);

semilogy(SNR, BER_SM, '-^', 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'Color', [0 0.45 0.74]);
hold on;
semilogy(SNR, BER_GSM, '-p', 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'Color', [0.85 0.33 0.10]); 
hold on;
semilogy(SNR, BER_EVGSM, '-h', 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'Color', [0.47 0.67 0.19]); hold on;
xlim([0 25]);
ylim([1e-2 1]);
set(gca, 'FontSize', 12);  
xlabel('SNR (dB)', 'FontSize', 14);
ylabel('Bit Error Rate (BER)', 'FontSize', 14);
title('BER Performance of SM, GSM, and EV-GSM vs SNR for 5 bits/s/Hz Transmission, N_r=1')
legend({'SM (8QAM, N_t = 4)', 'GSM (8QAM, N_t = 4)', 'EV-GSM (4QAM, N_t = 4)'}, 'Location','southwest', 'FontSize',12);
xticks(0:5:25);
yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1]);
grid on; grid minor;

function [BER_SM] = SM_BERvsSNR_Simulation(SNR, N, Nt, Nr, M_SM)
    BER_SM = zeros(1, length(SNR));
    antennaBits_BitNumber = log2(Nt);
    bits_Spatial = randi([0,1], N, antennaBits_BitNumber);
    bits_SpatialMask = bi2de(bits_Spatial, 'left-msb') + 1;
    APMBits_BitNumber = log2(M_SM);
    % 8QAM bits
    bits_Data = randi([0,1], N, APMBits_BitNumber);
    
    constellation = [
         1+1i;  1-1i; -1+1i; -1-1i;    % 4 outer QPSK-like points
         3+1i;  3-1i; -3+1i; -3-1i     % 4 extra amplitude-shifted points
    ];
    constellation = constellation / sqrt(mean(abs(constellation).^2));
    
    bits_8QAMSymbols = de2bi(0:7, APMBits_BitNumber, 'left-msb');
    
    bpcu = antennaBits_BitNumber + APMBits_BitNumber;
    symbolIndices = bi2de(bits_Data, 'left-msb') + 1; 
    x = constellation(symbolIndices);
    Es = mean(abs(x).^2);

    for i_SNR = 1:length(SNR)
        gamma = 10.^(SNR(i_SNR)/10);
        N0  = Es / gamma;
        
        % Channel Matrix
        H = sqrt(0.5) * (randn(Nr,Nt,N) + 1i*randn(Nr,Nt,N));
       
        n_awgn = sqrt(N0/2) * (randn(Nr,N) + 1i*randn(Nr,N));
        
        y = complex(zeros(Nr, N));

        for symbol = 1:N
            Index_Tx = bits_SpatialMask(symbol);
            h = H(:, Index_Tx, symbol); 

            xs = x(symbol);
            y(:, symbol) = h*xs + n_awgn(:,symbol);
        end
        
        decoded_SpatialBits = zeros(N, antennaBits_BitNumber);
        decoded_APMBits = zeros(N, APMBits_BitNumber);

        for symbol = 1:N
            minVal = inf; 
            for i_Tx = 1:Nt
                h = H(:, i_Tx, symbol);  

                y_Expected = h * constellation.';
                metric = sum(abs(y(:,symbol) - y_Expected).^2, 1);   
                [val, symbolIdx] = min(metric);

                if val < minVal
                    minVal = val;
                    minIndex_Tx = i_Tx;  
                    minIndex_Constellation = symbolIdx; 
                end
                
            end
            
            decoded_SpatialBits(symbol,:) = de2bi(minIndex_Tx-1, antennaBits_BitNumber, 'left-msb');
            decoded_APMBits(symbol,  :) = bits_8QAMSymbols(minIndex_Constellation, :); 
        end   

        Error_Count = sum(decoded_APMBits ~= bits_Data, 'all') + sum(decoded_SpatialBits ~= bits_Spatial,'all');
        BER_SM(i_SNR) = Error_Count / (bpcu * N);
        fprintf('SM -> SNR=%d dB | BER_SM=%d\n', SNR(i_SNR), BER_SM(i_SNR));
    end 
end

function [BER_GSM] = GSM_BERvsSNR_Simulation(SNR, N, Nt, Nr, M_GSM)
    BER_GSM = zeros(1, length(SNR));
    antennaPattern_BitNumber = 2;
    APMBits_BitNumber = log2(M_GSM);
    
    activeAntennaPatterns = [[1 2]; [1 3]; [1 4]; [2 3]];
    
    bits_Data = randi([0,1], N, APMBits_BitNumber);
    
    constellation = [
         1+1i;  1-1i; -1+1i; -1-1i;    % 4 outer QPSK-like points
         3+1i;  3-1i; -3+1i; -3-1i     % 4 extra amplitude-shifted points
    ];
    constellation = constellation / sqrt(mean(abs(constellation).^2));
    
    bits_8QAMSymbols = de2bi(0:7, APMBits_BitNumber, 'left-msb');
    
    bpcu = antennaPattern_BitNumber + APMBits_BitNumber;
    symbolIndices = bi2de(bits_Data, 'left-msb') + 1; 
    x = constellation(symbolIndices);
    Es = mean(abs(x).^2);
    
    bits_Spatial = randi([0, 1], N, antennaPattern_BitNumber);
    bits_SpatialMask = bi2de(bits_Spatial, 'left-msb')+1; 
    
    for i_SNR = 1:length(SNR)
        gamma = 10^(SNR(i_SNR)/10);
        N0 = Es / (gamma);
        
        % Channel Matrix 
        H = sqrt(0.5) * (randn(Nr,Nt,N) + 1i*randn(Nr,Nt,N));
    
        n_awgn = sqrt(N0/2) * (randn(Nr, N) + 1i*randn(Nr, N));
        
        y = complex(zeros(Nr, N));
        
        % Transmission
        for symbol = 1:N
            Index_Tx = bits_SpatialMask(symbol);
            antennas = activeAntennaPatterns(Index_Tx, :);
            activeAntennaNumber = numel(antennas);

            hActive = H(:, antennas, symbol);
            hNorm = sqrt(sum(sum(abs(hActive).^2)) / Nr);  
            y(:, symbol) = (hActive * ones(activeAntennaNumber,1)) * (x(symbol)/hNorm) + n_awgn(:, symbol);

        end
    
        decoded_SpatialBits = zeros(N, antennaPattern_BitNumber);
        decoded_APMBits = zeros(N, APMBits_BitNumber);
        % Decoding
        for symbol = 1:N
            min_val = inf; 
            for i_ActiveAntenna = 1:4
                antennas = activeAntennaPatterns(i_ActiveAntenna, :);
                activeAntennaNumber = numel(antennas);
                
                hActive = H(:, antennas, symbol);
                hNorm = sqrt(sum(sum(abs(hActive).^2)) / Nr);
        
                y_Expected = (hActive * ones(activeAntennaNumber,1)) * (constellation.' / hNorm);
        
                metric = sum(abs(y(:,symbol) - y_Expected).^2, 1);
                [val, symbol_idx] = min(metric);
    
                if val < min_val
                    min_val = val;
                    minIndex_Tx= i_ActiveAntenna;  
                    minIndex_Constellation = symbol_idx;  
                end            
            end
    
            % Decode the spatial bits and BPSK symbols
            decoded_SpatialBits(symbol, :) = de2bi(minIndex_Tx-1, antennaPattern_BitNumber, 'left-msb');
            decoded_APMBits(symbol, :) = bits_8QAMSymbols(minIndex_Constellation, :); 
        end
    
        % Calculate BER
        Error_Count_GSM = sum(decoded_SpatialBits ~= bits_Spatial, 'all') + sum(decoded_APMBits ~= bits_Data, 'all');
        BER_GSM(i_SNR) = Error_Count_GSM / (bpcu * N);
        fprintf('GSM -> SNR=%d dB| BER_GSM=%d\n', SNR(i_SNR), BER_GSM(i_SNR));
    end
end

function [BER_EVGSM] = EVGSM_BERvsSNR_Simulation(SNR, N, Nt, Nr, M_EVGSM)
    BER_EVGSM = zeros(1, length(SNR));
    activeAntenna_Selection = randi([1,2], N, 1);
    activeAntenna_BitNumber = 1;
    antennaPattern_Selection = randi([1,4], N, 1);
    antennaPattern_BitNumber = 2;
    APMBits_BitNumber = log2(M_EVGSM);
    bpcu = activeAntenna_BitNumber+antennaPattern_BitNumber+APMBits_BitNumber; 

    activeAntenna_Patterns = {
        {[1 2], [1 3], [1 4], [2 3], [2 4], [3 4]}, ...
        {[1 2 3], [1 2 4], [1 3 4], [2 3 4]}
    };
    
    bits_ActiveAntenna = de2bi(activeAntenna_Selection - 1, activeAntenna_BitNumber, 'left-msb');
    bits_AntennaPattern = de2bi(antennaPattern_Selection - 1, antennaPattern_BitNumber, 'left-msb');
    
    % 4QAM bits
    bits_Data = randi([0,1], N, APMBits_BitNumber);
    I = 1 - 2*bits_Data(:,2);   % bit2 -> I
    Q = 1 - 2*bits_Data(:,1);   % bit1 -> Q
    x = (I + 1i*Q)/sqrt(2);     % 00->++, 01->-+, 11->--, 10->+-  
    constellation = [1+1i; -1+1i; -1-1i; 1-1i]/sqrt(2);
    bits_4QAMSymbols = [0 0; 0 1; 1 1; 1 0];
    Es = mean(abs(x).^2);
    
    for i_SNR = 1:length(SNR) 
        gamma = 10^(SNR(i_SNR)/10);
        N0 = Es/(gamma);
    
        % Channel Matrix 
        H = sqrt(0.5) * (randn(Nr,Nt,N) + 1i*randn(Nr,Nt,N));
    
        n_awgn = sqrt(N0/2) * (randn(Nr, N) + 1i*randn(Nr, N));
        
        y = complex(zeros(Nr, N));
    
        % Transmission
        for symbol = 1:N
            antennas = activeAntenna_Patterns{activeAntenna_Selection(symbol)}{antennaPattern_Selection(symbol)};
            activeAntennaNumber = numel(antennas);

            hActive = H(:, antennas, symbol);
            
            y(:, symbol) = (hActive * ones(activeAntennaNumber,1)) * (x(symbol)) + n_awgn(:, symbol);
        end

        decodedBits_ActiveAntenna = zeros(N, activeAntenna_BitNumber);
        decodedBits_AntennaPattern = zeros(N, antennaPattern_BitNumber);
        decodedBits_APM = zeros(N, APMBits_BitNumber);

        % Decoding
        for symbol = 1:N
            minVal = inf;
            for antennaSelection = 1:2
                for antennaPattern = 1:4
                    antennas = activeAntenna_Patterns{antennaSelection}{antennaPattern};
                    activeAntennaNumber = numel(antennas);
                    
                    hActive = H(:, antennas, symbol);
        
                    y_Expected = (hActive * ones(activeAntennaNumber,1)) * (constellation.');   
                    metric = sum(abs(y(:,symbol) - y_Expected).^2, 1);       
        
                    [val, symbolIdx] = min(metric);
                    if val < minVal
                        minVal = val;
                        minIndex_AntennaSelection = antennaSelection;
                        minIndex_ActivePattern = antennaPattern;
                        minIndex_APM = symbolIdx;
                    end
                end
            end
     
            decodedBits_ActiveAntenna(symbol,:) = de2bi(minIndex_AntennaSelection-1, activeAntenna_BitNumber, 'left-msb');
            decodedBits_AntennaPattern(symbol,:) = de2bi(minIndex_ActivePattern-1, antennaPattern_BitNumber, 'left-msb');
            decodedBits_APM(symbol,:) = bits_4QAMSymbols(minIndex_APM,:);
        end
        errorCount = sum(decodedBits_ActiveAntenna ~= bits_ActiveAntenna, 'all')...
            + sum(decodedBits_AntennaPattern ~= bits_AntennaPattern, 'all')...
            + sum(decodedBits_APM ~= bits_Data, 'all');
        BER_EVGSM(i_SNR) = errorCount / (bpcu*N);
        fprintf('EV-GSM -> SNR=%2d dB | BER_EVGSM=%d\n', SNR(i_SNR), BER_EVGSM(i_SNR));
    end
end

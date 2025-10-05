SNR = 0:1:16;
N = 5e4; % # of symbols
Nt = 6; % # of Tx Antennas
Nr = 4; % # of Rx Antennas

BER_EVGSM7  = EVGSM_BERvsSNR_Simulation(SNR, N, Nt, Nr, 4, 7);
BER_EVGSM8  = EVGSM_BERvsSNR_Simulation(SNR, N, Nt, Nr, 8, 8);
BER_EVGSM9  = EVGSM_BERvsSNR_Simulation(SNR, N, Nt, Nr, 32, 9);
BER_EVGSM10  = EVGSM_BERvsSNR_Simulation(SNR, N, Nt, Nr, 64, 10);

semilogy(SNR, BER_EVGSM7,  '-h', 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'Color', [1.00 0.41 0.70]); hold on;  
semilogy(SNR, BER_EVGSM8,  '-h', 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'Color', [0.00 0.80 0.80]); hold on;  
semilogy(SNR, BER_EVGSM9,  '-h', 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'Color', [1.00 0.00 0.00]); hold on;  
semilogy(SNR, BER_EVGSM10, '-h', 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'Color', [0.00 0.70 0.00]); hold on;  
xlim([0 16]);
ylim([1e-5 1]);
set(gca, 'FontSize', 12);  
xlabel('SNR (dB)', 'FontSize', 14);
ylabel('Bit Error Rate (BER)', 'FontSize', 14);
title('BER Performance of EV-GSM vs SNR for Different Bit Rates')
legend({'K = 7 bits/sn/Hz', 'K = 8 bits/sn/Hz', 'K = 9 bits/sn/Hz', 'K = 10 bits/sn/Hz'}, 'Location','southwest', 'FontSize',12);
xticks(0:2:16);
yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1]);
grid on; grid minor;

function [BER_EVGSM] = EVGSM_BERvsSNR_Simulation(SNR, N, Nt, Nr, M_EVGSM, R) 
    BER_EVGSM = zeros(1, length(SNR));

    if R == 7 || R == 8
        selectionMax = 5;
        patternMax = 4;
        
        activeAntenna_BitNumber = 3;
        antennaPattern_BitNumber = 2;
    else
        selectionMax = 4;
        patternMax = 4;

        activeAntenna_BitNumber = 2;
        antennaPattern_BitNumber = 2;
    end
    activeAntenna_Selection = randi([1,selectionMax], N, 1);
    antennaPattern_Selection = randi([1,patternMax], N, 1);
    APMBits_BitNumber = log2(M_EVGSM);

    activeAntenna_Patterns = {
        {1, 2, 3, 4, 5, 6}, ...
        {[1 2], [1 3], [1 4], [2 3], [2 4], [2 5]}, ...
        {[1 2 3], [1 2 4], [1 3 4], [2 3 4], [3 4 5], [4 5 6]}, ...
        {[1 2 3 4], [1 3 4 5], [2 3 4 5], [3 4 5 6], [2 3 5 6], [1 2 4 5]},...
        {[1 2 3 4 5], [1 3 4 5 6], [2 3 4 5 6], [1 2 4 5 6], [1 2 3 4 6]}
    };
    
    bits_ActiveAntenna = de2bi(activeAntenna_Selection - 1, activeAntenna_BitNumber, 'left-msb');
    bits_AntennaPattern = de2bi(antennaPattern_Selection - 1, antennaPattern_BitNumber, 'left-msb');
    
    bits_Data = randi([0,1], N, APMBits_BitNumber);

    constellation = qammod(0:M_EVGSM-1, M_EVGSM, 'gray', 'UnitAveragePower', true);

    bits_16QAMSymbols = de2bi(0:M_EVGSM-1, APMBits_BitNumber, 'left-msb');
    
    bpcu = activeAntenna_BitNumber+antennaPattern_BitNumber+APMBits_BitNumber; 
    symbolIndices = bi2de(bits_Data, 'left-msb') + 1; 
    x = constellation(symbolIndices);
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
            for antennaSelection = 1:selectionMax
                for antennaPattern = 1:patternMax
                    antennas = activeAntenna_Patterns{antennaSelection}{antennaPattern};
                    activeAntennaNumber = numel(antennas);
                    
                    hActive = H(:, antennas, symbol);
        
                    y_Expected = (hActive * ones(activeAntennaNumber,1)) * (constellation);
   
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
            decodedBits_APM(symbol,:) = bits_16QAMSymbols(minIndex_APM,:);
        end
        errorCount = sum(decodedBits_ActiveAntenna ~= bits_ActiveAntenna, 'all')...
            + sum(decodedBits_AntennaPattern ~= bits_AntennaPattern, 'all')...
            + sum(decodedBits_APM ~= bits_Data, 'all');
        BER_EVGSM(i_SNR) = errorCount / (bpcu*N);
        fprintf('EV-GSM -> SNR=%2d dB | BER_EVGSM=%d\n', SNR(i_SNR), BER_EVGSM(i_SNR));
    end
end
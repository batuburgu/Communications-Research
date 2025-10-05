SNR = 0:1:10;
N = 1e4; % # of symbols
Nt = 16; % # of Tx Antennas
Nr = 4; % # of Rx Antennas
M = 2; % M-ary Modulation

[BER_GSMJ, BER_GSSK] = GSMGSSK_BERvsSNR_Simulation(SNR, Nt, Nr, M, N);
BER_JMVSM = JMVSM_BERvsSNR_Simulation(SNR, Nt, Nr, N);

semilogy(SNR, BER_GSMJ, '-o', 'Color', 'green', 'MarkerFaceColor', 'white', 'MarkerSize', 6);
hold on;
semilogy(SNR, BER_GSSK, '-square', 'Color', 'magenta', 'MarkerFaceColor', 'white', 'MarkerSize', 6);
hold on;
semilogy(SNR, BER_JMVSM, '-v', 'Color', 'red', 'MarkerFaceColor', 'white', 'MarkerSize', 6);
xlim([0 10]);
ylim([1e-5 1]);
xlabel('E_s/N_0 [dB]', 'FontSize', 14);
ylabel('Bit Error Rate', 'FontSize', 14);
title('GSM, GSSK and JM-VSM BER Comparison: N_t = 16, N_r = 4', 'FontSize', 16);
legend('GSM (Joint)', 'GSSK', 'JM-VSM', 'Location', 'SouthWest', 'FontSize', 12);
grid on;
grid minor;  
set(gca, 'FontSize', 12);  


function [BER_JMVSM] = JMVSM_BERvsSNR_Simulation(SNR, Nt, Nr, N)
    BER_JMVSM = zeros(1, length(SNR));
    antennaSelections = [
        1   9;
        4  12;
        7  15;
        2  16
    ];                         
    
    Constellation = {
        [ 1  1;  1 -1;  -1  1;  -1 -1];
        [ 1  1;  1 -1;  -1  1;  -1 -1];
        [ 1  1;  1 -1;  -1  1;  -1 -1];
        [ 1  1;  1 -1;  -1  1;  -1 -1]
    };     
    Es = 1;
    
    bpcu = 4;
    antennaBits_Number = 2;
    bits_Spatial = randi([0, 1], N, antennaBits_Number);
    tupel_Selection = randi([0, 1], N, 2);
    antennaPatternIndex = bi2de(bits_Spatial, 'left-msb')+1; 
    tupelSelectionIndex = bi2de(tupel_Selection, 'left-msb')+1; 
    
    sequence_Idx = (antennaPatternIndex-1)*4 + (tupelSelectionIndex-1);   
    sequence_Transmitted = de2bi(sequence_Idx, bpcu, 'left-msb');
    
    for i_SNR = 1:length(SNR)
        gamma = 10^(SNR(i_SNR)/10);
        N0 = Es / gamma;
        
        % Channel Matrix 
        H = sqrt(0.5) * (randn(Nr,Nt,N) + 1i*randn(Nr,Nt,N));
    
        n_awgn = sqrt(N0/2) * (randn(Nr, N) + 1i*randn(Nr, N));
        
        y = complex(zeros(Nr, N));
        
        % Transmission
        for symbol = 1:N
            patternIdx = antennaPatternIndex(symbol);        
            tupelIdx = tupelSelectionIndex(symbol);
            
            activeAntennas = antennaSelections(patternIdx,:); 
            activeAntennas = activeAntennas(activeAntennas~=0); 
            Na = numel(activeAntennas);     
        
            HActive = H(:, activeAntennas, symbol); 
            s = Constellation{patternIdx}(tupelIdx,:).';
    
            y(:, symbol) = sqrt(1/Na) .* (HActive * s) + n_awgn(:, symbol);
        end
    
        sequence_Decoded = zeros(N, 4);
        % Decoding
        for symbol = 1:N
            min_val = inf;
    
            for i_Antennas = 1:size(antennaSelections, 1)
                activeAntennas = antennaSelections(i_Antennas, :); 
                activeAntennas = activeAntennas(activeAntennas~=0); 
                Na  = numel(activeAntennas);
                
                HActive = H(:, activeAntennas, symbol);
                tuples = Constellation{i_Antennas}.';
    
                y_Expected = sqrt(1/Na) * HActive * tuples;
    
                metric = sum(abs(y(:, symbol) - y_Expected).^2, 1);
                [val, symbol_idx] = min(metric(:));
    
                if val < min_val
                    min_val = val;
                    MinIndex_Tx = i_Antennas;  
                    MinIndex_Constellation = symbol_idx;  
                end            
            end    
            selectedIndex = (MinIndex_Tx - 1)*4 + (MinIndex_Constellation-1);     
            sequence_Decoded(symbol,:) = de2bi(selectedIndex, bpcu, 'left-msb');
        end
        Error_Count_JMVSM = sum(sequence_Decoded ~= sequence_Transmitted, 'all');
        BER_JMVSM(i_SNR) = Error_Count_JMVSM / (bpcu * N);
        fprintf('SNR=%2d dB | BER_JMVSM =%.4g\n', SNR(i_SNR), BER_JMVSM(i_SNR));
    end
end

function [BER_GSMJ, BER_GSSK] = GSMGSSK_BERvsSNR_Simulation(SNR, Nt, Nr, M, N)
    BER_GSMJ = zeros(1, length(SNR));
    BER_GSSK = zeros(1, length(SNR));
    antennaBits_Number = 3;
    bits_Spatial = randi([0, 1], N, antennaBits_Number);
    bits_SpatialMask = bi2de(bits_Spatial, 'left-msb')+1; 
    antennaTriplets = [ 1 2 3;
                        9 10 11;
                        8 12 13;
                        5 12 13;
                        13 14 15;
                        5 9 12;
                        3 7 11;
                        6 8 10];
    bpcu = 4;
    % BPSK bits
    APMBits_Number = log2(M);
    bits_Data = randi([0, 1], N, APMBits_Number);
    x = (1 - 2 * bits_Data);   % 0->+, 1->-
    constellation = [1, -1];   % BPSK Constellation
    Es = mean(abs(x).^2);
    bits_BPSKSymbols = [0; 1];

    sequence_idx = (bits_SpatialMask - 1) * 2 + bits_Data;   % 0..15
    sequence_Transmitted = de2bi(sequence_idx, bpcu, 'left-msb');
    
    for i_SNR = 1:length(SNR)
        gamma = 10^(SNR(i_SNR)/10);
        N0 = Es / gamma;
        
        % Channel Matrix 
        H = sqrt(0.5) * (randn(Nr,Nt,N) + 1i*randn(Nr,Nt,N));
    
        n_awgn = sqrt(N0/2) * (randn(Nr, N) + 1i*randn(Nr, N));
        
        y = complex(zeros(Nr, N));
        
        % Transmission
        for symbol = 1:N
            Index_Tx = bits_SpatialMask(symbol);
            antennaTriplet_Active = antennaTriplets(Index_Tx, :);
            
            y(:, symbol) = sqrt(1/3) * sum(H(:, antennaTriplet_Active, symbol), 2) * x(symbol) + n_awgn(:, symbol);
        end
    
        SpatialBits_Decoded = zeros(N, antennaBits_Number);
        APMBits_Decoded = zeros(N, APMBits_Number);
        sequence_Decoded = zeros(N, bpcu);
        % Decoding
        for symbol = 1:N
            min_val = inf;
            
            for i_Triplet = 1:size(antennaTriplets, 1)
                triplet = antennaTriplets(i_Triplet, :); 
    
                y_Expected = sqrt(1/3) * sum(H(:, triplet, symbol), 2) * constellation;
    
                metric = sum(abs(y(:, symbol) - y_Expected).^2);
                [val, symbol_idx] = min(metric(:));
    
                if val < min_val
                    min_val = val;
                    MinIndex_Tx= i_Triplet;  
                    MinIndex_Constellation = symbol_idx;  
                end            
            end
    
            % Decode the spatial bits and BPSK symbols
            SpatialBits_Decoded(symbol, :) = de2bi(MinIndex_Tx-1, antennaBits_Number, 'left-msb');
            APMBits_Decoded(symbol, :) = bits_BPSKSymbols(MinIndex_Constellation, :); 

            data_bit     = MinIndex_Constellation - 1;           % 0 for +1, 1 for -1
            selectedIndex = (MinIndex_Tx - 1)*2 + data_bit;     
            sequence_Decoded(symbol,:) = de2bi(selectedIndex, bpcu, 'left-msb');
        end
    
        % Calculate BER
        Error_Count_GSSK = sum(SpatialBits_Decoded ~= bits_Spatial, 'all');
        Error_Count_GSMJ = sum(sequence_Decoded ~= sequence_Transmitted, 'all');
        BER_GSSK(i_SNR) = Error_Count_GSSK / (antennaBits_Number * N);
        BER_GSMJ(i_SNR) = Error_Count_GSMJ / (bpcu * N);
        fprintf('-> SNR=%d | BER_GSMJ=%d -> BER_GSSK=%d\n', SNR(i_SNR), BER_GSMJ(i_SNR), BER_GSSK(i_SNR));
    end
end



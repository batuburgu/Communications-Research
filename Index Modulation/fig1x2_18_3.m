SNR_Fig1a = 12; % Signal to Noise Ratio (SNR)
SNR_Fig12 = 0:1:40;
N = 1e6; % # of symbols
Nt_Fig1 = 2; % # of Tx Antennas
Nr_Fig1a = 1;
Nr_Fig1b = [1, 2]; % # of Rx Antennas
M = 4; % M-ary Modulation

theta_deg_fig1a = 0:1:90;
theta_deg_fig1b = 0:1:22;
theta_fig1a = deg2rad(theta_deg_fig1a);
theta_fig1b = deg2rad(theta_deg_fig1b);

%BER_Fig1a = DESM_BERvsPhaseRotation_Simulation(SNR_Fig1a, Nt_Fig1, Nr_Fig1a, M, N, theta_fig1a);

%[Optimum_Theta, Optimum_BER] = DESM_OptThetaOptBERvsSNR_Simulation(SNR_Fig12, Nt_Fig1, Nr_Fig1b, M, N, theta_fig1b, theta_deg_fig1b);
SM_BER = SM_BERvsSNR_Simulation(SNR_Fig12, Nt_Fig1, Nr_Fig1b, M, N);

figure('Name', 'Fig1.a');
semilogy(theta_deg_fig1a, BER_Fig1a, '-v', ...
    'LineWidth', 1.6, 'MarkerSize', 5, ...
    'Color', 'b', ...             
    'MarkerEdgeColor', 'g');      
grid on; grid minor; box on;
xlim([0 90]);                                   
xticks(0:10:90);
ylim([5e-2 1e-1]);                               
xlabel('Phase rotation \theta (degrees)');
ylabel('Average BER');
title('Fig.1: The ABER Versus the Phase Rotation for DE-SM when SNR = 12 dB');
legend({'DE-SM, N_t = 2, N_r = 1, 4-QAM'}, 'Location','best');

%SM Sim 1e6, DE-SM Sim 1e5 Symbols
figure('Name', "Fig.2")
semilogy(SNR_Fig12, Optimum_BER(1,:), '-o', ...
    'Color','g', 'MarkerFaceColor','r', 'MarkerEdgeColor','g', ...
    'LineWidth',1.6, 'MarkerSize',5);
hold on;
semilogy(SNR_Fig12, SM_BER(1,:), '--o', ...
    'Color','g', 'MarkerFaceColor','r', 'MarkerEdgeColor','g', ...
    'LineWidth',1.6, 'MarkerSize',5);
hold on;
semilogy(SNR_Fig12, Optimum_BER(2,:), '-s', ...   
    'Color','c', 'MarkerFaceColor','r', 'MarkerEdgeColor','c', ...
    'LineWidth',1.6, 'MarkerSize',5);
hold on;
semilogy(SNR_Fig12, SM_BER(2,:), '--s', ...
    'Color','c', 'MarkerFaceColor','r', 'MarkerEdgeColor','c', ...
    'LineWidth',1.6, 'MarkerSize',5);
grid on; grid minor; box on;
xlim([0 40]); xticks(0:5:40);
ylim([1e-4 1]);

xlabel('SNR [dB]');
ylabel('Average BER');
title('DE-SM & SM: ABER vs SNR (\theta = \theta_{opt} for DE-SM)');

legend({'DE-SM, N_t = 2, N_r = 1, 4-QAM', ...
        'SM, N_t = 2, N_r = 1, 4-QAM',...
        'DE-SM, N_t = 2, N_r = 2, 4-QAM', ...
        'SM, N_t = 2, N_r = 2, 4-QAM'},...
        'Location','best');

function [BER] = SM_BERvsSNR_Simulation(SNR, Nt, Nr, M, N)

    BER = zeros(length(Nr), length(SNR));
    AntennaBits_Number = log2(Nt);
    Bits_Spatial = randi([0,1], N, AntennaBits_Number);
    Bits_Spatial_Mask = bi2de(Bits_Spatial, 'left-msb') + 1;
    
    % 4QAM bits
    APMBits_Number = log2(M);
    Bits_Data = randi([0,1], N, APMBits_Number);
    I = 1 - 2*Bits_Data(:,2);   % bit2 -> I
    Q = 1 - 2*Bits_Data(:,1);   % bit1 -> Q
    x = (I + 1i*Q) / Nt;     % 00->++, 01->-+, 11->--, 10->+-  
    Constellation = [1+1i, -1+1i, -1-1i, 1-1i] / Nt;
    Es = mean(abs(x).^2);
    Bits_4QAMSymbols = [0 0; 0 1; 1 1; 1 0];

    for i_Nr = 1:length(Nr)
        Current_Nr = Nr(i_Nr);
        for i_SNR = 1:length(SNR)
            % Channel Matrix
            H = complex(zeros(Current_Nr, Nt, N));
            for i_tx = 1:Nt
                H(:, i_tx, :) = sqrt(1/2) .* (randn(Current_Nr,N) + 1i.*randn(Current_Nr,N));    
            end
            gamma = 10.^(SNR(i_SNR)/10);
            N0  = Es / gamma;
            n_awgn = sqrt(N0/2) * (randn(Current_Nr,N) + 1i*randn(Current_Nr,N));
            
            y = complex(zeros(Current_Nr, N));

            for symbol = 1:N
                Index_Tx = Bits_Spatial_Mask(symbol);
                h_tx = H(:, Index_Tx, symbol); 

                xs = x(symbol);
                y(:, symbol) = h_tx*xs + n_awgn(:,symbol);
            end
            
            SpatialBits_Decoded = zeros(N, AntennaBits_Number);
            APMBits_Decoded = zeros(N, APMBits_Number);

            for symbol = 1:N
                min_val = inf; 
                for i_tx = 1:Nt
                    h = H(:, i_tx, symbol);  
                    
                    for i_Const = 1:length(Constellation)
                        metric = sum(abs(y(:,symbol) - h*Constellation(i_Const)).^2);
            
                        if metric < min_val
                            min_val = metric;
                            MinIndex_Tx = i_tx;  
                            MinIndex_Const = i_Const; 
                        end
                    end
                end
                
                SpatialBits_Decoded(symbol,:) = de2bi(MinIndex_Tx-1, AntennaBits_Number, 'left-msb');
                APMBits_Decoded(symbol,  :) = Bits_4QAMSymbols(MinIndex_Const, :); 
            end   
            
            Error_Count = sum(APMBits_Decoded ~= Bits_Data, 'all') + sum(SpatialBits_Decoded ~= Bits_Spatial,'all');
            BER(i_Nr,i_SNR) = Error_Count / ((APMBits_Number + AntennaBits_Number) * N);
            fprintf('Antenna Index=%d -> SNR=%d dB -> BER=%.4e\n', i_Nr, SNR(i_SNR), BER(i_Nr, i_SNR));
        end
    end
end

function [Optimum_Theta, Optimum_BER] = DESM_OptThetaOptBERvsSNR_Simulation(SNR, Nt, Nr, M, N, theta, theta_deg)
    
    Optimum_Theta = zeros(length(Nr), length(SNR));
    Optimum_BER = zeros(length(Nr), length(SNR));
    AntennaBits_Number = log2(Nt);
    Bits_Spatial = randi([0,1], N, AntennaBits_Number);
    Bits_Spatial_Mask = bi2de(Bits_Spatial, 'left-msb') + 1;
    
    % 4QAM bits
    APMBits_Number = log2(M);
    Bits_Data = randi([0,1], N, APMBits_Number);
    I = 1 - 2*Bits_Data(:,2);   % bit2 -> I
    Q = 1 - 2*Bits_Data(:,1);   % bit1 -> Q
    x = (I + 1i*Q) / Nt;     % 00->++, 01->-+, 11->--, 10->+-  
    Constellation = [1+1i, -1+1i, -1-1i, 1-1i] / Nt;
    Es = mean(abs(x));
    Bits_4QAMSymbols = [0 0; 0 1; 1 1; 1 0];
    
    for i_Nr = 1:length(Nr)
        Current_Nr = Nr(i_Nr);

        for i_SNR = 1:length(SNR)
            BER = zeros(1,length(theta));
            % Channel Matrix
            H = complex(zeros(Current_Nr, Nt, N));
            for i_tx = 1:Nt
                H(:, i_tx, :) = sqrt(1/2) .* (randn(Current_Nr,N) + 1i.*randn(Current_Nr,N));    
            end
            gamma = 10.^(SNR(i_SNR)/10);
            N0  = Es / gamma;
            n_awgn = sqrt(N0/2) * (randn(Current_Nr,N) + 1i*randn(Current_Nr,N));
            
            for i = 1:length(theta)
                y = complex(zeros(Current_Nr, N));
                Constellation_Shifted = Constellation.*exp(1i*theta(i));
                
                % Possible Constellations at the receiver
                [Constellation_T1, Constellation_T2] = ndgrid(Constellation_Shifted,Constellation_Shifted);
                Constellation_PossibilitiesT1 = real(Constellation_T1) + 1i*imag(Constellation_T2);
                Constellation_PossibilitiesT2 = real(Constellation_T2) + 1i*imag(Constellation_T1);
            
                % Transmission
                for symbol = 1:2:N
                    Index_Tx = Bits_Spatial_Mask(symbol: symbol+1);
            
                    s1 = x(symbol).*exp(1i*theta(i));
                    s2 = x(symbol+1).*exp(1i*theta(i));
            
                    y(:, symbol) = H(:,Index_Tx(1),symbol).*(real(s1) + 1i*imag(s2)) + n_awgn(:,symbol); % First time period
                    y(:, symbol+1) = H(:,Index_Tx(2),symbol+1).*(real(s2) + 1i*imag(s1)) + n_awgn(:,symbol); % Second time period 
                end
            
                SpatialBits_Decoded = zeros(N, AntennaBits_Number);
                APMBits_Decoded = zeros(N, APMBits_Number);
                MinIndex_Tx = zeros(Nt, N);
            
                % Receiving and Decoding
                for symbol = 1:2:N 
                    min_val = inf;
            
                    if Nr(i_Nr) < 2
                        for i_txT1 = 1:Nt
                            metric_T1 = abs(y(:, symbol) - (H(:,i_txT1, symbol).*Constellation_PossibilitiesT1)).^2;
                
                            for i_txT2 = 1:Nt
                                metric_T2 = abs(y(:, symbol+1) - (H(:,i_txT2, symbol+1).*Constellation_PossibilitiesT2)).^2;
                                metric = metric_T1 + metric_T2;
                                
                                [val, idx] = min(metric(:));   
                                if val < min_val
                                    min_val = val;
                                    [x_T1, x_T2] = ind2sub(size(metric),idx);
                                    MinIndex_Tx(:,symbol) = [i_txT1, i_txT2];
                                end
                            end
                        end  
    
                    else
                        Constellation_PossibilitiesT1 = reshape(Constellation_PossibilitiesT1,1,M,M);
                        Constellation_PossibilitiesT2 = reshape(Constellation_PossibilitiesT2,1,M,M);
                        for i_txT1 = 1:Nt
                            metric_T1 = sum(abs(y(:, symbol) - (H(:,i_txT1, symbol).*Constellation_PossibilitiesT1)).^2);
                            metric_T1 = squeeze(metric_T1);
                
                            for i_txT2 = 1:Nt
                                metric_T2 = sum(abs(y(:, symbol+1) - (H(:,i_txT2, symbol+1).*Constellation_PossibilitiesT2)).^2);
                                metric_T2 = squeeze(metric_T2);
                                metric = metric_T1 + metric_T2;
                                
                                [val, idx] = min(metric(:));   
                                if val < min_val
                                    min_val = val;
                                    [x_T1, x_T2] = ind2sub(size(metric),idx);
                                    MinIndex_Tx(:,symbol) = [i_txT1, i_txT2];
                                end
                            end
                        end 
                    end
            
                    APMBits_Decoded(symbol,  :) = Bits_4QAMSymbols(x_T1, :); 
                    APMBits_Decoded(symbol+1,:) = Bits_4QAMSymbols(x_T2, :); 
                
                    SpatialBits_Decoded(symbol,:) = de2bi(MinIndex_Tx(1,symbol)-1, AntennaBits_Number, 'left-msb');
                    SpatialBits_Decoded(symbol+1,:) = de2bi(MinIndex_Tx(2,symbol)-1, AntennaBits_Number, 'left-msb');
                end
                
                Error_Count = sum(APMBits_Decoded ~= Bits_Data, 'all') + sum(SpatialBits_Decoded ~= Bits_Spatial,'all');
                BER(i) = Error_Count / ((APMBits_Number + AntennaBits_Number) * N);
                fprintf('SNR=%d dB -> ThetaIndex=%d -> BER=%.4e\n',SNR(i_SNR), i, BER(i));
            end
            [~, index] = min(BER);
            Optimum_BER(i_Nr, i_SNR) = BER(index);
            Optimum_Theta(i_Nr, i_SNR) = theta_deg(index);
            fprintf('ThetaOptDegree=%d -> BEROpt=%d -> SNR=%.4e\n',Optimum_Theta(i_Nr, i_SNR), Optimum_BER(i_Nr, i_SNR), SNR(i_SNR));
        end
    end
end

function [BER] = DESM_BERvsPhaseRotation_Simulation(SNR, Nt, Nr, M, N, theta)
    BER = zeros(1,length(theta));

    AntennaBits_Number = log2(Nt);
    Bits_Spatial = randi([0,1], N, AntennaBits_Number);
    Bits_Spatial_Mask = bi2de(Bits_Spatial, 'left-msb') + 1;
    
    % 4QAM bits
    APMBits_Number = log2(M);
    Bits_Data = randi([0,1], N, APMBits_Number);
    I = 1 - 2*Bits_Data(:,2);   % bit2 -> I
    Q = 1 - 2*Bits_Data(:,1);   % bit1 -> Q
    x = (I + 1i*Q) / Nt;     % 00->++, 01->-+, 11->--, 10->+-  
    Constellation = [1+1i, -1+1i, -1-1i, 1-1i] / Nt;
    Es = mean(abs(x).^2);
    Bits_4QAMSymbols = [0 0; 0 1; 1 1; 1 0];
    
    % Channel Matrix
    H = complex(zeros(Nr, Nt, N));
    
    for i = 1:length(theta)
        gamma = 10.^(SNR/10);
        N0  = Es / gamma;
        y = complex(zeros(Nr, N));
        Constellation_Shifted = Constellation.*exp(1i*theta(i));
        
        % Possible Constellations at the receiver
        [Constellation_T1, Constellation_T2] = ndgrid(Constellation_Shifted,Constellation_Shifted);
        Constellation_PossibilitiesT1 = real(Constellation_T1) + 1i*imag(Constellation_T2);
        Constellation_PossibilitiesT2 = real(Constellation_T2) + 1i*imag(Constellation_T1);
    
        % Channel Matrix
        for i_tx = 1:Nt
            H(:, i_tx, :) = sqrt(1/2) .* (randn(Nr,N) + 1i.*randn(Nr,N));    
        end
        
        % Transmission
        for symbol = 1:2:N
            Index_Tx = Bits_Spatial_Mask(symbol: symbol+1);
    
            s1 = x(symbol).*exp(1i*theta(i));
            s2 = x(symbol+1).*exp(1i*theta(i));
    
            y(:, symbol) = H(:,Index_Tx(1),symbol).*(real(s1) + 1i*imag(s2)) + sqrt(N0/2) * (randn(Nr,1) + 1i*randn(Nr,1)); % First time period
            y(:, symbol+1) = H(:,Index_Tx(2),symbol+1).*(real(s2) + 1i*imag(s1)) + sqrt(N0/2) * (randn(Nr,1) + 1i*randn(Nr,1)); % Second time period 
        end
    
        SpatialBits_Decoded = zeros(N, AntennaBits_Number);
        APMBits_Decoded = zeros(N, APMBits_Number);
        MinIndex_Tx = zeros(Nt, N);
    
        % Receiving and Decoding
        for symbol = 1:2:N 
            min_val = inf;
    
            for i_txT1 = 1:Nt
                metric_T1 = abs(y(:, symbol) - (H(:,i_txT1, symbol).*Constellation_PossibilitiesT1)).^2;
    
                for i_txT2 = 1:Nt
                    metric_T2 = abs(y(:, symbol+1) - (H(:,i_txT2, symbol+1).*Constellation_PossibilitiesT2)).^2;
                    metric = metric_T1 + metric_T2;
                    
                    [val, idx] = min(metric(:));   
                    if val < min_val
                        min_val = val;
                        [x_T1, x_T2] = ind2sub(size(metric),idx);
                        MinIndex_Tx(:,symbol) = [i_txT1, i_txT2];
                    end
                end
            end       
    
            APMBits_Decoded(symbol,  :) = Bits_4QAMSymbols(x_T1, :); 
            APMBits_Decoded(symbol+1,:) = Bits_4QAMSymbols(x_T2, :); 
        
            SpatialBits_Decoded(symbol,:) = de2bi(MinIndex_Tx(1,symbol)-1, AntennaBits_Number, 'left-msb');
            SpatialBits_Decoded(symbol+1,:) = de2bi(MinIndex_Tx(2,symbol)-1, AntennaBits_Number, 'left-msb');
        end
        
        Error_Count = sum(APMBits_Decoded ~= Bits_Data, 'all') + sum(SpatialBits_Decoded ~= Bits_Spatial,'all');
        BER(i) = Error_Count / ((APMBits_Number + AntennaBits_Number) * N);
        fprintf('ThetaIndex=%d -> BER=%.4e\n', i, BER(i));
    end
end

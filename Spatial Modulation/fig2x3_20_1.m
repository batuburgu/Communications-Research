SNR_NMSE = [15,20];
SNR = 10:0.5:18; 
N = 1e6; % # of symbols
Nt = 4; % # of Tx Antennas
Nr = 4; % # of Rx Antennas
M = 4; % M-ary Modulation

% Rayleigh Fading with alpha-mu distribution
alpha = 2;
mu = 1;

omega = 1; % Expected Value of Fading Channel

% k^2 Values
k_NMSE = linspace(-20, -10, 21); 
k_sq = 10^(-16/10);

% Path Loss 
h_l = 1; 

% Antenna Misalignment Fading
h_mis = 1; 


[NMSE_LS, NMSE_MLS] = Simulation_NMSE(SNR_NMSE, k_NMSE, h_mis, h_l, Nt, Nr, N, alpha, mu, omega);
[BER_PerfectCSI] = Simulation_BER(SNR, k_sq, h_mis, h_l, Nt, Nr, N, M, alpha, mu, omega, 1);
[BER_ActualCSI] = Simulation_BER(SNR, k_sq, h_mis, h_l, Nt, Nr, N, M, alpha, mu, omega, 2);

% Channel Estimation NMSE Plots
figure;
plot(k_NMSE, NMSE_MLS(1,:), '-', 'Color', 'b', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); 
hold on;
plot(k_NMSE, NMSE_LS(1,:), '-s', 'Color', 'm', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm');
plot(k_NMSE, NMSE_MLS(2,:), '--', 'Color', 'g', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
plot(k_NMSE, NMSE_LS(2,:), '-o', 'Color', 'r', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
ylim([5*10^-3, 0.14]);
xlim([-20, -10]);
title('Figure 2: Channel Estimation Accuracy Comparison Between Proposed MLS and LS Estimators', ...
    'FontSize', 14, 'FontWeight', 'bold');
xlabel('k^2 [dB]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Normalized MSE', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Proposed MLS Estimator SNR = 15dB',...
        'Classical LS Estimator SNR = 15dB', ...
        'Proposed MLS Estimator SNR = 20dB', ...
        'Classical LS Estimator SNR = 20dB'}, ...
        'Location', 'northwest', 'FontSize', 12, 'FontWeight', 'bold', 'Box', 'off');
grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.5);
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

% BER Plots
figure;
semilogy(SNR, BER_PerfectCSI, 'b-diamond', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
hold on;
semilogy(SNR, BER_ActualCSI, 'r-o', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
ylim([10^-5, 10^-1]);
xlim([10, 17]);
title('Figure 3: Performance Comparison between THz-SM systems with proposed and conventional designs', ...
    'FontSize', 16, 'FontWeight', 'bold');
xlabel('SNR_r [dB]', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('BER', 'FontSize', 14, 'FontWeight', 'bold');
legend({'Proposed Detector, Perfect CSI, k_t=k_r=-16 dB', ...
        'Proposed Detector, Actual CSI, k_t=k_r=-16 dB'}, ...
        'Location', 'northeast', 'FontSize', 12, 'FontWeight', 'bold', 'Box', 'off');
grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.6);
set(gca, 'FontSize', 12, 'FontWeight', 'bold');


function [NMSE_LS, NMSE_MLS] = Simulation_NMSE(SNR_NMSE, k_dB, h_mis, h_l, Nt, Nr, N, alpha, mu, omega)
    linear_SNR = 10.^(SNR_NMSE/10);
    NMSE_LS  = zeros(length(SNR_NMSE), length(k_dB));
    NMSE_MLS = zeros(length(SNR_NMSE), length(k_dB));

    % Pilot Symbol
    sp = 1;  
    P = abs(sp)^2;  

    antenna_bits = log2(Nt);
    spatial_bits = randi([0,1], N, antenna_bits);
    mask = bi2de(spatial_bits, 'left-msb') + 1;  

    h = complex(zeros(Nr, Nt, N));

    for SNR_i = 1:length(SNR_NMSE)
        gamma = linear_SNR(SNR_i);
        N0  = P / gamma;         
        for k_i = 1:length(k_dB)
            k_sq = 10^(k_dB(k_i)/10);

            for i_tx = 1:Nt
                h_f = AlphaMuGenerator(alpha, mu, omega, Nr, N);  
                h(:, i_tx, :) = 10.^((h_f * h_mis * h_l)/10);     
            end

            y_pilot = complex(zeros(Nr, N));
            for symbol = 1:N
                tx_index = mask(symbol);                           
                h_tx = h(:, tx_index, symbol);      
                n_tilda = sqrt((k_sq*P)/2) * (randn + 1i*randn);  

                y_pilot(:,symbol) = h_tx * (sp + n_tilda) + sqrt(N0 / 2) .* (randn(Nr, 1) + 1i * randn(Nr, 1));
            end

            % Channel Estimation and NMSE
            for symbol = 1:N
                tx_index = mask(symbol);
                h_true = h(:, tx_index, symbol);

                % LS and MLS Channel Estimator
                H_est_LS  = y_pilot(:, symbol) / sp;
                H_est_MLS = (conj(sp) .* y_pilot(:, symbol)) / ((k_sq + 1) * P);

                % per-symbol NMSE contribution
                normalizer = mean(abs(h_true).^2);
                NMSE_LS(SNR_i,  k_i) = NMSE_LS(SNR_i,  k_i) + mean(abs(H_est_LS  - h_true).^2) / normalizer / N;
                NMSE_MLS(SNR_i, k_i) = NMSE_MLS(SNR_i, k_i) + mean(abs(H_est_MLS - h_true).^2) / normalizer / N;
            end
        end
    end
end

% CSI = 1 -> Perfect CSI ///=========\\\\ CSI = 2 -> Actual CSI
function [BER] = Simulation_BER(SNR, k_sq, h_mis, h_l, Nt, Nr, N, M, alpha, mu, omega, CSI)
    BER = zeros(1,length(SNR));

    AntennaBits_Number = log2(Nt);
    Bits_Spatial = randi([0,1], N, AntennaBits_Number);
    Bits_Spatial_Mask = bi2de(Bits_Spatial, 'left-msb') + 1;
    
    % QPSK bits
    APMBits_Number = log2(M);
    Bits_Data = randi([0,1], N, APMBits_Number);
    I = 1 - 2*Bits_Data(:,2);   % bit2 -> I
    Q = 1 - 2*Bits_Data(:,1);   % bit1 -> Q
    x = (I + 1i*Q)/sqrt(2);     % 00->++, 01->-+, 11->--, 10->+-  
    voltages = [1+1i; -1+1i; -1-1i; 1-1i] / sqrt(2);
    P = mean(abs(x).^2);
    Bits_QPSKSymbols = [0 0; 0 1; 1 1; 1 0];
    
    % Channel Matrix
    h = complex(zeros(Nr, Nt, N));
    H = complex(zeros(Nr, Nt, N));
    
    % Pilot Symbol
    sp = (1 + 1i)/sqrt(2);  
    P_pilot = abs(sp)^2;  

    for i = 1:length(SNR)
        gamma = 10.^(SNR(i)/10);
        N0  = P / gamma;
        y = complex(zeros(Nr, N));
        y_pilot = complex(zeros(Nr, Nt, N));
    
        for i_tx = 1:Nt
            h_f = AlphaMuGenerator(alpha, mu, omega, Nr, N);  
            h(:, i_tx, :) = h_f * h_mis * h_l;    
            H(:, i_tx, :) = abs(h(:, i_tx, :)).^2;
        end
    
        for symbol = 1:N
            Index_Tx = Bits_Spatial_Mask(symbol);
            h_tx = h(:, Index_Tx, symbol); 
            H_tx = abs(h_tx).^2;
    
            h_tx_pilot = h(:, :, symbol); 
            H_tx_pilot = abs(h_tx_pilot).^2;
    
            xs = x(symbol);
    
            n_t = sqrt(k_sq * P/2) * (randn(1,1) + 1i * randn(1,1));  % Transmitter noise
            n_r = sqrt(k_sq * P .* H_tx/2) .* (randn(Nr,1) + 1i * randn(Nr,1));  % Receiver Noise for information signal
            n_r_pilot = sqrt(k_sq * P .* H_tx_pilot/2) .* (randn(Nr,1) + 1i * randn(Nr,1)); % Receiver Noise for plot signal
    
            y_pilot(:,:, symbol) = h_tx_pilot*sp + h_tx_pilot*n_t + n_r_pilot + sqrt(N0/2) .* (randn(Nr, Nt, 1) + 1i * randn(Nr, Nt, 1));
            y(:, symbol) = h_tx*xs + h_tx*n_t + n_r + sqrt(N0/2) * (randn(Nr,1) + 1i*randn(Nr,1)); 
            
        end
    
        Metric = zeros(Nt, length(voltages));
        C = zeros(Nr, Nr);
        SpatialBits_Decoded = zeros(N, AntennaBits_Number);
        QPSKBits_Decoded = zeros(N, APMBits_Number);
    
        for symbol = 1:N
            Min_Value = inf;
    
            if (CSI == 1) % Perfect CSI
                for i_tx = 1:Nt
                    h_tx = h(:,i_tx, symbol);
                    
                    h_x_hherm = h_tx * h_tx'; % Nr x Nr
        
                    Hi = diag(abs(h_tx).^2); % I_Nr
        
                    Ci = (k_sq*(h_x_hherm) + k_sq*Hi) * P + N0*eye(Nr);
        
                    E = y(:,symbol) - h_tx * (voltages.'); %E = (y-h_i*x)
                    Z = Ci \ E; % Z = (y-h_i*x)^H * C_i^(-1)  
                    Metric_FirstTerm = real(log(det(Ci)));
                    Metric_SecondTerm = real(sum(conj(E).*Z, 1));      % Second Term = Z*E         
        
                    Metric = Metric_FirstTerm + Metric_SecondTerm; % Metric results for QPSK symbol                  
        
                    % Save the indexes of current Tx if needed
                    [temp, index] = min(Metric);
                    if temp < Min_Value
                        MinIndex_Voltage = index;
                        MinIndex_Tx = i_tx;
                        Min_Value = temp;
                    end
                end
            elseif (CSI == 2) % Actual CSI
                H_est_MLS = (conj(sp) .* y_pilot(:, :, symbol)) / ((k_sq+1)*P_pilot);
                for i_tx = 1:Nt
                    h_tx = H_est_MLS(:, i_tx);
    
                    h_x_hherm = h_tx * h_tx'; % Nr x Nr
        
                    Hi = diag(abs(h_tx).^2); % I_Nr
        
                    Ci = (k_sq*(h_x_hherm) + k_sq*Hi) * P + N0*eye(Nr);
        
                    E = y(:,symbol) - h_tx * (voltages.'); %E = (y-h_i*x)
                    Z = Ci \ E; % Z = (y-h_i*x)^H * C_i^(-1)  
                    Metric_FirstTerm = real(log(det(Ci)));
                    Metric_SecondTerm = real(sum(conj(E).*Z, 1));      % Second Term = Z*E         
        
                    Metric = Metric_FirstTerm + Metric_SecondTerm; % Metric results for QPSK symbol                  
        
                    % Save the indexes of current Tx if needed
                    [temp, index] = min(Metric);
                    if temp < Min_Value
                        MinIndex_Voltage = index;
                        MinIndex_Tx = i_tx;
                        Min_Value = temp;
                    end
                end              
            end
            
            if (MinIndex_Voltage == 1) 
                QPSKBits_Decoded(symbol,:) = [0 0];
            elseif (MinIndex_Voltage == 2)
                QPSKBits_Decoded(symbol,:) = [0 1];
            elseif (MinIndex_Voltage == 3)
                QPSKBits_Decoded(symbol,:) = [1 1];
            elseif (MinIndex_Voltage == 4)
                QPSKBits_Decoded(symbol,:) = [1 0];
            end
            SpatialBits_Decoded(symbol,:) = de2bi(MinIndex_Tx-1, APMBits_Number, 'left-msb');
            QPSKBits_Decoded(symbol,:) = Bits_QPSKSymbols(MinIndex_Voltage,:);
        end
    
        Error_Count = sum(QPSKBits_Decoded ~= Bits_Data, 'all') + sum(Bits_Spatial ~= SpatialBits_Decoded, 'all');
        BER(i) = Error_Count / (4 * N);
    end


end

function [alpha_mu_sample,theta] = AlphaMuGenerator(alpha, mu, omega, row, N)
    m=mu/2;
    sn=[-1 1];

    H = zeros(row,N);
    theta = zeros(row,N);
    for i = 1:row
        X = sn((rand(1,N)>0.5)+1).*sqrt(gamrnd(m,omega,1,N)/mu);
        Y = sn((rand(1,N)>0.5)+1).*sqrt(gamrnd(m,omega,1,N)/mu);
        
        H(i,:) = X + 1i*Y ;
        theta(i,:) = atan2(Y,X);
    
    end
    
    alpha_mu_sample = abs(H).^(2/alpha).*exp(1i*theta);
end
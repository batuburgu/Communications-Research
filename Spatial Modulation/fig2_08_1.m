SNR = 0:1:25; % Signal to Noise Ratio (SNR)
N = 1e6; % # of symbols
A = sqrt(2); % Maximum aplitude of the cos wave being sent
Nt = 4; % # of Tx Antennas
Nr = 4; % # of Rx Antennas
M = 2; % M-ary Modulation

ber_conv_suboptimal = SM_SubOptimalDetection_Simulation(SNR, N, A, Nt, Nr, M, 1);
ber_const_suboptimal = SM_SubOptimalDetection_Simulation(SNR, N, A, Nt, Nr, M, 2);
ber_conv_optimal = SM_OptimalDetection_Simulation(SNR, N, A, Nt, Nr, M, 1);
ber_const_optimal = SM_OptimalDetection_Simulation(SNR, N, A, Nt, Nr, M, 2);
analytical_result = SM_OptimalDetection_Analytical(SNR, Nt, Nr);

figure;
semilogy(SNR, ber_conv_suboptimal, '-d', 'LineWidth',1.2,'MarkerSize',6); hold on;
semilogy(SNR, ber_const_suboptimal, '--d', 'LineWidth',1.2,'MarkerSize',6);
semilogy(SNR, ber_conv_optimal, '-*', 'LineWidth',1.2,'MarkerSize',6);
semilogy(SNR, ber_const_optimal, '--*', 'LineWidth',1.2,'MarkerSize',6);
semilogy(SNR, analytical_result, '--','LineWidth',1.5);
ylim([1e-5 1e0]);
xlim([0 25]);
grid on;
title('BER Performance of Spatial Modulation versus SNR, for m = 3, Nt = 4 and Nr = 4');
legend('SM (Mesleh), Conventional H', ...
       'SM (Mesleh), Constrained H', ...
       'SM (Optimal), Conventional H', ...
       'SM (Optimal), Constrained H', ...
       'SM (Optimal), Analytical Bound', ...
       'Location','southwest');
xlabel('\rho (dB)');
ylabel('Bit Error Rate (BER)');

function [ber] =  SM_SubOptimalDetection_Simulation(SNR, N, A, Nt, Nr, M, ChannelSelection)
    AbsoluteAmplitude = A;
    ber = zeros(1, length(SNR));

    Eb = (AbsoluteAmplitude^2)/2; 
    Es = Eb * log2(M);

    rho = 10.^(SNR / 10);

    voltages = zeros(1,M);
    for SNR_i = 1:M
        theta = (SNR_i-1) *  180;
        voltages(SNR_i) = -Es * cosd(theta);
    end
   
    % Bits to Send with Spatial Modulation
    bits = randi([0, 1], N, 3);

    % BPSK Coded Third Bit Voltages
    x_l = ones(1,N) * voltages(1);
    x_l(bits(:,3) > 0) = voltages(2);
    
    for SNR_i = 1:length(SNR)
        H = zeros(Nr, Nt, N) + 1i * zeros(Nr, Nt, N); 
        y = zeros(Nr,N);
            
        if (ChannelSelection == 1)
            % Conventional H
            for i_tx = 1:Nt
                h = (1/sqrt(2)) * (randn(Nr, N) + 1i * randn(Nr, N));  % Column of H (size: Nrx1)
       
                H(:, i_tx, :) = h;  
            end
        
        elseif (ChannelSelection == 2)
            % Constrained H
            for i_tx = 1:Nt
                h = (1/sqrt(2)) * (randn(Nr, N) + 1i * randn(Nr, N));  % Column of H (size: Nrx1)
                h_norm = zeros(1,N);
                for symbol = 1:N
                    h_norm(1, symbol) = norm(h(:,symbol), 'fro');
                end
    
                H(:, i_tx, :) = h ./ h_norm;  
            end
        end
        
        % Tx to 4 Rx Channel Matrices
        h_tx1 = reshape(H(:,1,:),Nr,N);
        h_tx2 = reshape(H(:,2,:),Nr,N);
        h_tx3 = reshape(H(:,3,:),Nr,N);
        h_tx4 = reshape(H(:,4,:),Nr,N);

        % Transmission
        for symbol = 1:N
            % Tx 1 - 00
            if (bits(symbol, 1) == 0 && bits(symbol, 2) == 0)
                y(:, symbol) = sqrt(rho(SNR_i)) * h_tx1(:, symbol) .* x_l(symbol) + sqrt(1/2) .* (randn(4, 1) + 1i * randn(4, 1));
            
            % Tx 2 - 01
            elseif (bits(symbol, 1) == 0 && bits(symbol, 2) == 1)
                y(:, symbol) = sqrt(rho(SNR_i)) * h_tx2(:, symbol) .* x_l(symbol) + sqrt(1/2) .* (randn(4, 1) + 1i * randn(4, 1));
            
            % Tx 3 - 11
            elseif (bits(symbol, 1) == 1 && bits(symbol, 2) == 1)
                y(:, symbol) = sqrt(rho(SNR_i)) * h_tx3(:, symbol) .* x_l(symbol) + sqrt(1/2) .* (randn(4, 1) + 1i * randn(4, 1));
            
            % Tx 4 - 10
            elseif (bits(symbol, 1) == 1 && bits(symbol, 2) == 0)
                y(:, symbol) = sqrt(rho(SNR_i)) * h_tx4(:, symbol) .* x_l(symbol) + sqrt(1/2) .* (randn(4, 1) + 1i * randn(4, 1));
            end
        end

        % Antenna Number Estimation and Spatial Demodulation
        estimated_bits = zeros(N, 3);
        g = zeros(Nt, N);
        alpha_opt = 0;
        
        for i_tx = 1:Nt
            h_j = squeeze(H(:, i_tx, :));  
            g(i_tx, :) = sum(conj(h_j) .* y, 1);
        end

        for symbol = 1:N
            [~, l_hat] = max(g(:, symbol));
            
            % Tx 1 - 00
            if (l_hat == 1)
                estimated_bits(symbol,1) = 0;
                estimated_bits(symbol,2) = 0;
            
            % Tx 2 - 01
            elseif (l_hat == 2)
                estimated_bits(symbol,1) = 0;
                estimated_bits(symbol,2) = 1;

            % Tx 3 - 11
            elseif (l_hat == 3)
                estimated_bits(symbol,1) = 1;
                estimated_bits(symbol,2) = 1;
            
            % Tx 4 - 10
            elseif (l_hat == 4)
                estimated_bits(symbol,1) = 1;
                estimated_bits(symbol,2) = 0;
            end

            g_selected = g(l_hat, symbol);  
            estimated_bits(symbol, 3) = real(g_selected) > alpha_opt;
        end
        
        num_errors = sum(bits ~= estimated_bits, 'all');
        ber(1, SNR_i) = num_errors / (3 * N);
        fprintf('SNR=%d dB -> BER=%.4e\n',SNR(SNR_i), ber(SNR_i));
    end   
end

function [ber] =  SM_OptimalDetection_Simulation(SNR, N, A, Nt, Nr, M, ChannelSelection)
    AbsoluteAmplitude = A;
    ber = zeros(1, length(SNR));

    Eb = (AbsoluteAmplitude^2)/2; 
    Es = Eb * log2(M);

    rho = 10.^(SNR / 10);

    voltages = zeros(1,M);
    for SNR_i = 1:M
        theta = (SNR_i-1) *  180;
        voltages(SNR_i) = -Es * cosd(theta);
    end
   
    % Bits to Send with Spatial Modulation
    bits = randi([0, 1], N, 3);

    % BPSK Coded Third Bit Voltages
    x_l = ones(1,N) * voltages(1);
    x_l(bits(:,3) > 0) = voltages(2);
    
    for SNR_i = 1:length(SNR)
        H = zeros(Nr, Nt, N) + 1i * zeros(Nr, Nt, N); 
        y = zeros(Nr,N);
            
        if (ChannelSelection == 1)
            % Conventional H
            for i_tx = 1:Nt
                h = (1/sqrt(2)) * (randn(Nr, N) + 1i * randn(Nr, N));  % Column of H (size: Nrx1)
       
                H(:, i_tx, :) = h;  
            end
        
        elseif (ChannelSelection == 2)
            % Constrained H
            for i_tx = 1:Nt
                h = (1/sqrt(2)) * (randn(Nr, N) + 1i * randn(Nr, N));  % Column of H (size: Nrx1)
                h_norm = zeros(1,N);
                for symbol = 1:N
                    h_norm(1, symbol) = norm(h(:,symbol), 'fro');
                end
    
                H(:, i_tx, :) = h ./ h_norm;  
            end
        end
        
        % Tx to 4 Rx Channel Matrices
        h_tx1 = reshape(H(:,1,:),Nr,N);
        h_tx2 = reshape(H(:,2,:),Nr,N);
        h_tx3 = reshape(H(:,3,:),Nr,N);
        h_tx4 = reshape(H(:,4,:),Nr,N);

        % Transmission
        for symbol = 1:N
            % Tx 1 - 00
            if (bits(symbol, 1) == 0 && bits(symbol, 2) == 0)
                y(:, symbol) = sqrt(rho(SNR_i)) * h_tx1(:, symbol) .* x_l(symbol) + sqrt(1 / 2) .* (randn(4, 1) + 1i * randn(4, 1));
            
            % Tx 2 - 01
            elseif (bits(symbol, 1) == 0 && bits(symbol, 2) == 1)
                y(:, symbol) = sqrt(rho(SNR_i)) * h_tx2(:, symbol) .* x_l(symbol) + sqrt(1 / 2) .* (randn(4, 1) + 1i * randn(4, 1));
            
            % Tx 3 - 11
            elseif (bits(symbol, 1) == 1 && bits(symbol, 2) == 1)
                y(:, symbol) = sqrt(rho(SNR_i)) * h_tx3(:, symbol) .* x_l(symbol) + sqrt(1 / 2) .* (randn(4, 1) + 1i * randn(4, 1));
            
            % Tx 4 - 10
            elseif (bits(symbol, 1) == 1 && bits(symbol, 2) == 0)
                y(:, symbol) = sqrt(rho(SNR_i)) * h_tx4(:, symbol) .* x_l(symbol) + sqrt(1 / 2) .* (randn(4, 1) + 1i * randn(4, 1));
            end
        end

        % Antenna Number Estimation and Spatial Demodulation
        alpha_opt = 0;
        estimated_bits = zeros(N, 3);  
        j_estimated = zeros(1, N);  
        q_estimated = zeros(1, N);  

        for symbol = 1:N
            min_val = inf; 
            for i_tx = 1:Nt
                h_j = squeeze(H(:, i_tx, symbol));  
                
                for i_voltage = 1:length(voltages)
                    g = h_j .* voltages(i_voltage);  
                    right_side = sum(conj(y(:, symbol)) .* g, 1);  
                    
                    g_norm = norm(g, 'fro');  
                    
                    metric = sqrt(rho(SNR_i))*g_norm^2 - 2 * real(right_side); 
        
                    if metric < min_val
                        min_val = metric;
                        j_estimated(symbol) = i_tx;  
                        q_estimated(symbol) = voltages(i_voltage) > alpha_opt; 
                    end
                end
            end
            
            if (j_estimated(symbol) == 1)
                estimated_bits(symbol, 1) = 0;
                estimated_bits(symbol, 2) = 0;
            
            elseif (j_estimated(symbol) == 2)
                estimated_bits(symbol, 1) = 0;
                estimated_bits(symbol, 2) = 1;
        
            elseif (j_estimated(symbol) == 3)
                estimated_bits(symbol, 1) = 1;
                estimated_bits(symbol, 2) = 1;
        
            elseif (j_estimated(symbol) == 4)
                estimated_bits(symbol, 1) = 1;
                estimated_bits(symbol, 2) = 0;
            end
            
            estimated_bits(symbol, 3) = q_estimated(symbol);
        end   
        num_errors = sum(bits ~= estimated_bits, 'all'); 
        ber(1, SNR_i) = num_errors / (3 * N);   
        fprintf('SNR=%d dB -> BER=%.4e\n',SNR(SNR_i), ber(SNR_i));
    end  
end

function [analytical_error] =  SM_OptimalDetection_Analytical(SNR, Nt, Nr)
    analytical_error = zeros(1,length(SNR));
    
    for SNR_i = 1:length(SNR)
        rho = 10^(SNR(SNR_i)/10);   
        
        mu_alpha = 0.5 * (1 - sqrt((rho/2)/(1 + rho/2)));
        
        sum_term = 0;
        for k = 0:Nr-1
            sum_term = sum_term + nchoosek(Nr-1+k, k) * (1 - mu_alpha)^k;
        end
        
        analytical_error(SNR_i) = Nt * (mu_alpha^Nr) * sum_term;
    end
end


SNR = 0:1:30; % Signal to Noise Ratio (SNR)
N = 1e6; % # of symbols
A = sqrt(2); % Maximum aplitude of the cos wave being sent
Nt = 8; % # of Tx Antennas
Nr_Arr = [4 2 1]; % # of Rx Antennas
M = 8; % M-ary Modulation

ssk_ber_simulation = SSK_Simulation(SNR, N, Nt, Nr_Arr, M);
[ssk_ber_analytical1, ssk_ber_analytical2] = SSK_Analytical(SNR, Nt, Nr_Arr(1));

figure();
semilogy(SNR, ssk_ber_simulation(1, :), '-*', 'LineWidth',1.2,'MarkerSize',6); hold on;
semilogy(SNR, ssk_ber_analytical1, '--', 'LineWidth',1.2,'MarkerSize',6);
semilogy(SNR, ssk_ber_analytical2, '--.', 'LineWidth',1.2,'MarkerSize',6);
ylim([1e-5 1e0]);
xlim([0 20]);
grid on;
title('Figure 4: BER Performance of SSK versus SNR, for m = 3, Nt = 8 and Nr = 4');
legend('SSK (Nt = 8), Simulation',...
    'SSK Analytical (Nt = 8), Eq. (8)', ...
    'SSK Analytical (Nt = 8), Eq. (9)',...
       'Location','southwest');
xlabel('\rho (dB)');
ylabel('Bit Error Rate (BER)');

figure();
semilogy(SNR, ssk_ber_simulation(1, :), '-^', 'LineWidth',1.2,'MarkerSize',6); hold on;
semilogy(SNR, ssk_ber_simulation(2, :), '-square', 'LineWidth',1.2,'MarkerSize',6);
semilogy(SNR, ssk_ber_simulation(3, :), '-o', 'LineWidth',1.2,'MarkerSize',6);
ylim([1e-5 1e0]);
xlim([0 30]);
grid on;
title('Figure 6: BER Performance of SSK versus SNR for varying Nr, for m = 3');
legend('SSK (Nt = 8, Nr = 4), Simulation',...
    'SSK (Nt = 8, Nr = 2), Simulation', ...
    'SSK (Nt = 8, Nr = 1), Simulation',...
       'Location','southwest');
xlabel('\rho (dB)');
ylabel('Bit Error Rate (BER)');


function [ber] =  SSK_Simulation(SNR, N, Nt, Nr_Arr, M)
    ber = zeros(length(Nr_Arr), length(SNR));
    apm_bits = log2(M);
    antenna_bits = log2(Nt);

    rho = 10.^(SNR / 10);

    % Bits to Send with Spatial Modulation
    bits = randi([0, 1], N, apm_bits + antenna_bits);
    for Nr_i = 1: length(Nr_Arr)
        Nr = Nr_Arr(Nr_i);
        for SNR_i = 1:length(SNR)
            % Channel Matrix
            H = zeros(Nr, Nt, N) + 1i * zeros(Nr, Nt, N); 
            for i_tx = 1:Nt
                h = (1/sqrt(2)) * (randn(Nr, N) + 1i * randn(Nr, N));  % Column of H (size: Nrx1)
                H(:, i_tx, :) = h;  
            end
           
            % Tx to 4 Rx Channel Matrices
            h_tx1 = reshape(H(:,1,:),Nr,N);
            h_tx2 = reshape(H(:,2,:),Nr,N);
            h_tx3 = reshape(H(:,3,:),Nr,N);
            h_tx4 = reshape(H(:,4,:),Nr,N);
            h_tx5 = reshape(H(:,5,:),Nr,N);
            h_tx6 = reshape(H(:,6,:),Nr,N);
            h_tx7 = reshape(H(:,7,:),Nr,N);
            h_tx8 = reshape(H(:,8,:),Nr,N);
           
            % Transmission
            y = zeros(Nr, N);
    
            for symbol = 1:N
                % Tx 1 - 000
                mask = bi2de(bits(symbol,1:antenna_bits), 'left-msb');
                if (mask == 0)
                    y(:, symbol) = sqrt(rho(SNR_i)) * h_tx1(:, symbol) + sqrt(1/2) .* (randn(Nr, 1) + 1i * randn(Nr, 1));
                
                % Tx 2 - 001
                elseif (mask == 1)
                    y(:, symbol) = sqrt(rho(SNR_i)) * h_tx2(:, symbol) + sqrt(1/2) .* (randn(Nr, 1) + 1i * randn(Nr, 1));
                
                % Tx 3 - 010
                elseif (mask == 2)
                    y(:, symbol) = sqrt(rho(SNR_i)) * h_tx3(:, symbol) + sqrt(1/2) .* (randn(Nr, 1) + 1i * randn(Nr, 1));
                
                % Tx 4 - 011
                elseif (mask == 3)
                    y(:, symbol) = sqrt(rho(SNR_i)) * h_tx4(:, symbol) + sqrt(1/2) .* (randn(Nr, 1) + 1i * randn(Nr, 1));
                
                % Tx 5 - 100
                elseif (mask == 4)
                    y(:, symbol) = sqrt(rho(SNR_i)) * h_tx5(:, symbol) + sqrt(1/2) .* (randn(Nr, 1) + 1i * randn(Nr, 1));
                
                % Tx 6 - 101
                elseif (mask == 5)
                    y(:, symbol) = sqrt(rho(SNR_i)) * h_tx6(:, symbol) + sqrt(1/2) .* (randn(Nr, 1) + 1i * randn(Nr, 1));
                
                % Tx 7 - 110
                elseif (mask == 6)
                    y(:, symbol) = sqrt(rho(SNR_i)) * h_tx7(:, symbol) + sqrt(1/2) .* (randn(Nr, 1) + 1i * randn(Nr, 1));
                
                % Tx 8 - 111
                elseif (mask == 7)
                    y(:, symbol) = sqrt(rho(SNR_i)) * h_tx8(:, symbol) + sqrt(1/2) .* (randn(Nr, 1) + 1i * randn(Nr, 1));
                end
            end
    
            % Antenna Number Estimation and Spatial Demodulation
            estimated_bits = zeros(N, 3);
            metric = zeros(Nt, N);
            
            for i_tx = 1:Nt
                h_j = reshape(H(:,i_tx,:),Nr,N);
                metric(i_tx, :) = real(sum(conj(y - 0.5*sqrt(rho(SNR_i))*h_j) .* h_j, 1));
            end
    
            for symbol = 1:N
                [~, j_hat] = max(metric(:, symbol));
                
                % Tx 1 - 000
                if (j_hat == 1)
                    estimated_bits(symbol,1) = 0;
                    estimated_bits(symbol,2) = 0;
                    estimated_bits(symbol,3) = 0;
                
                % Tx 2 - 001
                elseif (j_hat == 2)
                    estimated_bits(symbol,1) = 0;
                    estimated_bits(symbol,2) = 0;
                    estimated_bits(symbol,3) = 1;
    
                % Tx 3 - 010
                elseif (j_hat == 3)
                    estimated_bits(symbol,1) = 0;
                    estimated_bits(symbol,2) = 1;
                    estimated_bits(symbol,3) = 0;
                
                % Tx 4 - 011
                elseif (j_hat == 4)
                    estimated_bits(symbol,1) = 0;
                    estimated_bits(symbol,2) = 1;
                    estimated_bits(symbol,3) = 1;
    
                % Tx 5 - 100
                elseif (j_hat == 5)
                    estimated_bits(symbol,1) = 1;
                    estimated_bits(symbol,2) = 0;
                    estimated_bits(symbol,3) = 0;
                
                % Tx 6 - 101
                elseif (j_hat == 6)
                    estimated_bits(symbol,1) = 1;
                    estimated_bits(symbol,2) = 0;
                    estimated_bits(symbol,3) = 1;
    
                % Tx 7 - 110
                elseif (j_hat == 7)
                    estimated_bits(symbol,1) = 1;
                    estimated_bits(symbol,2) = 1;
                    estimated_bits(symbol,3) = 0;
                
                % Tx 8 - 111
                elseif (j_hat == 8)
                    estimated_bits(symbol,1) = 1;
                    estimated_bits(symbol,2) = 1;
                    estimated_bits(symbol,3) = 1; 
                end
            end
            
            num_errors = sum(bits(:, 1:3) ~= estimated_bits, 'all');
            ber(Nr_i,SNR_i) = num_errors / (3 * N);
            fprintf('SNR=%d dB -> BER=%.4e\n',SNR(SNR_i), ber(Nr_i, SNR_i));
        end   
    
    end
     
end

function [ssk_eq8, ssk_eq9] =  SSK_Analytical(SNR, Nt, Nr)
    ssk_eq8 = zeros(1,length(SNR));
    ssk_eq9 = zeros(1,length(SNR)); 
    n_neigh = 2;
    N_sigma = 0;
    im_bits = log2(Nt);

    natural_code = de2bi(0:2^im_bits-1, im_bits, 'left-msb'); % Antenna Bits
        
    for j = 1:Nt
        for j_hat = j+1:Nt
            num_errors = sum(natural_code(j,:) ~= natural_code(j_hat,:), 'all');
            N_sigma = N_sigma + 2*num_errors; 
        end
    end

    for SNR_i = 1:length(SNR)
        rho = 10^(SNR(SNR_i)/10);   
        sum_term = 0;

        gamma_alpha = 0.5 * (1 - sqrt((rho/2)/(1 + rho/2)));

        for k = 0:Nr-1
            sum_term = sum_term + nchoosek(Nr-1+k, k) * (1 - gamma_alpha)^k;
        end

        ssk_eq8(SNR_i) = N_sigma * (gamma_alpha^Nr) * sum_term / Nt;
        ssk_eq9(SNR_i) = n_neigh * (N_sigma / (Nt * (Nt -1))) * (gamma_alpha^Nr) * sum_term;
    end
end
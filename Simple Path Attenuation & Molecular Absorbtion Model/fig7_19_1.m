% alpha-mu distribution parameters
mu_array = [1,3,8];

% Outage Threshold
lower_gamma_th_array = [10^(1/10), 10^(15/10)]; 

arg = 0:2:100;

outage_probability = analytical_outage_probability(mu_array, lower_gamma_th_array);

figure;
hold on;

colors = [
    0.1216 0.4667 0.7059;  % blue
    1.0000 0.7647 0.0000;  % gold
    0.0000 0.5882 0.5333;  % cyan 
    0.8392 0.1529 0.1569;  % red
    0.5804 0.4039 0.7412;  % purple
    0.5490 0.3373 0.2941   % brown
];

% gamma_th = 1
semilogy(arg, squeeze(outage_probability(1,1,:)), '--square', 'Color', colors(1,:), 'LineWidth', 1.4, 'MarkerSize', 6);
semilogy(arg, squeeze(outage_probability(2,1,:)), '--o', 'Color', colors(2,:), 'LineWidth', 1.4, 'MarkerSize', 6);
semilogy(arg, squeeze(outage_probability(3,1,:)), '--*', 'Color', colors(3,:), 'LineWidth', 1.4, 'MarkerSize', 6);

% gamma_th = 15
semilogy(arg, squeeze(outage_probability(1,2,:)), '-square', 'Color', colors(4,:), 'LineWidth', 1.4, 'MarkerSize', 6);
semilogy(arg, squeeze(outage_probability(2,2,:)), '-o', 'Color', colors(5,:), 'LineWidth', 1.4, 'MarkerSize', 6);
semilogy(arg, squeeze(outage_probability(3,2,:)), '-*', 'Color', colors(6,:), 'LineWidth', 1.4, 'MarkerSize', 6);

% Axis and Labels
ylim([1e-5, 1])
xlim([0, 100])
set(gca, 'YScale', 'log')
set(gca, 'FontSize', 12)

xlabel('$\frac{P \cdot |h_l|^2}{N_0}$ [dB]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Outage Probability', 'FontSize', 14)
title('Figure 7: Outage Probability vs SNR for Different Values of \mu and \gamma_{th}, Ideal Front End and d = 30m', 'FontSize', 16)

legend({'$\mu = 1$, $\gamma_{th} = 1$', ...
        '$\mu = 3$, $\gamma_{th} = 1$', ...
        '$\mu = 8$, $\gamma_{th} = 1$', ...
        '$\mu = 1$, $\gamma_{th} = 15$', ...
        '$\mu = 3$, $\gamma_{th} = 15$', ...
        '$\mu = 8$, $\gamma_{th} = 15$'}, ...
        'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 12, 'Box', 'off');

grid on;

function [outage_probability] = analytical_outage_probability(mu, lower_gamma_th)
    alpha = 2;
    
    arg = 0:2:100;
    lin_arg = 10.^(arg/10);
    
    % Pointing Error Parameters
    lower_gamma = 1;
    A0 = 1;
    outage_probability = zeros(length(mu), length(lower_gamma_th), length(lin_arg));

    for i_mu = 1:length(mu)
        for i_lowergamma = 1:length(lower_gamma_th)
            for i_snrargument = 1:length(lin_arg)
    
                % Probability Functions
                h_hat_f = sqrt(mu(i_mu).^(2/alpha) .* gamma(mu(i_mu)) ./ gamma(mu(i_mu) + 2/alpha));
                f_hfp = @(x) lower_gamma^2 * A0^(-lower_gamma^2) * mu(i_mu)^(lower_gamma^2 / alpha) .* x.^(lower_gamma^2 - 1) .*...
                        igamma((alpha.*mu(i_mu)-lower_gamma.^2)./alpha, mu(i_mu) .* (x.^alpha) .* ...
                        (A0.^(-alpha)) / h_hat_f^alpha) ./ (h_hat_f.^alpha .* gamma(mu(i_mu))); 
    
    
                snr_argument = sqrt(lower_gamma_th(i_lowergamma)./lin_arg(i_snrargument));
                outage_probability(i_mu,i_lowergamma,i_snrargument) = integral(f_hfp, 1e-6, snr_argument);    
            end   
        end
    end
end


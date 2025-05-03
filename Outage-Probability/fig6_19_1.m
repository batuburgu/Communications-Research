% System Parameters
c = 3*10^8; % Speed of Light (m/s)
Gt = 10^(55/10); % LinearTransmit Antenna Gain 
Gr = 10^(55/10); % Linear Receive Antenna Gain 
f = linspace(270e9, 400e9, 67);
d = 30; % Distance in meters

% Outage Threshold
threshold = [0.5 1 1.5 2]; % Linear

% P/N0
power_ratio_db = [10 20];

% alpha-mu Distribution Parameters
alpha = 2;
mu = 4;

% Pointing Error Parameters
lower_gamma = 3.6;
A0 = 1;

% Propagation Parameters
% g parameters
c1 = 10.835; % cm^-1
c2 = 12.664; % cm^-1
p1 = 5.54 * 10^-37; % Hz^-3
p2 = -3.94 * 10^-25; %Hz^-2
p3 = 9.06 * 10^-14; % Hz^-1
p4 = -6.36 * 10^-3; 

% A,B,C,D Parameters
g1 = 0.2205;
g2 = 0.1303;
g3 = 0.0294;
g4 = 0.4093;
g5 = 0.0925;
g6 = 2.014;
g7 = 0.1702;
g8 = 0.0303;
g9 = 0.537;
g10 = 0.0956;

% v parameters
phi = 0.5; 
p = 101325; %Pa

% p_w parameters
q1 = 6.1121;
q2 = 1.0007;
q3 = 3.46 * 10^-6;
q4 = 17.502;
q5 = 240.97;
T = 296;

% Molecular Absorbtion Parameters
p_w = q1 * (q2 + q3 * p) *  exp(q4*T/(q5 + T)); % Buck eq. 1
v = (phi / 100) * p_w/p;
A = g1 * v * (g2*v + g3);
B = (g4*v + g5)^2;
C = g6*v * (g7*v + g8);
D = (g9*v + g10)^2;

y1 = A ./ (B + (f/(100 * c) - c1).^2);
y2 = C ./ (D + (f/(100 * c) - c2).^2);
g = p1.*(f.^3) + p2.*(f.^2) + p3.*f + p4;
k_a = (y1 + y2) + g; % Absorbtion Coefficient 

analytical_result = outage_prob_analytical(threshold, power_ratio_db, c, Gt, Gr, f, d, alpha, mu, lower_gamma, A0, k_a);

% Plot
figure;
hold on;
colors = lines(length(threshold)); 
markers = {'square', 'o', 'pentagram', 'hexagram'};    
for pri = 1:length(power_ratio_db)
    for ti = 1:length(threshold)
        y = analytical_result(pri, ti, :);
        y = y(:);
        
        if pri == 1
            semilogy(f / 1e9, y, ...
                     'LineWidth', 1.5, ...
                     'Marker', markers{ti}, ...
                     'Color', colors(ti, :), ...
                     'DisplayName', ['\gamma_{th} = ' num2str(threshold(ti))]);
        else
            semilogy(f / 1e9, y, ...
                     'LineWidth', 1.5, ...
                     'Marker', markers{ti}, ...
                     'Color', colors(ti, :), ...
                     'HandleVisibility', 'off'); 
        end
    end
end
xlim([270, 400]);
ylim([1e-7, 1]);
set(gca, 'YScale', 'log')
xlabel('f (GHz)', 'FontSize', 12);
ylabel('Outage Probability', 'FontSize', 12);
title('Fig.6: Outage Probability vs Frequency for Different Threshold Values', 'FontSize', 14);
legend('Location', 'southeast');
grid on;

function [outage_probability] = outage_prob_analytical(threshold, power_ratio_db, c, Gt, Gr, f, d, alpha, mu, lower_gamma, A0, k_a)
    h_hat_f = sqrt(mu^(2/alpha) * gamma(mu) / gamma(mu + 2/alpha));

    outage_probability = zeros(length(power_ratio_db), length(threshold), length(f));
    
    for pri = 1:length(power_ratio_db)
        P_lin = 10^(power_ratio_db(pri)/10);

        for ti = 1:length(threshold)
            for fi = 1:length(f)
                
                h_fl = (c./(4*pi*d.*f(fi))).^2 * Gt * Gr; % Propagation Gain
                h_al = exp(-(k_a * d));        % Absorption Gain
                h_l = h_fl * h_al;
            
                f_hfp = @(x) lower_gamma^2 * A0^(-lower_gamma^2) * mu^(lower_gamma^2 / alpha) .* x.^(lower_gamma^2 - 1) .*...
                            igamma((alpha.*mu-lower_gamma.^2)./alpha, mu .* (x.^alpha) .* ...
                            (A0.^(-alpha)) / h_hat_f^alpha) ./ (h_hat_f.^alpha .* gamma(mu)); 
                
                snr = sqrt(threshold(ti) ./ (abs(h_l).^2 * P_lin));
            
                outage_probability(pri,ti, fi) = integral(f_hfp, 0, snr(fi));
            end  
        end
    end
end










 
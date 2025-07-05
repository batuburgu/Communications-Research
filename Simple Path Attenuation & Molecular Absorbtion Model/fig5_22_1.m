c = 3*10^8; % Speed of Light (m/s)
Gt = 10^(55/10); % Linear Transmit Antenna Gain 
Gr = 10^(55/10); % Linear Receive Antenna Gain 
f = 300e9;
d = 50; % m
N = 1e6;

phi = [2.0437, 8.1748];
alpha = 2;
mu = 1.5;
omega = 1;
S0 = 0.39;

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
ro = 0.5; 
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
v = (ro / 100) * p_w/p;
A = g1 * v * (g2*v + g3);
B = (g4*v + g5)^2;
C = g6*v * (g7*v + g8);
D = (g9*v + g10)^2;

y1 = A ./ (B + (f/(100 * c) - c1).^2);
y2 = C ./ (D + (f/(100 * c) - c2).^2);
g = p1.*(f.^3) + p2.*(f.^2) + p3.*f + p4;
k_a = (y1 + y2) + g; % Absorbtion Coefficient 

gamma_bar_db = linspace(5, 50, 181); % Average SNR in dB
gamma_bar = 10.^(gamma_bar_db/10); % Average SNR Linear
gamma_th = 10^(2/10); % Linear =  2dB

h_fl = c./(4*pi*d*f) * sqrt(Gt*Gr); % Propagation Gain
h_al = exp(-(k_a * d)/2);        % Absorption Gain
h_l = h_fl * h_al;

h_hat_f = sqrt(mu^(2/alpha) * gamma(mu) * omega / gamma(mu + 2/alpha));

analytical_op = op_analytical(gamma_bar, gamma_th, phi, alpha, mu, h_l, h_hat_f, S0);
simulation_op = op_simulation(gamma_bar, gamma_th, phi, alpha, mu, h_l, omega, S0, N);

% Plot
figure;
semilogy(gamma_bar_db, analytical_op(1,:), '-*', 'LineWidth', 1.4, 'MarkerSize', 5); hold on;
semilogy(gamma_bar_db, analytical_op(2,:), '-*', 'LineWidth', 1.4, 'MarkerSize', 5);
semilogy(gamma_bar_db, simulation_op(1,:), '-o', 'LineWidth', 1.4, 'MarkerSize', 5);
semilogy(gamma_bar_db, simulation_op(2,:), '-o', 'LineWidth', 1.4, 'MarkerSize', 5);

xlim([5, 50]);
ylim([1e-6, 1]);
set(gca, 'YScale', 'log');

grid on;
box on;
set(gca, 'FontSize', 12); 

legend({'$\phi = 2.0437$ (Analytical)', '$\phi = 8.1748$ (Analytical)', ...
        '$\phi = 2.0437$ (Simulation)', '$\phi = 8.1748$ (Simulation)'}, ...
        'Interpreter','latex', 'FontSize', 11, 'Location', 'southwest');

xlabel('$\bar{\gamma}$ [dB]', 'Interpreter','latex', 'FontSize', 14);
ylabel('Outage Probability', 'Interpreter','latex', 'FontSize', 14);
title('Figure 5: Outage Probability for Single THz Links ($d = 50$ m)', ...
    'Interpreter','latex', 'FontSize', 15);

function [outage_probability] = op_analytical(gamma_bar, gamma_th, phi, alpha, mu, h_l, h_hat_f, S0)
    outage_probability = zeros(length(phi), length(gamma_bar));
    for phi_i = 1:length(phi)
        A = (phi(phi_i) * mu^(phi(phi_i)/alpha) * h_l ^(-phi(phi_i))) / (2*h_hat_f^(alpha) * S0^(phi(phi_i)) *gamma(mu));
        B = mu ./ (h_hat_f .* h_l * S0)^alpha;

        F_gamma = @(gamma) 2 * A .* gamma_bar.^(-phi(phi_i)/2) .* gamma.^(phi(phi_i)/2) .* meijerG(1-(phi(phi_i)/alpha), 1, [0, ((alpha*mu - phi(phi_i))/alpha)], (-phi(phi_i)/alpha), B.*(gamma./gamma_bar).^(alpha/2)) / alpha;

        outage_probability(phi_i,:) = F_gamma(gamma_th);
    end
    
end

function [op_sim] = op_simulation(gamma_bar, gamma_th, phi, alpha, mu, h_l, omega, S0, N)
    op_sim = ones(length(phi), length(gamma_bar));
    U = rand(1, N); % Uniform Distribution

    for phi_i = 1:length(phi)
        h_mis = U.^(1./phi(phi_i));

        for gammai = 1:length(gamma_bar)
            [h_f, ~] = AlphaMuGenerator(alpha, mu, omega, N); 
            channel_power = abs(h_l .* h_mis .* h_f * S0).^2;
        
            received_snr = gamma_bar(gammai) * channel_power;
            op_sim(phi_i, gammai) = sum(received_snr <= gamma_th) / length(received_snr);
        end
    end
end

function [alpha_mu_sample,theta] = AlphaMuGenerator(alpha, mu, omega, N)
    m=mu/2;
    sn=[-1 1];

    X = sn((rand(1,N)>0.5)+1).*sqrt(gamrnd(m,omega,1,N)/mu);
    Y = sn((rand(1,N)>0.5)+1).*sqrt(gamrnd(m,omega,1,N)/mu);
    
    H = X + 1i*Y ;
    theta = atan2(Y,X);
    
    alpha_mu_sample = abs(H).^(2/alpha).*exp(1i*theta);
end
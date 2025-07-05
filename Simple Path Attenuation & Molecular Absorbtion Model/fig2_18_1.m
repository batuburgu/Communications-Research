% System Parameters
c = 3*10^8; % Speed of Light (m/s)
Gt = 1; % Transmit Antenna Gain (linear)
Gr = 1; % Receive Antenna Gain (linear)
d = linspace(0.1, 100, 1001); % Distance in meters
power_coeff = [10, 50, 75, 100, 120]; % dB
power_coeff_lin = 10.^(power_coeff/10); %linear

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
phi = 0.23; 
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


BW = 125e6; % filter bandwidth
integral_result = zeros(1,length(d));
for di = 1:length(d)
    I_fb = @(f) exp(-d(di) .* (((A ./ (B + (f./(100 * c) - c1).^2)) + (C ./ (D + (f./(100 * c) - c2).^2)) + p1.*(f.^3) + p2.*(f.^2) + p3.*f + p4)))./ f.^2;
    integral_result(di) = integral(I_fb, 275e9 - BW/2, 275e9 + BW/2);
end

figure;
hold on;

colors = lines(length(power_coeff)); 
legendStrings = cell(1, length(power_coeff));

for pci = 1:length(power_coeff)
    snr = (c^2 ./ (4*pi*d).^2) .* (power_coeff_lin(pci) / BW) .* integral_result;
    snr_db = 10*log10(snr);
    semilogx(d, snr_db, 'LineWidth', 1.5, 'Color', colors(pci,:));
    
    legendStrings{pci} = sprintf('g = %d dB', power_coeff(pci));
end

set(gca, 'XScale', 'log')
set(gca, 'XTick', [0.1 1 10 100])                   
set(gca, 'XTickLabel', {'0.1', '1', '10', '100'})   
ylim([-120, 60])
xlim([0.1, 100])
title('Figure 2: SNR vs Distance')
xlabel('d (m)');
ylabel('SNR (dB)');
legend(legendStrings, 'Location', 'northeast');
grid on;
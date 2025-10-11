% System Parameters
c = 3*10^8; % Speed of Light (m/s)
Gt = 1; % Transmit Antenna Gain in mW
Gr = 1; % Receive Antenna Gain in mW
freq = linspace(270,400,1301); % Frequencies 270 GHz - 400 GHz
f = freq * 1e9;
d = [0.1, 1, 5, 10, 100, 500]; % Distance in meters

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
q1 = 6.1115;
q2 = 23.306;
q3 = 333.7;
q4 = 279.82;
T = 296;

% Molecular Absorbtion Parameters
p_w = q1 * exp(((q2 - T/q3)*T) / (q4 + T)); % Buck eq. 1
v = (phi / 100) * p_w/p;
A = g1 * v * (g2*v + g3);
B = (g4*v + g5)^2;
C = g6*v * (g7*v + g8);
D = (g9*v + g10)^2;

y1 = A ./ (B + (f/(100 * c) - c1).^2);
y2 = C ./ (D + (f/(100 * c) - c2).^2);
g = p1.*(f.^3) + p2.*(f.^2) + p3.*f + p4;
k_a = (y1 + y2) + g; % Absorbtion Coefficient 

figure;
hold on;

for di = 1:length(d)
    L = (c./(4*pi*d(di).*f)).^2 * Gt * Gr .* exp(-d(di).*k_a);   
    L_db = 10*log10(L);
    plot(freq, L_db, 'LineWidth', 2, 'DisplayName', ['d = ' num2str(d(di)) ' m']);
end

ylim([-350, -50])
xlim([270, 400])
title('Path Loss vs Frequency for Different Distances', 'FontSize', 14)
xlabel('Frequency (GHz)', 'FontSize', 12)
ylabel('Path Gain (dB)', 'FontSize', 12)
legend('Location', 'southwest')
grid on
box on
set(gca, 'FontSize', 11) 

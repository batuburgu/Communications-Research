% System Parameters
c = 3*10^8; % Speed of Light (m/s)
f = linspace(270e9,400e9,1201);
d = [100, 500, 1000];
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
phi = 1; 
p = 101325; %Pa

% p_w parameters
q1 = 6.1121;
q2 = 1.0007;
q3 = 3.46 * 10^-6;
q4 = 17.502;
q5 = 240.97;
T = 296;

% Molecular Absorbtion Parameters
p_w = q1 * (q2 + q3*p) *  exp(q4*T/(q5 + T)); % Buck eq. 1
v = (phi / 100) * p_w/p;
A = g1 * v * (g2*v + g3);
B = (g4*v + g5)^2;
C = g6*v * (g7*v + g8);
D = (g9*v + g10)^2;

y1 = A ./ (B + (f/(100 * c) - c1).^2); % (100*f/c)
y2 = C ./ (D + (f/(100 * c) - c2).^2); % (100*f/c)
g = p1.*(f.^3) + p2.*(f.^2) + p3.*f + p4;
k_a = y1 + y2 + g;

% Absorbtion Coefficients Plot
figure(1);
semilogy(f, k_a, 'LineWidth', 2); 
xlim([270e9, 400e9])
ylim([1e-4, 1])
xticks(1e9 * (270:10:400))
xticklabels(string(270:10:400))
title('Figure 1: Absorption Coefficient vs Frequency', 'FontSize', 14)
xlabel('Frequency (GHz)', 'FontSize', 12)
ylabel('Absorption Coefficient (cm^{-1})', 'FontSize', 12)
legend('d = 100m', 'FontSize', 10, 'Location', 'northeast')
grid on

% Absorbtion Loss plot
figure(2);
hold on;
pl_abs = zeros(length(d), length(f));
for di2 = 1:length(d)
    pl_abs(di2, :) = 10 * log10(exp(d(di2) * k_a));
    if di2 ~= 2
        plot(f, pl_abs(di2, :), 'LineWidth', 2); 
    end
end

ylim([-50, 120])
xlim([270e9, 400e9])
xticks(1e9 * (270:10:400))
xticklabels(string(270:10:400))
title('Figure 3: Absorption Loss vs Frequency', 'FontSize', 14)
xlabel('Frequency (GHz)', 'FontSize', 12)
ylabel('Absorption Loss (dB)', 'FontSize', 12)
legend('d = 100m', 'd = 1000m', 'FontSize', 10, 'Location', 'northeast')
grid on

figure(3);
hold on;
for di2 = 1:length(d)
    plot(f, pl_abs(di2, :), 'LineWidth', 2); 
end

ylim([0, 50])
xlim([320e9, 380e9])
xticks(1e9 * (270:10:400))
xticklabels(string(270:10:400))
title('Figure 5: Absorption Loss vs Frequency (Zoomed)', 'FontSize', 14)
xlabel('Frequency (GHz)', 'FontSize', 12)
ylabel('Absorption Loss (dB)', 'FontSize', 12)
legend('d = 100m', 'd = 500m', 'd = 1000m', 'FontSize', 10, 'Location', 'northeast')
grid on

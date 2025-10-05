% ITU-R P.676-13 -Table 1 - Spectroscopic Data for Oxygen Attenuation
table1_data=[50.474214,0.975,9.651,6.690,0.0,2.566,6.850,50.987745,2.529,8.653,7.170,0.0,2.246,6.800,51.503360,6.193,7.709,7.640,0.0,1.947,6.729,52.021429,14.320,6.819,8.110,0.0,1.667,6.640,52.542418,31.240,5.983,8.580,0.0,1.388,6.526,53.066934,64.290,5.201,9.060,0.0,1.349,6.206,53.595775,124.600,4.474,9.550,0.0,2.227,5.085,54.130025,227.300,3.800,9.960,0.0,3.170,3.750,54.671180,389.700,3.182,10.370,0.0,3.558,2.654,55.221384,627.100,2.618,10.890,0.0,2.560,2.952,55.783815,945.300,2.109,11.340,0.0,-1.172,6.135,56.264774,543.400,0.014,17.030,0.0,3.525,-0.978,56.363399,1331.800,1.654,11.890,0.0,-2.378,6.547,56.968211,1746.600,1.255,12.230,0.0,-3.545,6.451,57.612486,2120.100,0.910,12.620,0.0,-5.416,6.056,58.323877,2363.700,0.621,12.950,0.0,-1.932,0.436,58.446588,1442.100,0.083,14.910,0.0,6.768,-1.273,59.164204,2379.900,0.387,13.530,0.0,-6.561,2.309,59.590983,2090.700,0.207,14.080,0.0,6.957,-0.776,60.306056,2103.400,0.207,14.150,0.0,-6.395,0.699,60.434778,2438.000,0.386,13.390,0.0,6.342,-2.825,61.150562,2479.500,0.621,12.920,0.0,1.014,-0.584,61.800158,2275.900,0.910,12.630,0.0,5.014,-6.619,62.411220,1915.400,1.255,12.170,0.0,3.029,-6.759,62.486253,1503.000,0.083,15.130,0.0,-4.499,0.844,62.997984,1490.200,1.654,11.740,0.0,1.856,-6.675,63.568526,1078.000,2.108,11.340,0.0,0.658,-6.139,64.127775,728.700,2.617,10.880,0.0,-3.036,-2.895,64.678910,461.300,3.181,10.380,0.0,-3.968,-2.590,65.224078,274.000,3.800,9.960,0.0,-3.528,-3.680,65.764779,153.000,4.473,9.550,0.0,-2.548,-5.002,66.302096,80.400,5.200,9.060,0.0,-1.660,-6.091,66.836834,39.800,5.982,8.580,0.0,-1.680,-6.393,67.369601,18.560,6.818,8.110,0.0,-1.956,-6.475,67.900868,8.172,7.708,7.640,0.0,-2.216,-6.545,68.431006,3.397,8.652,7.170,0.0,-2.492,-6.600,68.960312,1.334,9.650,6.690,0.0,-2.773,-6.650,118.750334,940.300,0.010,16.640,0.0,-0.439,0.079,368.498246,67.400,0.048,16.400,0.0,0.000,0.000,424.763020,637.700,0.044,16.400,0.0,0.000,0.000,487.249273,237.400,0.049,16.000,0.0,0.000,0.000,715.392902,98.100,0.145,16.000,0.0,0.000,0.000,773.839490,572.300,0.141,16.200,0.0,0.000,0.000,834.145546,183.100,0.145,14.700,0.0,0.000,0.000];
table1_data = reshape(table1_data, 7, []).';
Oxygen_Attenuation = table( ...
    table1_data(:,1), ...
    table1_data(:,2), ...
    table1_data(:,3), ...
    table1_data(:,4), ...
    table1_data(:,5), ...
    table1_data(:,6), ...
    table1_data(:,7), ...
    'VariableNames', {'f0','a1','a2','a3','a4','a5','a6'} ...
);

% ITU-R P.676-13 -Table 2 - Spectroscopic Data for Water Vapour Attenuation
table2_data=[22.235080,.1079,2.144,26.38,.76,5.087,1.00,67.803960,.0011,8.732,28.58,.69,4.930,.82,119.995940,.0007,8.353,29.48,.70,4.780,.79,183.310087,2.273,.668,29.06,.77,5.022,.85,321.225630,.0470,6.179,24.04,.67,4.398,.54,325.152888,1.514,1.541,28.23,.64,4.893,.74,336.227764,.0010,9.825,26.93,.69,4.740,.61,380.197353,11.67,1.048,28.11,.54,5.063,.89,390.134508,.0045,7.347,21.52,.63,4.810,.55,437.346667,.0632,5.048,18.45,.60,4.230,.48,439.150807,.9098,3.595,20.07,.63,4.483,.52,443.018343,.1920,5.048,15.55,.60,5.083,.50,448.001085,10.41,1.405,25.64,.66,5.028,.67,470.888999,.3254,3.597,21.34,.66,4.506,.65,474.689092,1.260,2.379,23.20,.65,4.804,.64,488.490108,.2529,2.852,25.86,.69,5.201,.72,503.568532,.0372,6.731,16.12,.61,3.980,.43,504.482692,.0124,6.731,16.12,.61,4.010,.45,547.676440,.9785,.158,26.00,.70,4.500,1.00,552.020960,.1840,.158,26.00,.70,4.500,1.00,556.935985,497.0,.159,30.86,.69,4.552,1.00,620.700807,5.015,2.391,24.38,.71,4.856,.68,645.766085,.0067,8.633,18.00,.60,4.000,.50,658.005280,.2732,7.816,32.10,.69,4.140,1.00,752.033113,243.4,.396,30.86,.68,4.352,.84,841.051732,.0134,8.177,15.90,.33,5.760,.45,859.965698,.1325,8.055,30.60,.68,4.090,.84,899.303175,.0547,7.914,29.85,.68,4.530,.90,902.611085,.0386,8.429,28.65,.70,5.100,.95,906.205957,.1836,5.110,24.08,.70,4.700,.53,916.171582 8.400 1.441 26.73 .70 5.150 .78,923.112692 .0079 10.293 29.00 .70 5.000 .80,970.315022 9.009 1.919 25.50 .64 4.940 .67,987.926764 134.6 .257 29.85 .68 4.550 .90,1780.000000, 17506. .952, 196.3, 2.00, 24.15, 5.00];
table2_data = reshape(table2_data, 7, []).';
WaterVapour_Attenuation = table( ...
    table2_data(:,1), ...
    table2_data(:,2), ...
    table2_data(:,3), ...
    table2_data(:,4), ...
    table2_data(:,5), ...
    table2_data(:,6), ...
    table2_data(:,7), ...
    'VariableNames', {'f0','b1','b2','b3','b4','b5','b6'} ...
);

h_lower = 0;
h_upper = 100;
line_centres = Oxygen_Attenuation.f0;
f = linspace(1e-6, 1000, 1001); % 0 - 1000 GHz
f = unique([f, line_centres']); %Uncomment to include line centres

% For the Specific Attenuation only
p = 1013.25; % hPa
T = 288.15; % Kelvin
WV_Density = [0, 7.5]; % Water Vapour Density g/m^3
WV_PartialPres = WV_Density*T/216.7; % Water Vapour Partial Pressure

SpecificAttenuation_dB = SpecificAttenuation(Oxygen_Attenuation, WaterVapour_Attenuation, f, p, T, WV_PartialPres, 1);

ZenithAttenuation_dB = PositiveSlantPathAttenuation(Oxygen_Attenuation, WaterVapour_Attenuation, h_lower, h_upper, f);

figure;
semilogy(f, ZenithAttenuation_dB(1,:), 'b-', 'LineWidth', 1.8); hold on;
semilogy(f, ZenithAttenuation_dB(2,:), 'r--', 'LineWidth', 1.8);
ylim([1e-3, 1e5]);
xlim([0, 1000]);
title('ITU-R P.676-13, Figure 4: Zenith Attenuation due to Atmospheric Gasses, Including Line Centres', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Frequency (GHz)', 'FontSize', 12);
ylabel('Zenith Attenuation (dB)', 'FontSize', 12);
legend({'Dry', 'Standard'}, 'FontSize', 11, 'Location', 'northwest');
grid on;
set(gca, 'FontSize', 11);

figure;
semilogy(f, SpecificAttenuation_dB(1,:), 'b', 'LineWidth', 1.8); hold on;
semilogy(f, SpecificAttenuation_dB(2,:), 'r--', 'LineWidth', 1.8);
ylim([1e-3, 1e5]);
xlim([0, 1000]);
title('ITU-R P.676-13, Figure 1: Specific Attenuation due to Atmospheric Gasses', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Frequency (GHz)', 'FontSize', 12);
ylabel('Specific Attenuation (dB/km)', 'FontSize', 12);
legend({'Dry', 'Standard'}, 'FontSize', 11, 'Location', 'northwest');
grid on;
set(gca, 'FontSize', 11);

function [Attenuation_dB] = SpecificAttenuation(Oxygen_Attenuation, WaterVapour_Attenuation, f, p, T, WV_PartialPres, PathLength)
    theta = 300/T;
    Attenuation_dB = zeros(length(WV_PartialPres), length(f));

    for WV_PartialPres_i = 1:length(WV_PartialPres)
        d = 5.6 * (1e-4) * (p+WV_PartialPres(WV_PartialPres_i)) * theta^0.8; % Width Parameter Acc. to Debye Spectrum 
    
        % Correction Factors
        Dirag_Oxy = (Oxygen_Attenuation.a5 + Oxygen_Attenuation.a6.*theta) * 1e-4 .* (p+WV_PartialPres(WV_PartialPres_i)) .* theta^0.8;
        Dirag_WVapour = 0;
        
        % Line Widths 
        df_Oxy = Oxygen_Attenuation.a3 .* 1e-4 .* (p .* theta.^(0.8 - Oxygen_Attenuation.a4) + 1.1*WV_PartialPres(WV_PartialPres_i)*theta);
        df_WVapour = WaterVapour_Attenuation.b3 .* 10^-4 .* ...
            (p.*(theta.^WaterVapour_Attenuation.b4) + WaterVapour_Attenuation.b5.*WV_PartialPres(WV_PartialPres_i).*(theta.^WaterVapour_Attenuation.b6));
        
        % Line Strengths
        S_Oxygen = Oxygen_Attenuation.a1 * 10^-7 * p .* (theta^3) .* exp(Oxygen_Attenuation.a2*(1-theta)); 
        S_WaterVapour = WaterVapour_Attenuation.b1 * 10^-1 * WV_PartialPres(WV_PartialPres_i) * theta^3.5 .* exp(WaterVapour_Attenuation.b2*(1-theta));
        
        % Line Shape Factors
        F_Oxygen = zeros(length(Oxygen_Attenuation.f0), length(f));
        F_WaterVapour = zeros(length(WaterVapour_Attenuation.f0), length(f));
        for fi = 1:length(f)
            F_Oxygen(:,fi) = (f(fi) ./ Oxygen_Attenuation.f0) ...
            .* ( ((df_Oxy - Dirag_Oxy .* (Oxygen_Attenuation.f0 - f(fi))) ./ ((Oxygen_Attenuation.f0 - f(fi)).^2 + df_Oxy.^2)) ...
               + ((df_Oxy - Dirag_Oxy .* (Oxygen_Attenuation.f0 + f(fi))) ./ ((Oxygen_Attenuation.f0 + f(fi)).^2 + df_Oxy.^2)) ); 
            F_WaterVapour(:,fi) = (f(fi) ./ WaterVapour_Attenuation.f0) ...
            .* ( ((df_WVapour - Dirag_WVapour .* (WaterVapour_Attenuation.f0 - f(fi))) ./ ((WaterVapour_Attenuation.f0 - f(fi)).^2 + df_WVapour.^2)) ...
               + ((df_WVapour - Dirag_WVapour .* (WaterVapour_Attenuation.f0 + f(fi))) ./ ((WaterVapour_Attenuation.f0 + f(fi)).^2 + df_WVapour.^2)) );
        end

        % Dry Continuum for Oxygen
        N_D = f .* p * theta^2 .* ((6.14*10^-5)./(d*(1 + (f/d).^2)) + (1.4*(10^-12)*p*theta^1.5)./(1+ 1.9*(10^-5).*f.^1.5));
   
        N_Oxygen = sum(S_Oxygen.*F_Oxygen, 1) + N_D;  

        N_WaterVapour = sum(S_WaterVapour.*F_WaterVapour,1);

        Attenuation_dB(WV_PartialPres_i,:) = 0.1820 .* f .* (N_Oxygen + N_WaterVapour) * PathLength;
    
    end
end

function [T, p, WV_PartialPres] = GetAtmosphereConditions(Z)
    H = 6356.766 * Z / (6356.766 + Z); % Geopotential Height
    a0 = 95.571899;
    a1 = -4.011801;
    a2 = 6.424731e-2;
    a3 = -4.789660e-4;
    a4 = 1.340543e-6;

    if H < 84.852
        if 0 <= H && H <= 11
            T = 288.15 - 6.5*H;
            p = 1013.25 * (288.15/(288.15 - 6.5*H))^(-34.1632/6.5);
        elseif  11 < H && H <= 20
            T = 216.65;
            p = 226.3226 * exp (-34.1632*(H-11)/216.65);
        elseif 20 < H && H <= 32
            T = 216.65 + (H - 20);
            p = 54.74980 * (216.65/(216.65 + (H-20)))^34.1632;
        elseif 32 < H && H <= 47
            T = 228.65 + 2.8*(H-32);
            p = 8.680422 * (228.65/(228.65 + 2.8*(H-32)))^(34.1632/2.8);
        elseif 47 < H && H <=51
            T = 270.65;
            p = 1.109106 * exp(-34.1632*(H-47)/270.65);
        elseif 51 < H && H <= 71
            T = 270.65 - 2.8*(H-51);
            p = 0.6694167 * (270.65 / (270.65 - 2.8*(H-51)))^(-34.1632/2.8);
        elseif 71 < H && H <= 84.852
            T = 214.65 - 2*(H-71);
            p = 0.03956649 * (214.65/(214.65 - 2*(H-71)))^(-34.1632/2);
        end
   
    elseif 86 <= Z 
        if 86 <= Z && Z<=91
            T = 186.8673;     
        elseif 91 < Z && Z <= 100
            T = 263.1905 - 76.3232*sqrt(1 - ((Z-91)/19.9429)^2);
        end
        p = exp(a0 + a1*Z + a2*Z^2 + a3*Z^3 + a4*Z^4);
    end

    WV_Density = 7.5*exp(-Z/2);
    WV_PartialPres(1) = 0;
    WV_PartialPres(2) = WV_Density*T/216.7;

end

function [Attenuation_dB] = PositiveSlantPathAttenuation(Oxygen_Attenuation, WaterVapour_Attenuation, h_lower, h_upper, f)
    i_lower = floor(100 * log(1e4*h_lower*(exp(1/100) - 1) + 1) + 1);
    i_upper = ceil(100 * log(1e4*h_upper*(exp(1/100) - 1) + 1) + 1);
    m = (exp(2/100) - exp(1/100))/(exp(i_upper/100) - exp(i_lower/100)) * (h_upper - h_lower);
    
    r1 = 6371;  
    beta1 = 0;
    Dirag_CurrentLayer = 0;
    Dirag1 = m;
    [T1, p1, WV_PartialPres] = GetAtmosphereConditions(Dirag1/2);
    n1 = 77.6*p1/T1 - 5.6*WV_PartialPres/T1 + 3.75e5*WV_PartialPres/(T1^2);
    
    if i_upper > 922
        i_upper = 922;
    end
  
    Attenuation_dB = zeros(length(WV_PartialPres), length(f));
   
    for i = i_lower:i_upper
        Dirag_PrevLayer = Dirag_CurrentLayer;
        Dirag_CurrentLayer = m * (exp((i-1)/100)); 
        h = h_lower + (m * (exp((i-1)/100) - exp((i_lower-1)/100)) / (exp(1/100)-1));

        if i == i_lower
            r = r1 + h;
        else
            r = r + Dirag_PrevLayer;
        end
        
        [T, p, WV_PartialPres] = GetAtmosphereConditions(h + Dirag_CurrentLayer/2);

        n_i = 77.6*p/T - 5.6*WV_PartialPres/T + 3.75e5*WV_PartialPres/(T^2);
        beta_i = asin((n1*r1/(n_i*r))*sin(beta1));

        PathLength = -r*cos(beta_i) + sqrt((r^2)*cos(beta_i)^2 + 2*r*Dirag_CurrentLayer + Dirag_CurrentLayer^2);
        
        % Path Lenght will have the size (2x1) if beta1 is not zero, change accordingly if needed
        Attenuation_dB(1, :) = Attenuation_dB(1, :) + SpecificAttenuation(Oxygen_Attenuation, WaterVapour_Attenuation, f, p, T, WV_PartialPres(1), PathLength);
        Attenuation_dB(2, :) = Attenuation_dB(2, :) + SpecificAttenuation(Oxygen_Attenuation, WaterVapour_Attenuation, f, p, T, WV_PartialPres(2), PathLength);
    end
end



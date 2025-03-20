% 6. Compute and plot the average cell downlink capacity of a OFDMA system and the bit-rate provided to an edge user versus the reuse factor N in case of:
% -  Full Frequency Reuse with reuse factor N=1,3,7
% -  SIR evaluation through formula (no simulator)
% -  Capacity evaluation through Shannon formula with attenuation factor alpha =0.75
% -  Large scale fading dB = 4 dB;
% -  Coverage probability F = 0.98;
% -  Path Loss exponent gamma = 4
% -  Path loss at 1km L0_db = 125.13 dB
% -  Noise figure NF =3.5 dB
% -  Noise spectral density kT0 = -174 dBm Hz-1
% -  Total System Bandwidth Bw = 5 MHz
% -  Bandwidth for single users 180 kHz
% -  Sectoring = 3
% -  No power control
% -  Cell Radius R=2 km
% -  Base Station Transmit Power Ptx_dbm= 40 dBm
% -  Cell load ro = 0.7					
% ------------------------------------------------ ----------------------------------------------
		
% Downlink capacity - ?
% Rb edge user vs N - ?

clear all;
RF = [1 3 7];
N = 1;
alpha = 0.75;
F = 0.98;
gamma = 4; 
L0_db = 125.13; 
R0 = 1; % km 
NF_db = 3.5; 
kt0_dbm = -174; 
Bw_tot = 5e6; 
Bw_eu = 180e3 
sec = 3;
R = 2; % km
Ptx_dbm = 40; 
ro = 0.7;
sigma_db = 4;  

Bw_eu_db = 10 * log10(Bw_eu);
N_dbm = kt0_dbm + Bw_eu_db + NF_db;                                         % Noise Power as a combination of Noise Spectral Density kt0, Bandwidth Bw_eu + Noise Power F
M_db = fzero(@(x) 0.5*erfc(-x/(sigma_db*sqrt(2)))-F, 0.5)
L_db = L0_db + 10 * gamma * log10(R/R0);                                    % Log-Distance Path Loss Model, R0 - reference distance for which the L0_dB is known
S_dbm = Ptx_dbm - L_db - M_db;                                         % Minimum required Signal Power at the cell edge from Link Budget L_db = Pt_db + G_db - S_dbm - M; 

% SNR
SNR_db = S_dbm - N_dbm;                                                    
SNR = 10^(SNR_db/10);

for i = 1:length(RF)
    N = RF(i);

    % SIR
    SIR = (1/ro) * sec * (1/6) * sqrt(3 * N)^gamma; 

    % SINR
    SINR = inv(inv(SNR)+inv(SIR));                                              % Compute SINR based on SNR and SIR

    Reu(i) = alpha * Bw_eu * log2(1 + SINR);                                       % The bit rate for a user at the cell edge
    f = @(D,R,SINR,gamma) alpha * log2(1 + (((D./R).^-gamma).*SINR)).*(2 * D./(R^2));     % Capacity function, defines the formula inside integral, SINR_db converted to linear
    eta_c = integral(@(D) f(D,R*1000,SINR,gamma), 1, R*1000);                          % R is in meters (by default consider R = 1000 m

    % Network Capacity of a cell
    cap(i) = eta_c * Bw_tot / N;                                   

    disp(['For Reuse Factor N = ', num2str(N)]);
    disp(['  Downlink Capacity: ', num2str(cap(i) / 1e6), ' Mbit/s']);
    disp(['  Edge User Bit Rate: ', num2str(Reu(i) / 1e6), ' Mbit/s']);
end

% Plot 1: Cell Downlink Capacity vs Reuse Factor
plot(RF,cap/1e6); 
xlabel('Reuse Factor (N)');
ylabel('Downlink Capacity (Mbit/s)');
title('Average Cell Downlink Capacity vs Reuse Factor');
grid on;
xticks(RF);     % x-axis values
legend('Capacity');

% Plot 2: Edge User Bit Rate vs Reuse Factor
figure;
plot(RF, Reu/1e6);
xlabel('Reuse Factor (N)');
ylabel('Edge User Bit Rate (Mbit/s)');
title('Edge User Bit Rate vs Reuse Factor');
grid on;
xticks(RF);
legend('Edge User Bit Rate');
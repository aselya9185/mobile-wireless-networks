% compute and plot the bit-rate provided to an edge user versus the cell load ro [0.5 0.6 0.7 0.8 0.9]

% Uplink power control
clear all;
N = 1;                     % Frequency reuse factor
sec = 3;              % Number of sectors
alpha = 0.75;               % Attenuation factor
sigma_db = 4;               % Standard deviation for fading
F = 0.95;                   % Coverage probability
gamma = 4;                  % Path loss exponent
R0 = 1000;                  % Reference distance (m)
L0_db = 141.44;             % Path loss at reference distance (dB)
NF_db = 3.5;                % Noise figure (dB)
kt0_dbm = -174;             % Noise spectral density (dBm/Hz)
Bw_eu = 360e3;              % Bandwidth per user (Hz)
R = 600;                    % Cell radius (m)
Ptx_x_sc_dbm = 18;          % Base station transmit power per subcarrier (dBm)
ro_values = [0.5 0.6 0.7 0.8 0.9];  % Cell load
% ro = 0.5;
Rue = R;                    % Edge user distance
pcu = 1;                    % Power control
subcarriers = 24;           % Number of subcarriers

Reu_values = zeros(size(ro_values));  % Array to store edge user bit rates

% SNR
Bw_eu_db = 10 * log10(Bw_eu);
N_dbm = kt0_dbm + Bw_eu_db + NF_db;                                         % Noise Power as a combination of Noise Spectral Density kt0, Bandwidth Bw_eu + Noise Power F
M_db = fzero(@(x) 0.5*erfc(-x/(sigma_db*sqrt(2)))-F,0.5);                   % Shadowing margin (fading margin) to have a coverage probability F, fzero finds the root of the function f(x)=0, for single-variable equations.
L_db = L0_db + 10 * gamma * log10(R/R0);                                    % Log-Distance Path Loss Model, R0 - reference distance for which the L0_dB is known
Pt_db = Ptx_x_sc_dbm + 10*log10(subcarriers);                               % Total BS transmit power
S_dbm = Pt_db - L_db - M_db;                                                % Minimum required Signal Power at the cell edge from Link Budget L_db = Pt_db + G_db - S_dbm - M; 
SNR_db = S_dbm - N_dbm;                                                    
SNR = 10^(SNR_db/10);
fprintf('SNR: %.2f \n', SNR);


for i = 1:length(ro_values)
    ro = ro_values(i);
    % SIR
    [SIRup SIRdown]=InterferenceComputation(R,gamma,ro,N,sec,Rue,pcu);          % SIR for both uplink and downlink, when pcu is given, using the function 'InterferenceComputation'
    SIR = SIRup;
    fprintf('ro = %.1f\n', ro);
    fprintf('SIR: %.2f \n', SIR);
    
    % SINR
    SINR = inv(inv(SNR)+inv(SIR));                                              % Compute SINR based on SNR and SIR
    fprintf('SINR: %.2f \n', SINR);
    
    % Alpha-Shannon capacity (with attenuation): C = alpha * B * log2(1 + SINR)
    Reu_values(i) = alpha * Bw_eu * log2(1 + SINR);                             % Bit rate for edge user
    fprintf('Reu bit-rate provided to an edge user: %.2f bps\n\n', Reu_values(i));

end

% Plot bit-rate vs. cell load (ro)
figure;
plot(ro_values, Reu_values, '-o');
xlabel('Cell Load (ro)');
ylabel('Bit-rate at Edge User (bps)');
title('Bit-rate vs. Cell Load (ro)');
grid on;
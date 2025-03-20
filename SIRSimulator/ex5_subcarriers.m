% Cell load ro - ? for R=700m

R = 700;    %m
N = 1;
sec = 3;
F = 0.95;
alpha = 0.75;
sigma_db = 4;
gamma = 4;
L0_db = 141.44;
R0 = 1000;  %m
NF_db = 3.5;
kt0_dbm = -174;
Br_eu = 3.5e6;    % bit rate for edge user [bit/s]
subcarriers = 144;
scs = 15e3;     % Subcarrier spacing Hz
Bw_eu = scs * subcarriers;  % Bandwidth for edge user
% no power control
Ptx_x_sc_dbm = 18;  % BS tx power per subcarrier

Bw_eu_db = 10 * log10(Bw_eu);
N_dbm = kt0_dbm + Bw_eu_db + NF_db;                                         % Noise Power as a combination of Noise Spectral Density kt0, Bandwidth Bw_eu + Noise Power F
M_db = fzero(@(x) 0.5*erfc(-x/(sigma_db*sqrt(2)))-F,0.5);                   % Shadowing margin (fading margin) to have a coverage probability F, fzero finds the root of the function f(x)=0, for single-variable equations.
L_db = L0_db + 10 * gamma * log10(R/R0);                                    % Log-Distance Path Loss Model, R0 - reference distance for which the L0_dB is known
Pt_db = Ptx_x_sc_dbm + 10*log10(subcarriers);                               % Total BS transmit power
S_dbm = Pt_db - L_db - M_db;                                                % Minimum required Signal Power at the cell edge from Link Budget L_db = Pt_db + G_db - S_dbm - M;

SNR_db = S_dbm - N_dbm;  
SNR = 10^(SNR_db/10);
fprintf('SNR: %.2f \n', SNR);

% Alpha-Shannon capacity (with attenuation): C = alpha * B * log2(1 + SINR)
SINR = 2^(Br_eu/(alpha * Bw_eu)) - 1;                                         % SINR through Shannon Capacity Reu = alpha * Bw * log2(1 + SINR)
fprintf('SINR: %.2f \n', SINR);

SIR = inv(inv(SINR) - inv(SNR));                                              % Compute SIR based on SNR and SIR
fprintf('SIR: %.2f \n', SIR);

ro = sec * sqrt(3 * N)^gamma / 6 / SIR;
fprintf('cell load: %.2f \n', ro);

% Output:
% SNR: 3.10 
% SINR: 3.47 
% SIR: -29.16 
% cell load: -0.15 
% 5. Compute the Tx power of the base station (dBm) to provide user at cell edge with a bit rate of 250 kbit/sec, assuming that:					
% - Soft Frequency Reuse (SFR) with 3 carrier frequencies (N=3) and external/central Tx power ratio omega = 4 (omega: Ptx_e / Ptx_c)
% - Cell load ro = 0.5
% - SIR evaluation through formula (no simulator)
% - Cell radius = 2 km
% - Capacity evaluation through Shannon formula with attenuation factor alpha = 0.75 
% - User Bandwidth Bw_eu = 180 kHz, assigned for the whole time (no TDMA)
% - Large scale fading sigma_db = 4 dB;
% - Coverage probability F = 0.98; 
% - Path Loss exponent gamma = 4 
% - Path loss at 1km = 125.13 dB
% - Noise figure NF =3.5 dB
% - Noise spectral density kT0 = -174 dBm Hz-1
% 
% ------------------------------------------------ ----------------------------------------------
% 
% Solution Pt_dbm= 33.242247


clear all;
Br_eu = 0.25e6; 
N = 3;          % in SFR only in the cell edge Frequency is reused, in the center Bw_tot
omega = 4;      % Ptx_e/Ptx_c 
ro = 0.5; 
R = 2; %km
alpha = 0.75; 
Bw_eu = 180e3 
sigma_db = 4; 
F = 0.98; 
gamma = 4; 
L0_db = 125.13; 
R0 = 1; %km 
NF_db = 3.5; 
kt0_dbm = -174; 


% Ptx - ?

Bw_eu_db = 10 * log10(Bw_eu);
N_dbm = kt0_dbm + Bw_eu_db + NF_db; 

% SINR
SINR = 2^(eta_c/alpha) - 1;                                                 % SINR through Shannon Capacity Reu = alpha * Bw * log2(1 + SINR) = eta_c * Bw
SINRe = 2^(Br_eu/(alpha * Bw_eu)) - 1;                                      % SINR through Shannon Capacity Reu = alpha * Bw * log2(1 + SINR)
SINRe_db = 10 * log10(SINRe);

% SIR
SIRe = (1/ro) * (1 / (((6/omega) * sqrt(3)^(- gamma)) + (6 * sqrt(3 * N)^(-gamma))));         % SIR for SFR systems at the cell edge (d=R -> d/R = 1), and considering a cell load
Im_db = 10 * log10(inv(1 - SINRe / SIRe));                                  % Interference margin in dB

S_dbm = SINRe_db + N_dbm + Im_db;                                          % Minimum required Signal Power at the cell edge Smin FOR SFR
M_db = fzero(@(x) 0.5*erfc(-x/(sigma_db*sqrt(2)))-F, 0.5);
Pr_dbm = S_dbm + M_db; 

L_db = L0_db + 10 * gamma * log10(R/R0);                                    % Log-Distance Path Loss Model, R0 - reference distance for which the L0_dB is known

Pt_dBm = Pr_dbm + L_db;

fprintf('Pt_dbm: %f', Pt_dBm);      % Pt_dbm: 33.242247
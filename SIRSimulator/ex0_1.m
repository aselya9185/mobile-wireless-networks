% Find:
% power received (Prx) at the cell edge in order to achieve a cell load of ro=0.7 

clear all;
N = 1;
R = 600;    % km
ro = 0.7;
sec = 3;
alpha = 0.75;
delta = 10e6;   % bit/s/km2
sigma_db = 6;
F = 0.95;
NF_db = 3.5;
kt0_dbm = -174;
gamma = 4;
Bw_tot = 5e6;
Bw_eu = 180e3;

% Prx - ?
% 
% MY SOLUTION
SIR = (1/ro) * sec * (1/6) * sqrt(3 * N)^gamma;                             % SIR for hexagon shaped cell Ni = 6, K = sectors, N = reuse factor

eta_c = delta * pi * (R)^2 / ro / Bw_tot; 
SINR = 2^(eta_c/alpha) - 1;                                                 % SINR through Shannon Capacity Reu = alpha * Bw * log2(1 + SINR) = eta_c * Bw

SNR = inv(inv(SINR)-inv(SIR));
SNR_db = 10 * log10(SNR);

Bw_eu_db = 10 * log10(Bw_eu);
N_dbm = kt0_dbm + Bw_eu_db + NF_db;                                         % Noise Power as a combination of Noise Spectral Density kt0, Bandwidth Bw_eu + Noise Power F
S_dbm = SNR_db + N_dbm;                                                     % Minimum required Signal Power at the cell edge

M_db = fzero(@(x) 0.5*erfc(-x/(sigma_db*sqrt(2)))-F, 0.5);
Pr_dbm = S_dbm + M_db;                                                      % Received signal power at the cell edge 
fprintf('Received power: %f', Pr_dbm);

% check solution ex 0

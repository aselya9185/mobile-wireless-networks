% Find:
% minimum number of cells (Nc) 
% maximum transmission power (Ptx_dbm_cell) of the base station

clear all;

A_city = 1285;      % City area in km²
delta_u = 2232;     % Users per km²
C_u = 1e5;          % User bit rate (bps)
Bw_tot = 100e6;     % Total system bandwidth (Hz)
Bw_eu = 360e3;      % Bandwidth for edge user (Hz)
ro = 0.7;           % Cell load
sec = 3;            % Sectoring
G_db = 8;           % Base Station Antenna Gain (dB)
sigma_db = 6;       % Large Scale Fading (dB)
F = 0.95;           % Coverage Probability
NF_db = 3.5;        % Noise Figure (dB)
kt0_dbm = -174;     % Noise Spectral Density (dBm/Hz)
L0_db = 141.43;     % Path Loss at 1 km (dB)
R0 = 1;             % km
gamma = 4;          % Path Loss exponent
N = 1;              % Reuse Factor
alpha = 0.75;       % Shannon attenuation factor
SINR_min_db = 8;    % SINR min
% Note: The maximum TX power of the base station (linear) Ptx_cell is the
% TX power towards users Ptx_eu * Bw_tot/Bw_eu

% Ncell - ?
% Ptx_dbm_cell - ?      maximum transmission power of the base station

% City capacity
Ccity = C_u * delta_u * A_city;

% Cell capacity
% Ccell = C_u * delta * pi * Rcell^2
% eta_c = Ccell / ro / Bw_cell ---> 
ec = eta_c(R0*1000, SINR_min_db, gamma, alpha);
Ccell = ec * ro * Bw_tot / N;

% Number of cells
Ncell = ceil(Ccity / Ccell);                                                % Rounding up to cover city capacity
fprintf('Number of cells: %.2f \n', Ncell);

Bw_eu_db = 10 * log10(Bw_eu);
N_dbm = kt0_dbm + Bw_eu_db + NF_db;                                         % Noise Power as a combination of Noise Spectral Density kt0, Bandwidth Bw_eu + Noise Power F

SIR = (1/ro) * sec * (1/6) * sqrt(3 * N)^gamma; 

SINR = 2^(C_u/(alpha * Bw_eu)) - 1;                                         % SINR through Shannon Capacity Reu = alpha * Bw * log2(1 + SINR)


SNR = inv(inv(SINR)-inv(SIR));
SNR_db = 10 * log10(SNR);

M_db = fzero(@(x) 0.5*erfc(-x/(sigma_db*sqrt(2)))-F, 0.5);

S_dbm = SNR_db + N_dbm;                                                     % Minimum required Signal Power at the cell edge
Pr_dbm = S_dbm + M_db;                                                      % Received signal power at the cell edge 

R = sqrt(A_city/(pi * Ncell));	 %km                                        % Cell Radius from A_city / Ncell = pi * Rcell^2
L_db = L0_db + 10 * gamma * log10(R/R0);                                    % Log-Distance Path Loss Model, R0 - reference distance for which the L0_dB is known

Pt_dbm = Pr_dbm - G_db + L_db;

% 
% Ptx_dbm_cell 

fprintf('Transmit power: %.2f \n', Pt_dbm);

% Output:
% Number of cells: 973.00 
% Transmit power: 15.70 
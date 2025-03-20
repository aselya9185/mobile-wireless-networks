% Find:
% 1) cell radius R 
% 2) maximum number of users per cell Nu_per_cell
% 3) number of cells Ncell (overlap factor of 1.2) 
% 4) maximum number of users Nu_per_city in the city

clear all;

% Given parameters
gamma = 4;                    % Path loss exponent
ro = 0.7;                     % Cell load
sec = 3;                      % Number of sectors
G_db = 8;                     % Antenna gain (dB)
SINR_min_sub_db = 0;          % Minimum SINR at cell edge (dB)
sigma_db = 6;                 % Shadow fading standard deviation (dB)
F = 0.98;                     % Probability of SINR threshold
NF_db = 3.5;                  % Noise figure (dB)
kt0_dbm = -174;               % Noise spectral density (dBm/Hz)
R0 = 1000;                    % Reference distance (m)
L0_db = 141.43;               % Path loss at reference distance (dB)
Bw_tot = 100e6;               % Total bandwidth (Hz)
N = 1;                        % Frequency reuse factor
Acity = 1285;                 % City area (km^2)
C_u = 1e5;                    % User's average bit rate (bps)
alpha = 0.75;                 % Attenuation factor for spectral efficiency
Bw_sub = 15e3;                % Subcarrier bandwidth (Hz)
Ptx_sub_dbm = 10*log10((8000/5e6) * Bw_sub);  % Transmit power per subcarrier (dBm)
overlap_area_factor = 1.2;    % Overlap factor for cells


% Find Prx = Smin + M for Link Budget Lp
% Smin:
SIR = (1/ro) * sec * (1/6) * sqrt(3 * N)^gamma;                             % SIR for hexagon shaped cell Ni = 6, K = sectors, N = reuse factor
SINR = 10^(SINR_min_sub_db / 10);                                           % SINR linear
Im_db = 10 * log10(inv(1 - SINR / SIR));                                    % Interference Margin
Bw_sub_db = 10 * log10(Bw_sub);                                             % Bw_sub in dB
N_dbm = kt0_dbm + Bw_sub_db + NF_db;                                        % Noise Power as a combination of Noise Spectral Density kt0, Bandwidth Bw_eu + Noise Power F
S_dbm = SINR_min_sub_db + N_dbm + Im_db;                                    % Minimum required Signal Power at the cell edge Smin

M_db = fzero(@(x) F-0.5*erfc(-x/(sigma_db*sqrt(2))),10);                    % FADING MARGIN M

Pr_dbm = S_dbm + M_db;                                                      % FADING MARGIN M

subcarriers = Bw_tot / Bw_sub;                                              % Number of subcarriers
Pt_dbm = Ptx_sub_dbm + 10*log10(subcarriers);                               % Total BS transmit power
L_db = Pt_dbm + G_db - Pr_dbm;                                              % Allowed path loss from a link budget

% Cell Radius
R = R0 * 10^((L_db - L0_db) / (10 * gamma));                                % in dB, maximum cell radius from Log-Distance Path Loss Model                           


% Find Max Users per Cell
f = @(D,R,SINR_db,gamma) alpha * log2(1 + (((D./R).^-gamma).*10^(SINR_db/10))).*(2 * D./(R^2));     % Capacity function, defines the formula inside integral, SINR_db converted to linear
eta_c = integral(@(D) f(D,R,SINR_min_sub_db,gamma), 1, R);                  % Average Spectral Efficiency, R is in meters
% fprintf('eta_c: %.2f \n', eta_c);

Ccell = eta_c * Bw_tot / N;                                                 % Cell Capacity
% Max Users per Cell
Nu_per_cell = Ccell / C_u * ro;

% Number of cells
Ncell = overlap_area_factor * Acity / (pi * (R/1000)^2);  

% Max Number of Users in the City
Nu_per_city = Nu_per_cell * Ncell;

% Display results
fprintf('Cell Radius R: %.2f m\n', R);
fprintf('Max Users per Cell: %.2f users\n', Nu_per_cell);
fprintf('Number of Cells in the City: %.2f cells\n', Ncell);
fprintf('Max Number of Users in the City: %.2f users\n', Nu_per_city);

% Output:
% Cell Radius R: 7200.08 m
% Max Users per Cell: 1714.74 users
% Number of Cells in the City: 9.47 cells
% Max Number of Users in the City: 16235.23 users
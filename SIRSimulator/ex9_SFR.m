sigma_db = 6;
F = 0.95;
ro = 0.7;
N = 3;  % SFR
omega = 2;
sec = 3;
gamma = 4;
NF_db = 3.5;
kt0_dbm = -174;
Bw_eu = 512e3;
P0_dbm = 10 * log10(0.1 * 1000); % Transmit power at 10m in dBm (converted from 0.1 W). This is a reference point for computing how much power is received over distance using the path loss formula.
R0 = 10;
SINR_min_db = 0;

% Rmax - ?

Bw_eu_db = 10 * log10(Bw_eu);
N_dbm = kt0_dbm + Bw_eu_db + NF_db; 

% one way of finding SINR, output maximum cell radius: 5620.20 
f = @(D,R,SINR_db,gamma) alpha * log2(1 + (((D./R).^-gamma).*10^(SINR_db/10))).*(2 * pi * D./(pi * R^2));     % Capacity function, defines the formula inside integral, SINR_db converted to linear
eta_c_req = integral(@(D) f(D,1000,SINR_min_db,gamma), 1, 1000);                          % R is in meters (by default consider R = 1000 m
eta_c = @(SINR_db,gamma) integral(@(D) f(D,1000,SINR_db,gamma), 1, 1000);         % anonymous function to call eta_c(SINR_db,gamma)
SINRe_db = fzero(@(SINR_db) eta_c(SINR_db,gamma) - eta_c_req, 3);            % MWN_ex.7


% % SINR another way, output maximum cell radius: 4153.73 
% f = @(D,R,SINR_db,gamma) 0.75 * log2(1 + (((D./R).^-gamma).*10^(SINR_db/10))).*(2 * pi * D./(pi * R^2));     % Capacity function, defines the formula inside integral, SINR_db converted to linear
% eta_c = integral(@(D) f(D,1000,SINR_min_db,gamma), 1, R);                          % R is in meters (by default consider R = 1000 m 
% SINRe = 2^(eta_c/alpha) - 1;                                                 % SINR through Shannon Capacity Reu = alpha * Bw * log2(1 + SINR) = eta_c * Bw
% SINRe_db = 10 * log10(SINRe);

% SIR
SIRe = (1/ro) * (1 / (((6/omega) * sqrt(3)^(- gamma)) + (6 * sqrt(3 * N)^(-gamma))));         % SIR for SFR systems at the cell edge (d=R -> d/R = 1), and considering a cell load
Im_db = 10 * log10(inv(1 - SINRe / SIRe));                                  % Interference margin in dB

S_dbm = SINRe_db + N_dbm + Im_db;                                          % Minimum required Signal Power at the cell edge Smin FOR SFR
M_db = fzero(@(x) 0.5*erfc(-x/(sigma_db*sqrt(2)))-F, 0.5);
Pr_dbm = S_dbm + M_db; 

% R = 10 * ((Pr/P0)^(-1/gamma))                                               % linear, maximum cell radius using received power (Pr​) and reference power (P0)
Rm = 10 * 10^((Pr_dbm - P0_dbm) / (-10 * gamma));                             % in dB, maximum cell radius using received power (Pr​) and reference power (P0)

fprintf('maximum cell radius: %.2f \n', Rm);
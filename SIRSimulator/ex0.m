% find Prx ?

clear all;
N=1;
ro=0.7;
s = 3;
delta = 10e6; % bit/s/ km2 
sigma_db = 6;
F = 0.95;
NF_db = 3.5;
kT0_dbm = -174;
gamma = 4;
Btot = 5e6;    
R = 0.6;
alpha = 0.75;
Bw_eu = 180e3;

Ao = delta * pi * R^2;
eta_c_req = Ao / ro / Btot;
SINR_shan = 2^(eta_c_req/alpha) - 1;
fprintf('SINR Shannon: %.2f \n', SINR_shan);

SINR = fzero(@(x) eta_c1(R*1000, x, gamma) - eta_c_req, 3);
fprintf('SINR from fzero(): %.2f \n', SINR);

% SINR = 10^(SINR_db/10);

SINR_db = 10*log10(SINR);
N_dbm = kT0_dbm + 10*log10(Bw_eu) + NF_db;
SIR = 1/6 * 1/ro * s * sqrt(3*N)^gamma;
Im = inv(1-SINR/SIR);
Im_db = 10*log10(Im);
M_db = fzero(@(x) F-0.5*erfc(-x/(sigma_db*sqrt(2))),10);
S_dbm = SINR_db + N_dbm + Im_db;
Prx_dbm = S_dbm + M_db


% Output:
% SINR Shannon: 18.81 
% SINR from fzero(): 2.36 
% 
% Prx_dbm =
% 
%  -102.3498

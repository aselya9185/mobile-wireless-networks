% Compute the SINR at the cell edge in order to obtain an average cell capacity equal to 15 Mbit/s in case of:
% -  Capacity evaluation through Shannon formula with attenuation factor alpha =0.75
% -  Path Loss exponent = 4
% -  Total System Bandwidth Bw = 5 MHz
% -  Cell radius R=1000
% ----------------------------------------------------------------------------------------
% Solution SINR_db=2.6753 dB
 							
clear all;
alpha = 0.75;
gamma = 4; 
Bw_tot = 5e6; 
Ccell = 15e6; 	% bit rate / Cu / Rcell
R = 1000;

% SINR - ?
% Required Spectral Efficiency
eta_c_req = Ccell / Bw_tot;
SINR_req = 2^(eta_c_req/alpha) - 1;                                         % SINR through Shannon Capacity Reu = alpha * Bw * log2(1 + SINR) = eta_c * Bw
SINR_req_db = 10 * log10(SINR_req);
fprintf('eta_c_req: %.2f \n', eta_c_req);    % 3.000000
fprintf('SINR_req: %f \n', SINR_req_db);     % 11.760913 

% Actual Average Spectral Efficiency
f = @(D,R,SINR_db,gamma) 0.75*log2(1+(((D./R).^- gamma).*10^(SINR_db/10))).*(2*pi*D./(pi*R^2)); %Shannon as a function of SINR at edge 
eta_c=@(SINR_db,gamma) integral(@(D) f(D,R,SINR_db,gamma),1,R); % function of average spectral efficiency 

% Required SINR at the cell edge, such that the Actual Average Spectral Efficiency (eta_c) matches the Target Spectral Efficiency (eta_c_req).
SINR_db = fzero(@(SINR_db) eta_c(SINR_db,gamma) - eta_c_req,3);
% fprintf('eta_c: %.2f', eta_c);
fprintf('SINR: %f', SINR_db);                 % 2.675303

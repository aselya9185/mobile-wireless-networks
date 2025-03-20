% 4. Compute the maximum cell radius R (km) to obtain a cell load of ro=0.5, assuming that:
% -  Uniform distribution of user in the cell
% -  Capacity evaluation through Shannon formula with attenuation factor alpha =0.75
% -  SINR at cell edge SINR_db = 4 dB;
% -  System Bandwidth Bw_tot = 4.5 MHz
% -  Average traffic demand per km2 delta = 1 Mbit/s/ km2
% -  Path Loss exponent gamma = 4
% 
% ------------------------------------------------ ----------------------------------------------
% 
% Solution R=1.5350 km

% Rmax - ?

clear all;
ro = 0.5;
alpha = 0.75;
SINR_db = 4;
Bw_tot = 4.5e6;
delta = 1e6; %(bit/s/km^2)
gamma = 4;

% eta_c = @(R) integral(@(D) alpha * log2(1 + (((D./R).^-gamma).*10^(SINR_db/10))).*(2 * pi * D./(pi * R^2)), 1, R);                          % R is in meters (by default consider R = 1000 m
% R = fzero(@(R) (sqrt(eta_c(R) * Bw_tot * ro / (pi * delta)) - R), 1500);

f = @(D,R,SINR_db,gamma) 0.75 * log2(1+(((D./R).^- gamma).*10^(SINR_db/10))).* (2*pi*D./(pi*R^2));
eta_c = integral(@(D) f(D,1800,SINR_db,gamma),1, 1800)
R = sqrt(ro * eta_c * Bw_tot / (pi * delta)); % km, because delta per km2

fprintf('radius %f', R); % radius 1.534950
% Find:
% Radius R


clear all;
SINR_db = 4;
ro = 0.5;
delta_t = 1e6; % Mbit/s/km^2
Bcell = 4.6e6; % Mbit/s
gamma = 4;


f = @(D,R,SINR_db,gamma) 0.75*log2(1+(((D./R).^-gamma).*10^(SINR_db/10))).*(2*pi*D./(pi*R^2));
eta_c = integral(@(D) f(D,1000,SINR_db,gamma), 1,1000);

R = sqrt(ro*eta_c*Bcell / (delta_t*pi))

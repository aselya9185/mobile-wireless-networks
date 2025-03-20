% from 7 a.m. until 10 p.m. (15 hours in tot) br_eu1 = 20e6, the rest of
% the hours (9 hours) br_eu2 = 1e6

% Find:
% 1) Ptx_bs_1 and Ptx_bs_2
% 2) Ptx_day


Bw_tot = 18e6;
Ptx_max = 80; %W
R = 500; %m
delta = 1000; %user per km^2
br_eu1 = 20e6; % bytes per hour
br_eu2 = 1e6; %bytes per hour
SINR_db = 4;
gamma = 4;
% Ptx_day - ? Average transmit power in watt during the whole day

% bit per sec per m2 delta
delta1 = (br_eu1 * 8 / 3600) * delta / 1e6;
delta2 = (br_eu2 * 8 / 3600) * delta / 1e6;

f = @(D,R,SINR_db,gamma) 0.75*log2(1+(((D./R).^-gamma).*10^(SINR_db/10))).*(2*pi*D./(pi*R^2));
eta_c=integral(@(D) f(D,1000,SINR_db,gamma),1,1000)

ro1 = ((delta1*pi*R^2)/Bw_tot)/(eta_c)
ro2 = ((delta2*pi*R^2)/Bw_tot)/(eta_c)

Ptx_bs_1 = Ptx_max * ro1
Ptx_bs_2 = Ptx_max * ro2

Ptx_day = 15*(Ptx_bs_1)/24 + 9*(Ptx_bs_2)/24

% In case of an OFDMA cellular system, compute the cell radius R for 
% reuse factor N=1,3,7 assuming:
% - Full Frequency Reuse
% - Cell load ro=0.5
% - Sectoring = 3
% - SIR evaluation through approximated downlink formula (no simulator)
% - Large scale fading margin MdB = 4 dB; (note this is the required
% margin and not the standard deviation of the lognormal)
% - Path Loss exponent gamma= 4
% - Path loss at 1km = 125.13 dB
% - Noise figure NF =3.5 dB
% - Noise spectral density kT0 = -174 dBm Hz-1
% - Bandwidth for single users 180 kHz
% - No power control
% - Transmit Power Ptx_dbm= 23 dBm
% - SINR = 2.4 dB
% 
% Solution R (km) = [1.6298 1.7103 1.7173]

clear all;
M_db = 4;
gamma = 4;
NF_db = 3.5;
kt0_dbm = -174;
L0_db = 125.13;
Bw_eu = 180e3
N_v = [1 3 7];
Pt_dbm = 23;
SINR_db = 2.4;
R0 = 1e3;
ro = 0.5;
sec = 3;

for i = 1:3
    N = N_v(i);
    SIR = (1/ro) * sec * (1/6) * (sqrt(3 * N)^(gamma)); % SIR at cell edge from approx. downlink formula 
    SINR = 10^(SINR_db/10);
    Im_db=10 * log10(inv(1 - SINR/SIR)); % Interference margin
    Bw_eu_db = 10 * log10(Bw_eu)
    N_dbm = kt0_dbm + Bw_eu_db + NF_db;
    S_dbm = SINR_db + N_dbm + Im_db
    %required signal power at cell edge without fading margin
    Pr_dbm = S_dbm + M_db; %required signal power at cell edge with fading margin
    L_max = Pt_dbm-Pr_dbm;

    L0 = 10^(L0_db/10);

    L = @(R) L0 * (R/R0)^gamma;
    L_db = @(R) 10*log10(L(R)) % path loss
    res(i) = fzero(@(R) (L_db(R) - L_max),[1 5000])
end

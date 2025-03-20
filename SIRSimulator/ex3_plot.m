% Compute and plot the cell radius of an OFDMA system versus values of coverage probability F

clear all;

N = 1;
sector = 3;
% SIR evaluation through formula (no simulator)
alpha = 0.75;       % Shannon attenuation factor
sigma_db = 4;       % Large scale fading
gamma = 4;          % Path Loss exponent
L0_db = 141.44;     % path loss at 1 Km
R0 = 1000;          % m
NF_db = 3.5;        % Noise figure
kt0_dbm = -174;     % Noise spectral density
Br_eu = 500e3;      % bit rate for edge user [bit/s]
scs = 15e3;         % Subcarrier spacing Hz
subcarriers = 48;
Ptx_dbm = 18;       % BS tx power per subcarrier
load = 0.7;         % Cell Load
pcu=0;              % Uplink power control 
Bw_eu = 720e3       % Bandwidth for edge user

F_vector=[ 0.9 0.95 0.98] ; % Coverage probability


% solution
Bw_db = 10 * log10(Bw_eu);      %linear to dB


for i=1:length(F_vector)
   F = F_vector(i)

    SINR = 2 ^ (Br_eu/(alpha * Bw_eu)) -1
    
    N_dbm=kt0_dbm+Bw_db+NF_db %noise, bandwidth of user(if in question we have BW-tot & BW-eu), Noise figure, Noise spectral density
    
    SIR = (1/load)*sector*(1/6)*sqrt(3*N)^gamma;%sir formula, sectoring, cell load, reuse factor,Path Loss exponent
    
    SINR_db = 10*log10(SINR)
    Im_db=10*log10(inv(1-SINR/SIR)); % interference margin
    S_dbm = SINR_db+N_dbm+Im_db; %signal , fading margin
    M_db = fsolve(@(x) 0.5*erfc(-x/(sigma_db*sqrt(2)))-F,0.5)
    Pr_dbm = S_dbm + M_db
    
    
    Lmax_db = Ptx_dbm - Pr_dbm
    Lmax=10^(Lmax_db/10);
    L0=10^(L0_db/10);
    Rm(i)=R0*((Lmax/L0)^(1/gamma)) %all in linear
    Rm(i)/1000
    fprintf('SIR: %.2f \n', SIR);
    fprintf('SINR: %.2f \n', SINR);

end


plot(F_vector,Rm)

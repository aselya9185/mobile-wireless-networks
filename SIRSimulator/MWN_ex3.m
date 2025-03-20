% 3. Compute the maximum cell radius R (km), compute the Ptx transmission power (W) assuming that:
% -  Shannon attenuation factor =0.75
% -  Large scale fading dB = 4 dB;
% -  Coverage probability F = 0.98;
% -  Path Loss exponent gamma = 4
% -  Path loss at 1km = 125.13 dB
% -  Noise figure NF =3.5 dB
% -  Noise spectral density kT0 = -174 dBm Hz-1
% -  System Bandwidth Bw = 4.5 MHz
% -  Reuse factor 1
% -  Sectoring = 3
% -  No power control
% -  SINR at for user at cell edge = 4dB
% -  Bandwidth for single users 180 kHz
% -  Cell load ro = 0.7
% -  Average traffic demand per km^2 = 1 Mbit/s
% Plot R (km) and Ptx versus Cell load (0.4:0.1:0.8)
% ----------------------------------------------------------------------------------------------
% Solution R=1.86 km, Pt=1.55 W 



clear all

alpha = 0.75;
sigma_db = 4; 
F = 0.98;
gamma = 4; 
NF_db = 3.5; 
kt0_dbm = -174; 
R0 = 1; %km 
L0_db = 125.13; 
SINR_db = 4; 
Bw_tot = 300*15e3; 
Bw_eu = 180e3;
ro = 0.7;
delta = 1e6; %(bit/s/km^2)
sec = 3;

% search average spectral efficiency eta_c (bit/s/Hz)
f = @(D,R,SINR_db,gamma) 0.75*log2(1+(((D./R).^-gamma).*10^(SINR_db/10))).*(2*pi*D./(pi*R^2)); 
eta_c=integral(@(D) f(D,1000,SINR_db,gamma),1,1000) 

rop=0.4:0.1:0.8; 
for i=1:length(rop) 
    ro = rop(i); 						
    
    % Compute the cell radius
    R = sqrt(ro*eta_c*Bw_tot/(pi*delta));  %km
    
    % search for rx power at cell edge 
    %Pr_dbm=S_dbm+M_db
    %S = SINR+N+Im
    % Im = 1/(1-SINR/SIR)

    % Compute the SIR
    SIR = (1/ro)*3*(1/6)*sqrt(3)^gamma; 
    SINR = 10^(SINR_db/10);
    Im = 1/(1-SINR/SIR); 
    Im_db=10*log10(Im); 

    % Compute the noise power
    Bw_eu_db=10*log10(Bw_eu); 
    N_dbm=kt0_dbm+Bw_eu_db+NF_db; 

    % Compute the required signal power without fading margin
    S_dbm = SINR_db+N_dbm+Im_db; 
   
    % search for shadowing margin
    M_db = fzero(@(x) 0.5*erfc(-x/(sigma_db*sqrt(2)))-F,0.5); 
    
    % Compute the received power at the cell edge
    Pr_dbm = S_dbm + M_db; 
    %Pr_dbm = Pt_dbm-L_db(R) 
   
    % Compute the path loss
    L0=10^(L0_db/10);
    L= L0*(R/R0)^gamma;
    L_db = 10*log10(L);
    
    % Compute the transmission power
    Pt_dbm = Pr_dbm + L_db;
    Pt = 10^(Pt_dbm/10)/1000; % watt
    
    % Store the results for plotting
    Ptp(i)=Pt;
    Rp(i)=R;

    % Print the results for R and Pt
    fprintf('For cell load ro = %.1f:\n', ro);
    fprintf('Cell radius R = %.4f km\n', R);
    fprintf('Transmission power Pt = %.4f W\n\n', Pt);
end

% Plot the results
figure 
plot(rop,Rp); 
figure 
plot(rop,Ptp); 			
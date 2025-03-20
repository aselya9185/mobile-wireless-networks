% 2. Compute maximum uplink bit rate Bup achievable by a user at the cell edge assuming:
% The use of Shannon attenuation factor =0.75 for bit rate computation;
% Cell radius R = 5 km
% 	Large scale fading dB = 4 dB; 
% Coverage probability F = 0.90; 
% Path Loss exponent = 4 
% Noise figure NF = 2 dB
% Noise spectral density kT0 = -174 dBm Hz-1
% Assigned Bandwidth Bw = 180 kHz 
% Prx(10m) = 100 mW
% 	Reuse Factor = 1
% No sectoring
% Uplink with power control
% Use simulator for SIR computation 
% Cell load = 0.7		
% ----------------------------------------------------------------------------------------------		
% Solution Bup about equal to 300 kbit/s
% compute S (required signal power without fading margin)
% find SIRup
% find SNR
% find SNRup
% Find Bup

% bit rate Bup = alpha * 180000 * log2(1+SINRup)
% SINRup= inv(inv(SNR)+inv(SIRup)) 
% SIR from simulator
% SNR = S/N	
% S_dbm = Pr_dmb - M_db = average signal power received at the cell edge without fading margin (i.e. worst case in which fading margin M is completely consumed to cope with shadowing)
% Pr computed by link budget: Pr_dbm = P0_dbm-gamma*10*log10(Rm/10)



alpha = 0.75;            
sigma_db = 4;             
F = 0.90;                 
gamma = 4;                
NF_db = 2;                
kt0_dbm = -174;           
P0_dbm = 10 * log10(100); 
Bw = 180000;              
R = 5000;                 
load = 0.7;               
Rue = R;                  
pcu = 1;                  

% Noise power
Bw_db = 10 * log10(Bw);        
N_dbm = kt0_dbm + Bw_db + NF_db;

% Shadowing margin calculation
M_db = fsolve(@(x) 0.5 * erfc(-x / (sigma_db * sqrt(2))) - F, 0.5);

% Signal power at cell edge
Pr_dbm = P0_dbm - gamma * 10 * log10(R / 10); 
S_dbm = Pr_dbm - M_db;                        

% SIR computation for uplink
[SIRup, SIRdown] = InterferenceComputation(R, gamma, load, 1, 1, Rue, pcu);

% SNR and SINR calculation
S = 10^(S_dbm / 10);        
N = 10^(N_dbm / 10);        
SNR = S / N;                
SINRup = inv(inv(SNR) + inv(SIRup)); 

% Maximum uplink bit rate
Bup = alpha * 180000 * log2(1 + SINRup); 
disp(['Maximum Uplink Bit Rate: ', num2str(Bup), ' bps']);
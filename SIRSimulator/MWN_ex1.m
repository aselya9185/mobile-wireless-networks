% From MWN-Exercises.pdf
% 
% cell Radius questions: 1, 3, 4
% Ptx power: 3, 5
% 
% ex 1. Compute maximum cell radius R (Km) assuming:
% Large scale fading dB = 6 dB; 
% Coverage probability, F = 0.95; 
% Path Loss exponent = 4;
% Noise figure NF = 3.5 dB;
% Noise spectral density kT0 = -174 dBm Hz-1;
% 
% Assigned Bandwidth Bw = 180 kHz; 
% Prx(10m) = 0.1 W;
% SINR at cell edge = -2dB;
% SIR at cell edge = 11.3 dB;		
% ----------------------------------------------------------------------------------------------	
% Solution Rm = 17.65 km	
% 
% required signal at cell edge without fading: S_dbm = SINR_db+N_dbm+Im_db 
% required signal at cell edge with fading margin Pr_dbm = S_dbm + M_db
% link budget db: Pr_dbm = P0_dbm- gamma*10*log10(Rm/10)
% link budget linear: Pr = P0*(Rm/10)^(- gamma) --> Rm=10*((Pr/P0)^(-1/gamma))	

clear all
sigma_db = 6;         % Standard deviation of fading
F = 0.95;             % Coverage probability
gamma = 4;            % Path loss exponent
NF_db = 3.5;          % Noise figure in dB
kt0_dbm = -174;       % Noise spectral density in dBm/Hz
P0_dbm = 10*log10(0.1*1000); % Transmit power at 10m in dBm (converted from 0.1 W)
SINR_db = -2;         % SINR at cell edge in dB
SIR_db = 11.3;        % SIR at cell edge in dB
Bw = 180000;          % Bandwidth in Hz

% Calculate bandwidth in dB
Bw_db = 10*log10(Bw);

% Calculate noise power in dBm
N_dbm = kt0_dbm + Bw_db + NF_db;

% Convert SINR and SIR to linear scale
SINR = 10^(SINR_db/10);
SIR = 10^(SIR_db/10);

% Calculate interference margin in dB
Im_db = 10*log10(inv(1 - SINR/SIR));

% Calculate required signal power at cell edge without fading in dBm
S_dbm = SINR_db + N_dbm + Im_db;

% Search for shadowing margin (fading margin)
M_db = fsolve(@(x) 0.5*erfc(- x/(sigma_db*sqrt(2))) - F, 0.5);

% Required signal power at cell edge with fading margin
Pr_dbm = S_dbm + M_db;

% Convert power from dBm to milliwatts
Pr = 10^(Pr_dbm/10);    % in milliwatts
P0 = 0.1*1000;          % Transmit power in milliwatts (converted from 0.1 W)

% Calculate the maximum cell radius in meters
Rm = 10 * ((Pr / P0)^(-1/gamma));

% Convert radius to kilometers
R = Rm / 1000;

% Display the maximum radius
disp(['Maximum Cell Radius: ', num2str(R), ' km']);
% Find:
% max cell radius R true for conditions:
% 1) Received power without margin Srx_dbm >= -100 dBm (receiver sensitivity)
% 2) SINR_db >= 0 dB

clear all;
gamma=4;                % Path loss exponent
ro=0.7;
sec=3;
Ptx_dbm=26;             % Transmission power over Bw_eu
G_db = 8;               % Base station TX/RX antenna gain
Ssens_dbm = -100;       % Receiver sensitivity
SINR_min_db = 0;        % Minimum SINR
sigma_db = 6;           % Lage scale fading
%SM_db= 6;
F=0.95;                 % Coverage probability
NF_db=3.5;              % Noise figure
kt0_dbm=-174;           % Noise spectral density
R0=1000;                % Reference distance in m
L0_db=125.13;           % Path loss at 1 km
Bw_eu = 180e3;          % Bandwidth fpr the edge user Hz
N=1;

%M_db=SM_db(sigma_db,F);%
M_db = fsolve(@(x) 0.5*erfc(-x/(sigma_db*sqrt(2)))-F,0.5); %FADING MARGIN

% compute cell radius in case of sensitivity bound
Srx_dbm = @(R) Ptx_dbm + G_db - M_db - (L0_db + gamma * 10 * log10(R/R0));

% Uses fzero() to find the largest R where the received power equals the receiver sensitivity (−100 dBm).
Rsens = fzero(@(R) Srx_dbm(R)-Ssens_dbm,1000) % sensitivity radius

% SNR
SNR_db = @(R) Srx_dbm(R) - (kt0_dbm+10*log10(Bw_eu)+NF_db);
SNR = @(R) 10^(SNR_db(R)/10);

%SIR
SIR = (1/ro)*sec/6*sqrt(3)^gamma;

%SINR
SINR = @(R) inv(inv(SNR(R))+inv(SIR));
SINR_db = @(R) 10*log10(SINR(R));
SINR_res=SINR_db(Rsens); % 7.6549 --> SINR at Rsens > SINR_min, Rsens max cell radius
fprintf('Rsens: %.2f m\n\n', Rsens);
fprintf('SINR: %.2f \n\n', SINR_res);

% Output:
% Rsens: 944.11 m
% SINR: 7.65 
% Compute and plot the bit-rate provided to an edge user versus the number of sectors (1 3 6) in case of:

clear all;
N = 1;
% SIR evaluation through formula (no simulator)
alpha = 0.75;       % Capacity evaluation through Shannon formula
sigma_db = 4;       % Lage scale fading
F = 0.95;           % Coverage probability
gamma = 4;          % Path Loss exponent
L0_db = 141.43;     % Path Loss at R0 = 1km
R0 = 1000;          % Reference distance in meters
NF_db = 3.5;        % Noise Figure
kt0_dbm = -174;     % Noise spectral density
Bw_eu = 360e3;      % Bandwidth for edge user
% No power control
R = 500;            % Cell radius in meters
Pt_db = 26;         % Base Station Transmit Power
G_db = 8;           % Antenna gain (TX/RX)
sectors = [1 3 6];
% No antenna gain for mobile phones
ro = 0.7

Reu_values = zeros(size(sectors));  % Array to store edge user bit rates
% SNR
Bw_eu_db = 10 * log10(Bw_eu);
    
N_dbm = kt0_dbm + Bw_eu_db + NF_db;                                         % Noise Power as a combination of Noise Spectral Density kt0, Bandwidth Bw_eu + Noise Power F
    
M_db = fzero(@(x) 0.5*erfc(-x/(sigma_db*sqrt(2)))-F,0.5);                   % Shadowing margin (fading margin) to have a coverage probability F,, fsolve also finds the root of the equation g(x)=0, but it is designed for solving systems of nonlinear equations.
L_db = L0_db + 10 * gamma * log10(R/R0);                                    % Log-Distance Path Loss Model, R0 - reference distance for which the L0_dB is known
S_dbm = Pt_db + G_db - L_db - M_db;                                         % Minimum required Signal Power at the cell edge from Link Budget L_db = Pt_db + G_db - S_dbm - M; 
SNR_db = S_dbm - N_dbm; 
SNR = 10^(SNR_db/10);
fprintf('SNR: %.2f \n', SNR);

for i = 1:length(sectors)
    sec = sectors(i);
    fprintf('sector = %i\n', sec);

    % SIR
    SIR = (1/ro) * sec * (1/6) * sqrt(3 * N)^gamma;
    fprintf('SIR: %.2f \n', SIR);
    
    % SINR
    SINR = inv(inv(SNR)+inv(SIR));                                              % Compute SINR based on SNR and SIR
    fprintf('SINR: %.2f \n', SINR);
    
    % Shannon capacity
    Reu_values(i) = alpha * Bw_eu * log2(1 + SINR);                                       % The bit rate for a user at the cell edge
    fprintf('Reu bit-rate provided to an edge user: %.2f bps\n\n', Reu_values(i));

end

% Plot bit-rate vs. cell load (ro)
figure;
plot(sectors, Reu_values);
xlabel('Sectors');
ylabel('Bit-rate at Edge User (bps)');
title('Bit-rate vs. Sectors');
grid on;

% Output:
% ro = 0.7000
% 
% SNR: 19.81 
% sector = 1
% SIR: 2.14 
% SINR: 1.93 
% Reu bit-rate provided to an edge user: 419232.83 bps
% 
% sector = 3
% SIR: 6.43 
% SINR: 4.85 
% Reu bit-rate provided to an edge user: 688312.10 bps
% 
% sector = 6
% SIR: 12.86 
% SINR: 7.80 
% Reu bit-rate provided to an edge user: 846982.75 bps
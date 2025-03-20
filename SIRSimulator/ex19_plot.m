% Compute and plot the required System Bandwidth Bw_tot for serving a cell 
% with radius R=1.5 Km as a function of number of sectors (1,3,6) and assuming that:
% - Uplink power control
% - Reuse factor N=1
% - Cell load ro=0.7
% - Uniform distribution of users in the cell
% - Capacity evaluation through Shannon formula with attenuation factor alpha =0.75
% - Signal to Noise ratio (SNR) at cell edge = 13 dB
% - SIR at cell edge to be computed with simulator considering the uplink value;
% - Average traffic demand per km2 delta = 1 Mbit/s/ km2
% - Path Loss exponent gamma = 4


clear all
N = 1;                 % Reuse factor
ro = 0.7;              % Cell load
alpha = 0.75;          % Attenuation factor
SNR_db = 13;           % SNR in dB at cell edge
delta = 1;             % Traffic demand per km^2
s_v = [1 3 6];         % Number of sectors to evaluate
R = 1.5;               % Radius of the cell in km
gamma = 4;             % Path Loss exponent

% Step 1: Calculate cell area
Ao = delta * pi * R^2;

% Position of the user at the cell edge
Reu = R;               % distance of the reference UE from its BS
pcu = 1;               % Power control on uplink
SNR = 10^(SNR_db / 10); % Convert SNR to linear scale

Btot = [];             % Initialize array for total bandwidth

for i = 1:length(s_v)
    s = s_v(i);        % Number of sectors for current iteration
    
    % Step 2: Compute SIR using the function 'InterferenceComputation'
    [SIRu, SIRd] = InterferenceComputation(R * 1000, gamma, ro, N, s, Reu * 1000, pcu);
    
    % Step 3: Compute SINR based on SNR and SIR
    SINR = inv(inv(SNR) + inv(SIRu));
    
    % Step 4: Define the capacity function and integrate over the cell area
    f = @(D, R, SINR, gamma) 0.75 * log2(1 + (((D ./ R) .^ -gamma) .* SINR)) .* (2 * pi * D ./ (pi * R^2));
    
    eta_c = integral(@(D) f(D, R * 1000, SINR, gamma), 1, R * 1000);
    
    % Step 5: Compute total bandwidth for the current sector configuration
    Btot(i) = (Ao / ro) / eta_c;

    % Print the result for each sector value
    fprintf('System Bandwidth for %d sectors: %f Hz\n', s, Btot(i));
end

% Step 6: Plot the results
figure
bar(s_v, Btot);
xlabel('Number of Sectors');
ylabel('Total System Bandwidth Btot (Hz)');
title('Required System Bandwidth vs Number of Sectors');


% OUTPUT
% uplink S/I 5.848406 dB, downlink S/I -3.292635 dB
% System Bandwidth for 1 sectors: 2.857982 Hz
% uplink S/I 10.989962 dB, downlink S/I -0.236771 dB
% System Bandwidth for 3 sectors: 2.285195 Hz
% uplink S/I 13.515525 dB, downlink S/I Inf dB
% System Bandwidth for 6 sectors: 2.126595 Hz


% DESCRIPTION OF THE OUTPUT
% Uplink refers to the communication from the user equipment (UE) to the base station.
% S/I = 5.85 dB indicates that the signal received at the base station is stronger than the interference by about 5.85 dB
% Positive S/I values (like 5.85 dB) suggest that the signal strength is higher than the interference, which is generally good for communication quality.

% Downlink refers to the communication from the base station to the user equipment (UE).
% An S/I ratio of -3.292635 dB indicates that the signal is weaker than the interference by about 3.29 dB.
% Negative S/I values (like -3.29 dB) suggest that the interference is stronger than the signal, which can degrade communication quality, leading to potential errors or interruptions in the downlink communication.
% Interpretation:
% In this case, the uplink communication (from the user to the base station) is better (with an S/I ratio of 5.85 dB), while the downlink communication (from the base station to the user) is suffering from interference (with an S/I ratio of -3.29 dB).
% To improve downlink performance, you would need to either reduce interference or boost the signal strength.
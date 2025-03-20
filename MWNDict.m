N                      % Reuse factor, if the Reuse factor is 1, meaning every cell in the network uses the same frequency band. B_tot = B_cell
ro                     % Cell load, a cell load of 0.7 means that 70% of the cell s capacity is being used.
alpha = 0.75;          % Attenuation factor, A higher attenuation factor indicates more signal loss, reducing capacity and coverage.
                       % Example: If the attenuation factor is 0.75, it means that 75% of the theoretical maximum capacity is available, accounting for signal loss.
delta                  % User density [user/ km^2]. Traffic demand per km^2, meaning each user consumes 1 Mbit/s.
s_v                    % Number of sectors to evaluate
R                      % Radius of the cell in km
gamma                  % Path Loss exponent, Path loss is the reduction in signal strength as a signal travels from the transmitter to the receiver.
omega                  % Ptx_e/Ptx_c    external/central Tx power ratio
Reu                    % distance of the reference UE from its BS
pcu                    % Power control on uplink, Devices farther from the base station increase their transmission power, while those closer reduce it.
sigma_db               % Standard deviation of fading, long-term variations in the received signal power
F                      % Coverage probability, 95% of the users in the coverage area will have a good signal, accounting for fading and other variations in signal strength.
NF_db                  % Noise figure in dB, The noise figure measures how much the receiver amplifies the noise along with the signal.
                       % A low noise figure improves the receiver's sensitivity, allowing it to detect weaker signals.
kt0_dbm = -174;        % Noise spectral density in dBm/Hz, This parameter sets the baseline for the noise present in any system. Combined with bandwidth, it helps to calculate the total noise power in the system.
P0_dbm                 % Transmit power at 10m in dBm (converted from 0.1 W). This is a reference point for computing how much power is received over distance using the path loss formula. 
SNR_db                 % SINR at cell edge in dB
SIR_db                 % SIR at cell edge in dB
SNR_db                 % SNR in dB at cell edge
Bw                     % Bandwidth in Hz
Bw_eu                  % User bandwidth 
Br_ue                  % User bit rate 250 kbit/sec
Reu                    % The bit rate for a user at the cell edge
Rue                    % User distance from the base station
S_dbm                  % Required signal power at cell edge
Srx_dbm                % represents the actual received power after signal degradation due to path loss, shadowing, etc. (Prx_dbm - M_db)
Pr_dbm / Prx_dbm       % Required signal power at cell edge with fading margin  (S_dbm + M_db)
N_dbm                  % Noise Power in the bandwidth
M_db                   % Shadowing Margin, provides extra signal strength to maintain coverage despite obstacles that might block or reflect the signal.
Im_db                  % Interference Margin
L0                     % Path Loss
L                      % final path loss at the maximum distance (R) from the base station.
R0                     % Radius for path loss formulas
cap                    % capacity of the channel


% Convertings
% dB-linear & linear-dB:
y = 10^(y_db/10);
x_db = 10 * log10(x);   


% ---------------------------------------------------------- WIRELESS TRANSMISSION ----------------------------------------------------------
    % SPECTRAL EFFICIENCY is the actual achievable data rate per unit of bandwidth, considering all system factors (modulation, coding, interference, noise, etc.). Example, QPSK + FEC:
        % η_c = C/B = [Channel capacity] / [Bandwidth], [bps/Hz]
    % It is computed using the Shannon formula adjusted with the attenuation factor.
eta_c = cap / Bw_eu;                                                        % Spectral efficiency (ηc), From the Shannon capacity formula, [Bit rate / bandwidth] (bps/Hz)
eta_c = alpha * log2(1 + SINR);                                             % Spectral efficiency (ηc) through SINR
eta_c = alpha * log2(1 + 10^(SINR_db / 10));                                % Spectral efficiency (ηc) when SINR is given in dB
eta_c = alpha * log2(1 + SINR_min);                                         % Spectral efficiency (ηc) if no other SINR values are given, assuming SINR_min ensures the system meets the minimum performance requirement.
cap = eta_c * Bw_tot / N                                                    % capacity evaluation


    % Theoretical limit on data rate - SHANNON CAPACITY formula (Shannon-Hartley:): C = B * log2(1 + SINR)
    % Shannon capcity Real World (η - attenuation factor): C = η * B * log2(1 + SINR)       

    % Alpha-Shannon capacity (with attenuation): C = alpha * B * log2(1 + SINR)
Reu = alpha * Bw_eu * log2(1 + SINR);                                       % The bit rate for a user at the cell edge
Bup = alpha * Bw_eu * log2(1 + SINRup);                                     % Attenuated Shannon capacity Uplink, maximum possible data rate a channel can achieve without errors

% FFR FRACTIONAL FREQUENCY REUSE
fe = @(d,R) alpha*log2(1+SINRe(d)).*(2*pi*d./Ae);                           % capacity function (for Spectral efficiency) at the edge
fc = @(d,R) alpha*log2(1+SINRc(d)).*(2*pi*d./Ac);                           % capacity function (for Spectral efficiency) in the center of the cell
Ccell = (Ae/Atot)*(Be)*eta_ce + (Ac/Atot)*(Bc)*eta_cc;                      % Average Cell Capacity

c_fe = @(d,R) Be.*alpha.*log2(1+SINRe(d)).*(2*pi*d./Atot);                  % Capacity Contribution from the Edge Region
c_fc = @(d,R) Bc.*alpha.*log2(1+SINRc(d)).*(2*pi*d./Atot);                  % Capacity Contribution from the Central Region
Ccell=integral(@(d) c_fe(d,R),Rc,R)+integral(@(d) c_fc(d,R),1,Rc);          % Average Cell Capacity by separately evaluating the contribution of both the edge region and the central region of the cell

    % Modulation Spectral Efficiency (QPSK, QAM...): 
        % η = Br_eu/Bw_eu  = [bit rate without coding] / [Bandwidth] = log2M/(1 + beta) 
            % M = 2^b modulation order, example: QPSK -> M = 4, 16WAM -> M = 16
eta_c = Br_eu / Bw_eu;                                                      % Spectral efficiency (ηc), From the Shannon capacity formula, [Bit rate / bandwidth] (bps/Hz)    

    % Receiver Sensitivity Prs is the minimum average received signal power Pr to satisfy a certain BER: Prs = (Eb/N0)min * N0 * Rb


% ******* PATH LOSS *******
    % Modeling Path Loss: L(d,t)_dB = PathLoss(d)_dB + sigma_dB + theta_dB [dB]
    % where L(d,t)_dB - Attenuation at distance d and time t can be modeled as a random variable
    %       PathLoss(d)_dB - Average value of the attenuation over space and time
    %       sigma_dB - Zero mean random variable modeling large scale fading
    %       theta_dB - Zero mean random variable modeling small scale fading
L_max = Pt_dbm - Pr_dbm;                                                    % Maximum Path Loss as a difference between received signal and transmitted signal
L_dbm = Pt_dbm - Pr_dbm;                                                    % Path Loss as 
Pt_dbm = Pr_dbm + L_dbm;                                                    % Required transmission power of a base station Pt in dbm
Pt = 10^(Pt_dbm/10)/1000;                                                   % Required transmission power of a base station Pt in watt
Ptx_ue_dbm = Prx_min_ue_dbm + L_db;                                         % Required transmission power of the user equipment Ptx_ue

L_db = L0_db + 10 * gamma * log10(Rue / R0);                                % Log-Distance Path Loss Model, R0 - reference distance for which the L0_dB is known
L = L0 * (R/R0)^gamma;                                                      % path loss as a linear ratio instead of dB (equvalent og Log-Distance Path Loss Model). final path loss at the maximum distance (R) from the base station.
L = @(R) L0*(R/R0)^gamma;                                                   % This defines L as a function of R, not a single value. You can later compute L for different values of R by calling L(some_R_value).

% Subcarriers
subcarriers = Bw_tot / Bw_sub;                                              % Number of subcarriers
Pt_dbm = Ptx_sub_dbm + 10*log10(subcarriers);                               % Total BS transmit power



% *********************************************************** RADIO NETWORK PLANNING ***********************************************************
    % Network Spectral Efficiency: more cells, more capacity
    % Network Capacity = frequency reuse * Ccell
    % Network Spectral Efficiency = (frequency reuse * Ccell) / Btot
ec = N * Ccell / Bw_tot;                                                    % Network Spectral Efficiency 
Ccell = ec * Bw_tot / N;                                                    % cell capacity
Nu_per_cell = Ccell / C_u * ro;                                             % Maximum number of users per cell


% ---------------------------------------------------------- COVERAGE PLANNING ----------------------------------------------------------
% ----------------------------------------------------- No interference, single cell ----------------------------------------------------

    % to find a Cell edge (cell radius): SINR >= SINRmin (0 dB), and Prx >= Prs (receiver sensitivity = -110 dB)
    % Minimum required Signal Power at the cell edge Smin via receiver sensitivity Prs CONSIDERING ONLY NOISE, WITHOUT INTERFERENCE (in SNR): 
        % Smin = Prs = SNRmin * N = X * N
N_dbm = kt0_dbm + Bw_eu_db + NF_db;                                         % Noise Power as a combination of Noise Spectral Density kt0, Bandwidth Bw_eu + Noise Power F
S_dbm = SNR_db + N_dbm;                                                     % Minimum required Signal Power at the cell edge

% ******* INTERFERENCE MARGIN *******
    % Minimum required Signal Power at the cell edge Smin via receiver sensitivity Prs CONSIDERING INTERFERENCE MARGIN Im (in SINR not SNR): 
        % Smin = SINRmin * (N + I) = X * (N + I) = Prs * (N + I/N) = Prs * Im ---> 
            % Smin_dBm = Prs_dBm + Im_dB 
S_dbm = Prs_dBm + Im_dB;                                                    % Minimum required Signal Power at the cell edge Smin
S_dbm = SINR_db + Im_db;                                                    % Minimum required Signal Power at the cell edge Smin

    % Minimum required Signal Power at the cell edge Smin via SINRmin CONSIDERING NOISE N and INTERFERENCE MARGIN Im: 
        % Smin = Beu[dB*Hz] + kT0[dBm*Hz^-1] + NF[dB] + Im[dB] + SINRmin[dB]
S_dbm = SINR_db + N_dbm + Im_db;                                            % Minimum required Signal Power at the cell edge Smin FOR SFR
S_dbm_up = SINR_db + N_dbm + Im_db_up;                                      % Minimum required Signal Power uplink
S_dbm_down = SINR_db + N_dbm + Im_db_down;                                  % Minimum required Signal Power uplink


% ******* FADING MARGIN *******
    % Received Signal Power at the cell edge Pr WITH FADING MARGIN M should be:
        % Received power Pr = Smin [Minimum required Signal Power at the cell edge] + M [Fading Margin]
        % Pr >= Smin -> SINR >= SINRmin
Pr_dbm = S_dbm + M_db;                                                      % Received signal power at the cell edge 
Pr_dbm = Pt_db - L_db;
S_dbm = Pr_dbm - M_db;                                                      % Minimum required Signal Power at the cell edge
S_dbm = Pt_db + G_db - L_db - M_db;                                         % Minimum required Signal Power at the cell edge from Link Budget L_db = Pt_db + G_db - S_dbm - M; 
S_dbm = Pt_db - L_db - M;                                                   % Minimum required Signal Power at the cell edge from Link Budget L_db = Pt_db + G_db - S_dbm - M; 

    % Fading Margin M provides a coverage probability F at cell edge, and generally at a given distance d.
    % Radio Coverage at cell edge is guaranteed is the Fading Margin M is greater than fading losses:
        % M > sigma [loss due to Large Scale Fading] + theta [loss due to Small Scale Fading]
    % To simplify:
        % M = ML [Large Scale Fade Margin] + MS [Small Scale Fade Margin], ML >= sigma, MS >= theta
        % F = 0.5 * erfc(-M/sigma*sqrt(2))
        % 0.5 * erfc(-M/sigma*sqrt(2)) - F = 0
M_db = fzero(@(x) F-0.5*erfc(-x/(sigma_db*sqrt(2))), 10);                   % Shadowing margin (fading margin) to have a coverage probability F, fzero finds the root of the function f(x)=0, for single-variable equations.
M_db = fsolve(@(x) 0.5*erfc(-x/(sigma_db*sqrt(2)))-F, 0.5);                 % Shadowing margin (fading margin) to have a coverage probability F,, fsolve also finds the root of the equation g(x)=0, but it is designed for solving systems of nonlinear equations.
M_db = fzero(@(x) 0.5*erfc(-x/(sigma_db*sqrt(2)))-F, 0.5)

        % Transmit signal power at reference distance R0 from the base station: Pt = Pt + G0 - L0 + M [dBm]
        % Received signal power at the cell edge model: Pr = Pt + G0 - L0 - Lp + M [dBm]



% ******* LINK BUDGET AND CELL RADIUS *******
    % LINK BUDGET. At the cell edge we need: Pr(Dmax)_dBm = Smin_dBm + M_L*_dB + M_S*_dB = Pt_dBm + G0_dB - L0_dB - L(Dmax)_dB
% RECEIVED POWER
Pr_dbm = P0_dbm - gamma * 10 * log10(R/10);                                 % Link budget in dB, received power at distance R, using a reference power P0
Pr = P0 * (R/10)^(-gamma)                                                   % Link budget in linear, received power at distance R, using a reference power P0
Prx_dbm = Ptx_dbm - L_db                                                    % received power as a difference of transmit power and path loss in dB

    % There are few ways of finding a CELL RADIUS Dmax or R: 
        % 1) Using a PATH LOSS
        % Firstly, find ALLOWED PATH LOSS from a LINK BUDGET
            % L(Dmax) = Pt + G0 - L0 - Smin - M_L* - M_S*
L_db = Pt_db + G_db - Pr_db;                                                % Allowed path loss from a link budget
L_db = Pt_db + G_db - S_dbm - M;                                            % Allowed path loss from a link budget, Pr = Smin + M

            % L(Dmax) = Pt + G0 - L0 - Pr(Dmax),    Pr(Dmax) = Smin + M_L* + M_S*
                % Small scale fading M_S* can be neglected because of Receiver Intelligence, which help recover the signal, or take into account a small margin (2.0-5.0 dB)
                % Large scale fading M_L* computed by using a log-normal characterization
L_max = Ptx_ue_dbm + G_ue_db + G_bs_db - Prx_min_bs_dbm - Im_db - M_db      % Allowed Path loss
L_max = Ptx_bs_dbm + G_ue_db + G_bs_db - Prx_min_ue_dbm - Im_db - M_db

        % Secondly, from Log-Distance Path Loss Model
            % Lp = Lp(d0) + 10 * gamma * log(d/d0) [dB]
% PATH LOSS
L_db = L0_db + 10 * gamma * log10(R/R0);                                    % Log-Distance Path Loss Model, R0 - reference distance for which the L0_dB is known
L = L0 * (R/R0)^gamma                                                       % path loss as a linear ratio instead of dB (equvalent of Log-Distance Path Loss Model). final path loss at the maximum distance (R) from the base station.
L = @(R) L0 * (R/R0)^gamma;                                                 % This defines L as a function of R, not a single value. You can later compute L for different values of R by calling L(some_R_value).

        % find a d (maximum distance or cell radius):
            % d = d0 * 10^((Lp - Lp(d0)) / (10 * gamma))
% MAXIMUM DISTANCE (CELL RADIUS)
R = R0 * 10^((L_db - L0_db) / (10 * gamma));                                % in dB, maximum cell radius from Log-Distance Path Loss Model
R = R0 * ((L/L0)^(1/gamma));                                                % linear, maximum cell radius from Log-Distance Path Loss Model

        % 2) Power-based approach instead of path loss:
R = 10 * ((Pr/P0)^(-1/gamma))                                               % linear, maximum cell radius using received power (Pr​) and reference power (P0)
Rm = 10 * 10^((Pr_db - P0_db) / (-10 * gamma));                             % in dB, maximum cell radius using received power (Pr​) and reference power (P0)

R = fzero(@(R) (L_db(R)-L_max),[1 5000]);

        % # In case of receiver sensitivity bound
Srx_dbm = @(R) Ptx_dbm + G_db - M_db - (L0_db + gamma * 10 * log10(R/R0));  % Compute the received power at distance R
% Uses fzero() to find the largest R where the received power equals the receiver sensitivity (−100 dBm).
Rsens = fzero(@(R) Srx_dbm(R)-Ssens_dbm,1000)                               % Sensitivity bound radius, maximum distance at which the receiver can still function.




% ---------------------------------------------------------- COVERAGE PLANNING ----------------------------------------------------------
% ----------------------------------------------------- Interference, multple cells -----------------------------------------------------

% ******* RESOURCE REUSE *******
    % INTRA-CELL INTERFERENCE: happens inside the same cell
    % INTER-CELL INTERFERENCE: the interference generated by transmitters of near cells
        % FIXED CHANNEL ALLOCATION (FCA) to reduce inter-cell interference
            % B_cell = Bw_tot / N [Reuse factor]
Bcell = Bw_tot / N                                                          % Bandwidth allocated per cell

        % Reuse distance:
            % D = R * sqrt(3N)

        % SIR - FCA - Downlink - Omnidirectional Antenna (radiates equally around 360)
            % SIRinter = S [signal power at the cell edge] / Iinter [inter-cell interference) = 1/Ni [number of interfering cells] * (sqrt(3N))^gamma
            % SIRinter = 1/6 * (sqrt(3N))^gamma     for hexagon shaped cell
SIR = (1/Ni) * sqrt(3 * N)^gamma;   
SIR = (1/6) * sqrt(3 * N)^gamma;   

        % SIR - FCA - Uplink - Omnidirectional Antenna (radiates equally around 360)
            % SIRinter = 1 / Ni * (sqrt(3N) - 1) ^ gamma
SIR = (1/Ni) * sqrt(3 * N - 1)^gamma;   
SIR = (1/6) * sqrt(3 * N - 1)^gamma;   

        % POWER CONTROL
            % No (transmit) power control: 
                    % Ptx = const = Prx_t(R)^gamma, R - distance to cell edge
            % With (transmit) power control: 
                    % Ptx = Prx_t(r)^gamma, r - distance to mobile phone located in any position
            % Received power from reference mobile (at the cell edge): 
                    % Prx = Ptx(r)^-gamma
            % Received power from interfering mobile (in another cell): 
                    % Pri = Ptx(Dother)^-gamma
            % SIRinter = Prx / sum(Pri)
            % Pri = 0 if interefering and reference stations use different channels, sectors, or there is no load.

        % SITE SECTORING - Downlink and Uplink
            % SIR_inter-sec = SIR_inter-omni * K,   K - number of sectors
SIR = K * (1/Ni) * sqrt(3 * N)^gamma;  
SIR = K * (1/6) * sqrt(3 * N)^gamma;  

        % CELL LOAD:
            % SIR_inter-load = SIR_inter-full * (1/ro)
        % SIR DOWNLINK improved by cell load, sectoring, and capacity reuse in different cells:
            % SIR = 1/ro * sec * (1/Ni * (sqrt(3*N))^gamma)
SIR = (1/ro) * K * (1/Ni) * sqrt(3 * N)^gamma;                              % SIR for N Frequency Reuse, Ni = N-1 number of interfering cells, K = sectors
SIR = (1/ro) * K * (1/6) * sqrt(3 * N)^gamma;                               % SIR for hexagon shaped cell Ni = 6, K = sectors, N = reuse factor
SIR = (1/ro) * K * (1/6) * sqrt(3)^gamma;                                   % SIR for hexagon shaped cell Ni = 6, with N = 1 (all cells use the same frequiency, K = sectors

        % SIR UPLINK improved by cell load, sectoring, and capacity reuse in different cells:
            % SIR = 1/ro * sec * (1/Ni * (sqrt(3*N)-1)^gamma)

[SIRup SIRdown]=InterferenceComputation(R,gamma,ro,N,sectors,Rue,pcu);    % SIR for both uplink and downlink, when pcu is given, using the function 'InterferenceComputation'

        % SIR and INTERFERENCE MARGIN
            % Im = SNR / SINRmin
            % Im = 1 / (1 - (SINRmin / SIRinter))
        % ALWAYS NEED TO BE CONSIDERED FOR SFR
Im = inv(1 - SINR / SIR);                                                   % Interference margin linear (inv(x) is equivalent to 1/x)
Im_db = 10 * log10(inv(1 - SINR / SIR));                                    % Interference margin in dB
Im_db_up = 10 * log10(inv(1 - SINR / SIRup));                               % Interference margin Uplink in dB
Im_db_down = 10 * log10(inv(1 - SINR / SIRdown));                           % Interference margin Downlink in dB



% ******* FRACTIONAL / SOFT FREQUENCY REUSE *******
    % FFR: The center of the cell and the edge use completely different frequency sets. If Reuse Factor = 3, then the edge frequencies are split into 3 groups, and each cell gets only 1 of the 3 groups.
        % Downlink / Uplink:      
        % SIR_inter_e = (d/R)^(-gamma) * (1/6) * (sqrt(3*N))^gamma
SIRe = (1/ro) * K * (1/6) * sqrt(3 * N)^gamma;                              % SIR for FFR at the edge equal to the case of full frequency reuse with N reuse factor, d=R in the edge, considers cell load

        % SIR_inter_e = (d/R)^(-gamma) * (1/6) * (sqrt(3))^gamma
SIRc = (1/ro) * K * (1/6) * sqrt(3)^gamma;                                  % SIR for FFR in the central part equal to the case of full frequency reuse with N = 1


% Anonymous functions. To call them use SIRe(some_value for d) / SIRc(some_value for d)
SIRe = @(d) (d/R).^(-gamma). * (1/ro) * k/6 * sqrt(3 * N)^gamma            % SIR for FFR at the cell edge, for hexagonal cell
SIRc = @(d) (d/R).^(-gamma). * (1/ro) * k/6 * sqrt(3)^gamma                % SIR for FFR in the central part, for hexagonal cell

    % SFR: Most frequencies are used in the center (high capacity). Only a part (1/N) of frequencies is used at the edge. If a cell has 12 channels, The center gets 8 channels (used in all cells). The edge gets 4 channels, but each cell uses a different set of 4 channels to reduce interference.
        % Downlink / Uplink:
            % SIR in the central part
                % SIR_inter_c = (d/R)^(-gamma) * (sqrt(3))^(-gamma) / (3 * (1 + omega)))

            % SIR at the edge
                % SIR_inter_e = (d/R)^(-gamma) * (1 / ((6 / omega) * (sqrt(3))^(-gamma) + 6 * (sqrt(3 * N))^(-gamma)))
            % omega is power ratio between the transmission power at the cell edge and the transmission power at the cell center in Soft Frequency Reuse (SFR) schemes:
                % omega = Ptx_e / Ptx_c       (> 1)
SIRe = (1/ro) * (1 / (((6/omega) * sqrt(3)^(- gamma)) + (6 * sqrt(3 * N)^(-gamma))));         % SIR for SFR systems at the cell edge (d=R -> d/R = 1), and considering a cell load

            % Formula for deriving neccassary value of omega to have a SIR at cell edge D = R
                %omega = 6 * (sqrt(3))^(-gamma) / (SIR_inter_e^(-1) - 6 * (sqrt(3 * N))^(-gamma))
omega = (6 * sqrt(3)^(-gamma)) / ((SIRe^-1) - (6 * sqrt(3 * N)^(-gamma)))   % omega - 

% SNR
SNR = S/N                                                                   % Calculate SNR with the ratio of signal power S to noise power N
SNR_db = S_dbm - N_dbm;                                                    
SNR = inv(inv(SINR)-inv(SIR));

SNR_db = @(R) Srx_dbm(R) - (kt0_dbm+10*log10(Bw_eu)+NF_db);                 % Compute SNR at a given distance R
SNR = @(R) 10^(SNR_db(R)/10);

% SINR
    % SINR = P_signal / (P_interferference + P_noise)
        % SIR = P_signal / P_interferference    &   SNR = P_signal / P_noise
    % SINR = 1 / (1/SNR + 1/SIR)
SINR = inv(inv(SNR)+inv(SIR));                                              % Compute SINR based on SNR and SIR
SINRup = inv(inv(SNR) + inv(SIRup));                                        % SINR for the uplink

SINR = @(R) inv(inv(SNR(R))+inv(SIR));                                      % Compute SINR at a given distance 𝑅 R
SINR_db = @(R) 10 * log10(SINR(R));

SINR = 2^(eta_c/alpha) - 1;                                                 % SINR through Shannon Capacity Reu = alpha * Bw * log2(1 + SINR) = eta_c * Bw
SINR = 2^(Reu/(alpha * Bw_eu)) - 1;                                         % SINR through Shannon Capacity Reu = alpha * Bw * log2(1 + SINR)
SINR = @(Reu) 2^(Reu/(alpha * Bw_eu)) - 1;                                  % SINR anonymous function of Shannon Capacity Reu = alpha * Bw * log2(1 + SINR)

% REQUIRED SINR at the cell edge, such that the Actual Average Spectral Efficiency (eta_c) matches the Target Spectral Efficiency (eta_c_req).
    % eta_c(SINR) = eta_c_req
SINR_db = fzero(@(SINR_db) eta_c(SINR_db,gamma) - eta_c_req, 3);            % MWN_ex.7
SINR_db = fzero(@(SINR_db) eta_c(R, SINR_db, gamma, alpha) - eta_c_req, 3);    % for eta_c_req computed as eta_c_req = delta * pi * R^2 / ro / Btot;

SINRe=@(d) (d/R).^(-gamma)./((inv(SNR))+inv(SIRe(R)));                      % The SINR at the cell edge, d = R
SINRc=@(d) (d/R).^(-gamma)./((inv(SNR))+inv(SIRc(R)));                      % The SINR in the central part, d < R




% ---------------------------------------------------------- CAPACITY PLANNING ----------------------------------------------------------
    % Ao [bit/s] - TOTAL required average bit rate from ALL USERS
    % Ao_u [bit/s] - required average bit rate from EACH USER
    % delta [user/km^2] - user density
        % Ao(R) = Ao_u * delta * pi * R^2
C = C_u * delta * pi * R^2;
Ccity = C_u * delta * A_city;                                               % Network capacity of CITY from Capacity planning Ao(R) = Ao_u * delta * pi * R^2
Ccell = C_u * delta * pi * Rcell^2;                                         % Network capacity of CELL from Capacity planning Ao(R) = Ao_u * delta * pi * R^2
Ccell = eta_c * ro * Bw_cell;                                               % Network capacity of CELL 
Ccell = eta_c * ro * Bw_tot / N;                                            % Network capacity of CELL 

delta = N_users/Acity;                                                      % User density - Number of users per km2 (area of the city)
delta = N_users/(pi * R^2);                                                 % User density - Number of users per km2 (area of the cell)
delta = (ro * eta_c * Bcell)/(Ao_u * pi * R^2);                             % User density from Capacity planning (from eta_c = (Ao_u * delta * pi * R^2 / ro) / Bcell)
delta = (ro * eta_c * Bw_tot)/(pi * R^2);                                   % WHERE IS Ao_u??? User density from Capacity planning

    % Cell radio resource should provide Rcell = Ao/ro [bit/s]
    % Required AVERAGE CELL SPECTRAL EFFICIENCY of the cell:
        % eta_c >= (Ao(R) / ro) / Bcell,      Bcell - cell bandwidth
        % eta_c >= (Ao_u * delta * pi * R^2 / ro) / Bcell
eta_c_req = Ao / ro / Bcell;                                                % Required Average Spectral Efficiency (ηc) of a cell
eta_c_req = Ao / ro / Btot;                                                 % Required Average Spectral Efficiency (ηc) of a cell when N = 1
eta_c_req = Ao_u * delta * pi * R^2 / ro / Bcell;                           % Required Average Spectral Efficiency (ηc) of a cell
eta_c_req = delta * pi * R^2 / ro / Btot;                                   % Required Average Spectral Efficiency (ηc) of a cell when N = 1, WHERE IS Ao_u???
eta_c_req = Rcell / Bw_tot / N;                                             % Required Average Spectral Efficiency (ηc) of a cell, Rcell = Ao / ro
eta_c = Ccell / ro / Bw_cell;                                               % Required Average Spectral Efficiency (ηc) of a cell from Capacity planning eta_c_req


% NUMBER OF CELLS
Ncell = ceil(Ccity / Ccell);                                                % Rounding up to cover city capacity

    % AVERAGE CELL SPECTRAL EFFICIENCY can be computed by through the Shannon formula in the point x,y of the cell multiplied for the probability Pr(x,y) of having a mobile located at the x,y point
        % eta_c = integral(integral(alpha * log2(1 + SINR(x,y)) * Pr(x,y) dxdy

    % assuming uniform distribution of the mobiles in the cell, rewrite the formula:
        % eta_c = integral(alpha * log2(1 + SINR(Down, R)) * ((2 * pi * Down) / (pi * R^2)) dDown, 0, R)
            % but in matlab formulas the integral is computed from 1 to R, to avoid numerical instability

    % AVERAGE CELL SPECTRAL EFFICIENCY when SINR about equal to SIR
    % Almost independent of R
        % eta_c = alpha * integral(log2(1 + ((Down / R)^-gamma) * 1/6 * (sqrt(3 * N))^gamma) * ((2 * pi * Down) / (pi * R^2)) dDown, 0, R)
eta = @(D,R,N,gamma,ro) alpha * log2(1+((1/ro) * (1/6) * ((D./R).^-gamma).*(sqrt(3*N)^gamma))).*(2 * pi * D./(pi * R^2));
eta_c = integral(@(D) eta(D,R,N,gamma,ro), 1, R);   % [bps/Hz]

    % AVERAGE CELL SPECTRAL EFFICIENCY from SINR at the edge
    % SINR(R) is the value of the SINR that is guaranteed by link budget at cell edge
    % Almost independent of R
        % eta_c = alpha * integral(log2(1 + ((Down / R)^-gamma) * SINR(R)) * ((2 * pi * Down) / (pi * R^2)) dDown, 0, R)
eta_c = @(Reu) alpha * integral(@(D) log2(1+((D./R).^-gamma).*SINR(Reu)).*(2.*D./R.^2), 1, R);      % SINR depends on Reu, find SINR through shannon capacity formula

        
f = @(D,R,SINR_db,gamma) alpha * log2(1 + (((D./R).^-gamma).*10^(SINR_db/10))).*(2 * pi * D./(pi * R^2));     % Capacity function, defines the formula inside integral, SINR_db converted to linear
eta_c = integral(@(D) f(D,R,SINR_db,gamma), 1, R);                          % R is in meters (by default consider R = 1000 m
eta_c = integral(@(D) f(D,R*1000,SINR_db,gamma), 1, R*1000);                % R is in kilometers

eta_c = @(SINR_db,gamma) integral(@(D) f(D,R,SINR_db,gamma), 1, R);         % anonymous function to call eta_c(SINR_db,gamma)

eta_ce = integral(@(d) fe(d,R),Rc,R);                                       % Spectral efficiency at the edge of the cell
eta_cc=integral(@(d) fc(d,R),1,Rc);                                         % Spectral efficiency at the center of the cell
ec = eta_c(SINR_min_sub_db,gamma,alpha);                                    % average spectral efficiency

% eta_c.m
function res = eta_c(R, SINR_min_db, gamma, alpha)
    % Convert SINR from dB to linear scale
    SINR = 10^(SINR_min_db / 10);

    % Define the function f as per the updated formula
    f = @(D,R,SINR,gamma,alpha) alpha * log2(1 + ((D./R).^(-gamma) .* SINR)) .* (2 * pi * D ./ (pi * R^2));

    % Calculate the integral
    res = integral(@(D) f(D, R, SINR, gamma, alpha), 1, R);
end

% eta_c1.m
function res = eta_c(R,SINR,gamma)
    f = @(D,R,SINR,gamma) 0.75*log2(1+(((D./R).^-gamma).*SINR)).*(2*pi*D./(pi*R^2));
    res = integral(@(D) f(D,R,SINR,gamma), 1,R);
return;

% CELL RADIUS
R = sqrt((eta_c * Bcell * ro) / (Ao_u * delta * pi));                       % Cell Radius from Capacity planning (from eta_c = (Ao_u * delta * pi * R^2 / ro) / Bcell)
R = sqrt(eta_c * Bw_tot * ro / (pi * delta));                               % Cell Radius from Capacity planning, Reuse factor N = 1, (idk where is Ao_u)
R = sqrt(eta_c * Bcell * ro / (delta_t * pi);                               % Cell Radius from Capacity planning, delta_t = Ao_u * delta (delta_t directly represents that total traffic demand per km²) 
R = sqrt(A_city/(pi * Ncell));	                                            % Cell Radius from A_city / Ncell = pi * Rcell^2

% Cell Load from Capacity planning (from eta_c = (Ao_u * delta * pi * R^2 / ro) / Bcell)
ro = (Ao_u * delta * pi * R^2) / ( eta_c * Bcell);
ro = delta * pi * R^2 / (eta_c * Bw_tot);

% Bandwidth from Capacity planning (from eta_c = (Ao_u * delta * pi * R^2 / ro) / Bcell)
Bcell = (Ao_u * delta * pi * R^2) / (eta_c * ro);
Bw_tot=(delta * pi * R^2)/(eta_c * ro)

% AREA
Atot = pi * R^2;                                                            % Total geometric area of a single cell [km2]
delta = N_users/Atot;                                                       % user density
N_users = delta * pi * R^2;                                                 % Total number of users in a area, or active users area Ao, thus considered user density - delta
Ac = pi * Rc^2;                                                             % Total Central area, Rc - radius of central area
Ae = Atot - Ac;                                                             % Total Edge area

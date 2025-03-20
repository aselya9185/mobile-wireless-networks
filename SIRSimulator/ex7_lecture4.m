% 2 example from the lecture 4-RadioPlanning

N_user = 830000;    % Number of users aimed by the operator
Acity = 1285;     % Area of the city 1285km2
% Each userin a month (31 days) consumes on average:
%   a: 2GB of Internet = 2 * 8 * 1e9 bits
%   b: 500 min of calls (12.65 kbps = 12.65 * 1e3 bps) = 500 * 60 * 12.65 * 1e3 bps
%   c: All traffic generated during 12 hours of peak traffic time = 12 * 3600 sec
% During peak time the rate of each user is:
%   (a + b)/ (31 * c) = 12230 bps = 12 kbps:
Ao_u = 12e3;        % Required bit rate for each user
Btot = 15e6;        % Available System Bandwidth
gamma = 4;          % Propagation Exponent (Path Loss)
N = 3;              % Full Reuse Factor
% sec = 0;            No sectoring
ro = 0.5;           % Cell load 50%
alpha = 0.75        % Attenuated Shannon factor
% Assume SINR = SIR
% Determine Capacity Cell Radius Rmax - ?

% SOLUTION:
% 1    
% AVERAGE CELL SPECTRAL EFFICIENCY when SINR about equal to SIR
        % eta_c = alpha * integral(log2(1 + ((Down / R)^-gamma) * 1/6 * (sqrt(3 * N))^gamma) * ((2 * pi * Down) / (pi * R^2)) dDown, 0, R)
eta = @(D,R,N,gamma,ro) alpha * log2(1+((1/ro) * (1/6) * ((D./R).^-gamma).*(sqrt(3*N)^gamma))).*(2 * pi * D./(pi * R^2));
% assuming range R = 1000m for numerical computation of eta_c from 1 to 1000
eta_c = integral(@(D) eta(D,1000,N,gamma,ro), 1, 1000);

% 2
% User density:
delta = N_user / Acity;     % N_users / km2
Bcell = Btot / N;
% Total average bit rate:
    % Ao(R) = Ao_u * delta * pi * R^2 = 12e3 * delta * pi * R^2
    % eta_c = Ao(R) / ro / Bcell ---> eta_c = (Ao_u * delta * pi * R^2 / ro) / Bcell  --->
    % R = sqrt((eta_c * Bcell * ro) / (Aou * delta * pi));      
% Cell Radius in km
R = sqrt((eta_c * Bcell * ro) / (Ao_u * delta * pi));       % Output: 0.7679 km
% Compute and plot the average spectral efficiency of a cell versus the 
% bit-rate Reu provided to a user at the cell edge, in case of 
% Reu = 250 kbit/s, 500 kbit/s, 1 Mbit/s:
% - Capacity evaluation through Shannon formula with attenuation factor alpha =0.75
% - Radio bandwidth provided to the edge user Beu = 180 kHz
% - Path Loss exponent = 4
% - Cell radius R = 1 km

% nc/Reu ?
clear all;
Reu_v = [250e3 500e3 1e6]
gamma = 4;
alpha = 0.75;
Bw_eu = 180e3;
R = 1000;

res = zeros(size(Reu_v));

% Reu = alpha * Beu * log2(1+SINR))
SINR = @(Reu) 2.^(Reu./(alpha * Bw_eu)) - 1;                                  % SINR anonymous function of Shannon Capacity Reu = alpha * Bw * log2(1 + SINR)

f = @(D,R,SINR,gamma) alpha * log2(1 + (((D./R).^-gamma).*SINR)).*(2 * pi * D./(pi * R^2));     % Capacity function, defines the formula inside integral, SINR_db converted to linear
eta_c = @(Reu) integral(@(D) f(D,R,SINR(Reu),gamma), 1, R);

% eta_c = @(Reu) alpha*integral(@(D)log2(1+((D./R).^-gamma).*SINR(Reu)).*(2.*D./R.^2),1,R)

for i = 1:length(Reu_v)
    res(i) = eta_c(Reu_v(i));
end

fprintf('spectral efficiency: %.2f \n', res);
% Plot bit-rate vs. spectral efficiency (nc)
figure;
plot(Reu_v, res, '-o');
xlabel('Bit-rate at Edge User (bps)');
ylabel('Spectral efficiency (bps/Hz)');
title('Bit-rate vs. Spectral efficiency (nc)');
grid on;
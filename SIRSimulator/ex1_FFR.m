% Fractional Frequency Reuse (FFR)
% 
% Find:
% Edge Spectral Efficiency (η_ce)
% Central Spectral Efficiency (η_cc)
% Average Cell Capacity (Ccell)


clear all

gamma = 4; 
ro = 0.7; 
N = 3; 
sec = 3; 
R = 1000; 
SNR_db = 13;    % SNR at the cell edge
alpha = 0.75; 
Btot = 5e6;     % Hz
Bc = Btot/2;    % bandwidth used at the cell center
Be = Btot/6;    %  bandwidth used by a single cell at the edge
% Switching point (edge/central) at a distance Rc where the SIRc(Rc) == SIRe(Rc)

% SIR formulas for edge and central parts
SIRe= @(d) (d/R).^(-gamma).*(1/ro)*sec/6*sqrt(3*N)^gamma;
SIRc= @(d) (d/R).^(-gamma).*(1/ro)*sec/6*sqrt(3)^gamma;

% Find switching point Rc where SIRc(Rc) == SIRe(R)
Rc=fsolve(@(d) SIRc(d)-SIRe(R),700);

% Convert SNR from dB to linear scale
SNR = 10^(SNR_db/10);

% SINR calculations
SINRe=@(d) (d/R).^(-gamma)./(inv(SNR) + inv(SIRe(R)));
SINRc=@(d) (d/R).^(-gamma)./(inv(SNR) + inv(SIRc(R)));

% Total, central, and edge areas
Atot = pi*R^2;
Ac = pi*Rc^2;
Ae = Atot - Ac;

% Spectral efficiency at the edge (η_ce) and central part (η_cc)
fe = @(d,R) alpha*log2(1 + SINRe(d)) .* (2*pi*d./Ae);
eta_ce = integral(@(d) fe(d,R), Rc, R);

fc = @(d,R) alpha*log2(1 + SINRc(d)) .* (2*pi*d./Ac);
eta_cc = integral(@(d) fc(d,R), 1, Rc);

% Compute average cell capacity (Ccell)
Ccell = (Ae/Atot)*(Be)*eta_ce + (Ac/Atot)*(Bc)*eta_cc;

% Print results
fprintf('Edge Spectral Efficiency (η_ce): %.4f bps/Hz\n', eta_ce);
fprintf('Central Spectral Efficiency (η_cc): %.4f bps/Hz\n', eta_cc);
fprintf('Average Cell Capacity (Ccell): %.4f Mbps\n', Ccell/1e6);

% Output:
% Edge Spectral Efficiency (η_ce): 3.9280 bps/Hz
% Central Spectral Efficiency (η_cc): 6.2607 bps/Hz
% Average Cell Capacity (Ccell): 7.3995 Mbps
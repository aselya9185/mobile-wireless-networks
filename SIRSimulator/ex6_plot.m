R = 2000;
Ptx_user_dbm = 40;
N_v = [1 3 7];
cap=[];
Reu=[];
% kt0_dbm=-174;

No_dbm = kt0_dbm + NF_db + 10*log10(Bw_eu);
L_db = L0_db + 10*gamma*log10(R/R0);
Prx_user_dbm = Ptx_user_dbm - L_db;
M_db = fzero(@(x) F-0.5*erfc(-x/(sigma_db*sqrt(2))),10);
S_dbm = Prx_user_dbm - M_db;

for i=1:length(N_v)
    N = N_v(i); % reuse factor
    SIR = 1/6 * 1/ro * s * sqrt(3*N)^gamma;
    SNR_db = S_dbm - No_dbm;
    SNR = 10^(SNR_db/10);
    SINR = inv(inv(SNR)+inv(SIR));
    
    Reu(i) = alpha * Bw_eu * log2(1+SINR);

    f = @(D,R,SINR,gamma) 0.75*log2(1+(((D./R).^-gamma).*SINR)).*(2*pi*D./(pi*R^2));
    eta_c = integral(@(D) f(D,R,SINR,gamma), 1,R);
    cap(i) =  eta_c * Bw_tot/N;

end

figure
plot(N_v,Reu/1000)
xlabel('Reuse Factor \itN')
ylabel('Edge user bit rate (kbit/s')

figure
plot(N_v,cap/1000)
xlabel('Reuse Factor \itN')
ylabel('Average cell capacity (kbit/s')


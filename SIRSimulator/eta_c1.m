function res = eta_c(R,SINR,gamma)
    f = @(D,R,SINR,gamma) 0.75*log2(1+(((D./R).^-gamma).*SINR)).*(2*pi*D./(pi*R^2));
    res = integral(@(D) f(D,R,SINR,gamma), 1,R);
return;
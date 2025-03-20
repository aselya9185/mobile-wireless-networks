function res = eta_c(R, SINR_min_db, gamma, alpha)
    % Convert SINR from dB to linear scale
    SINR = 10^(SINR_min_db / 10);

    % Define the function f as per the updated formula
    f = @(D,R,SINR,gamma,alpha) alpha * log2(1 + ((D./R).^(-gamma) .* SINR)) .* (2 * pi * D ./ (pi * R^2));

    % Calculate the integral
    res = integral(@(D) f(D, R, SINR, gamma, alpha), 1, R);
end


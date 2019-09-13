function [E_D, lnLambda, collTime] ...
                                    = physical_parameters(n, Z0, T)
    e = 1.60217662e-19; % C
    eps0 = 8.8541878176e-12; % F/m
    me = 9.10938356e-31; % kg
    
    ne = dot(n, Z0);
    lnLambda = 14.9 - 0.5 * log(ne/1e20) + log(T/1e3);
    E_D = ne * e^3 * lnLambda / (4 * pi * eps0^2 * e * T);
    collTime = 8 * sqrt(2) * pi * (e * T) ^ (3/2) * eps0 ^ 2 * sqrt(me) / (lnLambda * e ^ 4 * ne);
end
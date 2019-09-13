function G = primary(n, Z, Z0, EoED, T)
% n, Z, Z0 are column vectors with same length as the number of species.
% n contains the densities (m^-3).
% Z contains the atomic numbers.
% Z0 contains the charges.
% EoED is the electric field normalized to the Dreicer field.
% T is the temperature in eV.
% When the inputs are matrices, the function will be applied to each column
% seperately.
% G = t_ee / n_e * (d/dt * n_re), where t_ee is the thermal electron
% collision time, n_e is the electron density, and n_re is the runaway
% density.
% t_ee = (4 * pi * eps0 ^ 2 * m_e ^ 2 * v_Te ^ 3) / (n_e * e ^ 4 * lnLambda)
% Where eps0 is the permittivity of vacuum, m_e is the electron mass, and
% v_Te = sqrt(2 * T / m_e) is the thermal speed.
% Calculation of t_ee, ED, and the coulomb logarithm can be found in
% physical_parameters.m.

nfree = dot(n, Z0);
ntot = dot(n, Z);

Zeff = dot(n, Z0 .^ 2) ./ nfree;
Zeff0 = dot(n, Z .^ 2 - Z0 .^ 2) ./ ntot;
Z0oZ = dot(n, Z0 ./ Z) ./ ntot;
ZZ0 = dot(n, Z0 .* Z) ./ ntot;
lnnfree = log(nfree);
nfreeontot = nfree ./ ntot;

G = primary_derived(EoED, T, Zeff, Zeff0, ZZ0, Z0oZ, lnnfree, nfreeontot);

G = G * 4/(3 * sqrt(pi)); % Normalize units.
end
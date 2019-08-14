function G = primary(n, Z, Z0, EoED, T)
% n, Z, Z0 are column vectors with same length as the number of species.
% n contains the densities (m^-3).
% Z contains the atomic numbers.
% Z0 contains the charges.
% EoED is the electric field normalized to the Dreicer field.
% T is the temperature in ev.
% When the inputs are matrices, the function will be applied to each column
% seperately.

nfree = dot(n, Z0);
ntot = dot(n, Z);

Zeff = dot(n, Z0 .^ 2) ./ nfree;
Zeff0 = dot(n, Z .^ 2 - Z0 .^ 2) ./ ntot;
Z0oZ = dot(n, Z0 ./ Z) ./ ntot;
ZZ0 = dot(n, Z0 .* Z) ./ ntot;
lnnfree = log(nfree);
nfreeontot = nfree ./ ntot;

G = primary_derived(EoED, T, Zeff, Zeff0, ZZ0, Z0oZ, lnnfree, nfreeontot);

end
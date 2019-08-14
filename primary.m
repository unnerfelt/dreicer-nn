function G = primary(n, Z, Z0, EoED, T)

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
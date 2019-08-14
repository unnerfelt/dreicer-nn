function G = primary_derived(EoED, T, Zeff, Zeff0, ZZ0, Z0oZ, lnnfree, nfreeontot)
% T is in eV.
  input = [Zeff; Zeff0; Z0oZ; ZZ0; lnnfree; nfreeontot; EoED; log(T ...
                                                    / 510998.946)];
  G = neural_network(input);
end
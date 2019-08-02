function G = primary(EoED, T, Zeff, Zeff0, ZZ0, Z0oZ, nfree, ntot)
% T is in eV.
  input = [Zeff; Zeff0; ZZ0; Z0oZ; log(nfree); nfree / ntot; EoED; log(T ...
                                                    / 510998.946)];
  G = neural_network(input);
end
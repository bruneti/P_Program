%%
%% This function computes the Fermi energy (EF) and the lower band edge
%% (BandEdge) from the eigen-energies enk, by bisection of the total charge.
%%
%% Written by J. Ferrer, October 2020
%%

function [ EF, BandEdge ] = util_FermiEnergy( enk, post, siesta, constants )

  Nk      = post.kpoints.Nk;
  wk      = post.kpoints.reciprocal_vectors(:,4);

  T      = siesta.T*constants.kB/constants.eV;
  Nspin  = siesta.spin.Nspin;
  Qtol   = 1.0e-5;
  Q      = siesta.Nelectrons;

  if Nspin == 1, Q  = Q/2; end

  eigen    = sort( reshape(enk,[],1), 'ascend' );
  BandEdge = min(eigen);
  %% Bracket EF
  if Q*Nk-1==0
      Emin=eigen(1);
  else 
      Emin     = eigen(floor(Q*Nk-1));
  end
  Qmin     = sum( (1./(exp((enk-Emin)./T)+1)) * wk, 'all' );
  ie = 0 ;
  while Qmin >= Q
    ie = ie - 1;
    Emin = eigen(floor(Q*Nk)+ie);
    Qmin = sum( (1./(exp((enk-Emin)./T)+1)) * wk, 'all' )
  end
  if Q*Nk==1
      Emax=eigen(2);
  else
  Emax = eigen(floor(Q*Nk)+2);
  end
  Qmax = sum( (1./(exp((enk-Emax)./T)+1)) * wk, 'all' );
  ie = 1;
  while Qmax <= Q
    ie = ie + 1;
    Emax = eigen(floor(Q*Nk)+ie);
    Qmax = sum( (1./(exp((enk-Emax)./T)+1)) * wk, 'all' );
  end

  %% Find EF
  Q0 = Qmax;
  while abs(Q0-Q) > Qtol
    EF0 = (Emax+Emin)/2;
    Q0 = sum( (1./(exp((enk-EF0)./T)+1)) * wk, 'all' );
    if Q0 > Q 
      Emax = EF0;
    else
      Emin = EF0;
    end
  end

  EF = EF0;

end

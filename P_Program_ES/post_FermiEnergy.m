%%
%% Calcula la energía de Fermi: diagonaliza en todos los k-points y llama a
%% util_FermiEnergy con las autoenergías resultantes.
%%
%% Escrito por J. Ferrer, octubre 2020
%%

function [ post ] = post_FermiEnergy( post, siesta, constants, flags )

  Nk      = post.kpoints.Nk;
  k       = 2.0*pi*post.kpoints.reciprocal_vectors(:,1:3);

  R      = siesta.lattice.supercell.R;
  nou    = siesta.Norbs.ucell;
  Nspin  = siesta.spin.Nspin;

  eikr   = exp( 1i* k*R );

  %% Calcula todas las autoenergías
  enk = zeros(nou,Nk);
  if Nspin == 2, enk = zeros(2*nou,Nk); end;
  for ik = 1:Nk
    if strcmp( flags.full_OHDM, 'yes')
      [ Hk, Ok ] = util_Hk( k(ik,:), siesta );
    else
      [ Hk, Ok ] = util_Hk_sparse( k(ik,:), siesta );  
    end
    if Nspin == 2
      enk0 = [];
      for ispin = 1:Nspin
        enk0 = [ enk0; eig(Hk(:,:,ispin),Ok,'vector') ];   
      end
      enk(:,ik) = enk0;
    else
      enk(:,ik) = eig(Hk,Ok,'vector');
    end
  end
  [ post.EF, post.BandEdge ] = util_FermiEnergy( enk, post, siesta, constants ); 

end

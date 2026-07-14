%%
%% Calcula las autoenergías de Kohn-Sham para un k-point dado.
%%
%% Escrito por J. Ferrer, octubre 2020
%

function [ post ] = post_eigenenergy( post, siesta, flags )

  k(1:3)  = 2.0*pi*siesta.kpoints.reciprocal_vectors(flags.kpoint.ik,1:3);
  post.eigen_k=k;
  if strcmp( flags.full_OHDM, 'yes')
    [ Hk, Ok ] = util_Hk( k, siesta );
  else
    [ Hk, Ok ] = util_Hk_sparse( k, siesta );
  end

  if siesta.spin.Nspin == 2
  post.eigen_energies = [ eig(Hk(:,:,1),Ok); eig(Hk(:,:,2),Ok)] ;
  [post.eigen_vectors1,~]=eig(Hk(:,:,1),Ok);
  [post.eigen_vectors2,~]=eig(Hk(:,:,2),Ok);
  post.eigen_vectors = [post.eigen_vectors1; post.eigen_vectors2];
  else
   [post.eigen_vectors, post.eigen_energies] = eig(Hk,Ok); 
  end
      
end
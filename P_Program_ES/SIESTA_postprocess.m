%% SIESTA_postprocess.m
%
% Post-procesado mínimo necesario para el pipeline de tasa de transición
% dipolar. Produce:
%
%     post.kpoints        - K-points de la simulación original.
%     post.EF             - Energía de Fermi.
%     post.BandEdge       - Borde inferior del espectro.
%     post.eigen_energies - Autoenergías al k-point flags.kpoint.ik
%                           (necesarias para identificar HOMO/LUMO en
%                           Calculo_Final).
%     post.eigen_vectors  - Autovectores asociados (usados en Calculo_Final
%                           cuando Nspin = 2; con Nspin = 1 se recalculan
%                           dentro del propio Calculo_Final).
%     post.bands          - (opcional) Bandas a lo largo de una línea de
%                           k-points.
%     post.k              - (opcional) Lista de k-points consumida por
%                           Reordenar cuando se trabaja con más de un
%                           k-point.
%
% Extraído de GROGU (J. Ferrer, 2020). Se han eliminado las ramas de
% exchange, DM, DOS/PDOS, fat-bands, Mulliken, populations,
% split-Hamiltonian e invariante Z2.

function post = SIESTA_postprocess( siesta, constants, flags )

  post.init    = 'yes';
  post.kpoints = siesta.kpoints;

%% ENERGÍA DE FERMI
  disp('Calculando energía de Fermi');
  tic
  post = post_FermiEnergy( post, siesta, constants, flags );
  toc

%% AUTOENERGÍAS PARA UN K-POINT (HOMO/LUMO en Calculo_Final)
  if strcmp(flags.post.eigenenergy,'yes')
    post = post_eigenenergy( post, siesta, flags );
  end

%% BANDAS A LO LARGO DE UNA LÍNEA
%  Sólo se calculan si la supercelda tiene más de una celda. Generan
%  'post.k', necesaria en Reordenar.m cuando hay más de un k-point.
  if strcmp(flags.post.bands,'yes') && siesta.lattice.supercell.Ncells > 1
    disp(' Calculando bandas');
    tic
    post = post_bands( post, siesta, flags );
    toc
  end

  clear post.init

end

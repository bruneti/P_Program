%% SIESTA_postprocess.m
%
% Minimal post-processing needed by the dipole transition-rate pipeline.
% Produces:
%
%     post.kpoints        - K-points of the original simulation.
%     post.EF             - Fermi energy.
%     post.BandEdge       - Lower edge of the spectrum.
%     post.eigen_energies - Eigen-energies at k-point flags.kpoint.ik
%                           (needed to identify HOMO/LUMO in Calculo_Final).
%     post.eigen_vectors  - Associated eigenvectors (used in Calculo_Final
%                           when Nspin = 2; with Nspin = 1 they are
%                           recomputed inside Calculo_Final itself).
%     post.bands          - (optional) Bands along a k-point line.
%     post.k              - (optional) List of k-points consumed by
%                           Reordenar when working with more than one
%                           k-point.
%
% Extracted from GROGU (J. Ferrer, 2020). The exchange, DM, DOS/PDOS,
% fat-bands, Mulliken, populations, split-Hamiltonian and Z2 branches have
% been removed.

function post = SIESTA_postprocess( siesta, constants, flags )

  post.init    = 'yes';
  post.kpoints = siesta.kpoints;

%% FERMI ENERGY
  disp('Computing Fermi energy');
  tic
  post = post_FermiEnergy( post, siesta, constants, flags );
  toc

%% EIGEN-ENERGIES AT ONE K-POINT (HOMO/LUMO in Calculo_Final)
  if strcmp(flags.post.eigenenergy,'yes')
    post = post_eigenenergy( post, siesta, flags );
  end

%% BANDS ALONG A LINE
%  Only computed if the super-cell has more than one cell. They produce
%  'post.k', needed by Reordenar.m when there is more than one k-point.
  if strcmp(flags.post.bands,'yes') && siesta.lattice.supercell.Ncells > 1
    disp(' Computing bands');
    tic
    post = post_bands( post, siesta, flags );
    toc
  end

  clear post.init

end

%% Cargar_SIESTA.m
%
% Loads the relevant data of a SIESTA simulation into MATLAB memory and
% leaves it ready to feed the radiative dipole transition-rate calculation.
%
% Full pipeline:
%     Cargar_SIESTA  ->  Pain_Program  ->  Reordenar  ->  Calculo_Final
%
% Structures left in the workspace when it finishes:
%
%     siesta     -  Geometry, species, orbitals (with their quantum
%                   numbers), PAOs, Hamiltonian and overlap in super-cell
%                   form.
%     post       -  Eigen-energies at the selected k-point (needed to
%                   identify HOMO and LUMO in Calculo_Final), Fermi energy
%                   and, optionally, bands along a k-point line.
%     constants  -  Physical constants (SI) used throughout the pipeline.
%
% Extracted and reduced from the original GROGU package
% (J. Ferrer & L. Oroszlany, 2021-2022). This version keeps only what is
% strictly needed for the dipole transition rate and removes the branches
% for exchange, DOS/PDOS, fat-bands, Mulliken, populations,
% split-Hamiltonian, Z2 invariant and Xatu export.

%% INITIALIZE
  clear all;
  runstart = num2cell(clock);
  disp( [ ' Run started  on ', datestr(datenum(runstart{:})) ] );
  tic

%% PHYSICAL CONSTANTS
  constants = SIESTA_constants;

%% FLAGS AND INPUT PATHS
  flags = SIESTA_flags;

%% READ THE SIESTA FILES
  disp('Reading SIESTA')
  tic
  siesta = SIESTA_reader( constants, flags );
  toc

%% MINIMAL POST-PROCESSING (Fermi, eigen-energies at the k-point, optional bands)
  post = SIESTA_postprocess( siesta, constants, flags );

%% SAVE THE WORKSPACE
%   siesta.mat   - compatible with Pain_Program, which expects the variable
%                  'siesta' at run time.
%   post.mat     - holds 'post' and 'constants' for Calculo_Final.
  save('siesta.mat','siesta')
  save('post.mat','post','constants')

%% RUN TIME
  runstop = num2cell(clock);
  pp      = toc;
  disp( [ ' Run finished on ', datestr(datenum(runstop{:})) ] );
  disp( [ ' Loading took ', num2str(pp), ' seconds' ]);
  clear pp directory runstart runstop

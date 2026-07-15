%% SIESTA_flags.m
%
% Minimal configuration needed to load a SIESTA simulation in the context of
% the dipole transition-rate calculation.
%
% Only the flags actually consumed by the following are kept:
%   - the SIESTA readers (geometry, H/O/DM, k-points, .EIG, ...)
%   - the Fermi-energy calculation
%   - the eigen-energies at the selected k-point
%   - the (optional) band calculation that produces post.k when working with
%     more than one k-point
%
% Extracted from GROGU (J. Ferrer, 2020).

function flags = SIESTA_flags

%% INPUT PATHS   >>> USER CHOICE: point these to your SIESTA simulation
  flags.input.directory = ['../Ejemplo_1/'];
  flags.input.file      = 'hBN.out';

%% READ FORMAT OF O, H, DM:  NetCDF ('yes') vs OLD BINARY ('no')
  flags.netcdf          = 'yes';

%% O, H, DM IN DENSE FORM  (not advisable for very large systems)
  flags.full_OHDM       = 'yes';

%% BINARY-FORMAT PRECISION (only used if flags.netcdf = 'no')
%  SIESTA writes DM in single precision.
%  For H, S: single (no spin-orbit) or double (with spin-orbit).
  if strcmp(flags.netcdf, 'no')
    flags.precision.HO  = 'single';
    flags.precision.DM  = 'double';
  end

%% TOLERANCE TO TRUNCATE SMALL ELEMENTS OF O, H, DM
  flags.tol             = 1.e-7;

%% EIGEN-ENERGIES AT ONE K-POINT  (needed to identify HOMO/LUMO)
  flags.post.eigenenergy = 'yes';
  flags.kpoint.ik        = 1;   % >>> USER CHOICE: k-point index for the eigen-energies

%% BANDS ALONG A K-POINT LINE
%  Produces the list 'post.k' that Reordenar.m consumes when there is more
%  than one k-point. For a Gamma-only system it can be left as 'no' without
%  affecting the dipole calculation.
  flags.post.bands          = 'yes';
  flags.bands.kline.vectors = [ [ 0 0 0 ]; [1 2 0]./3; [1,1,0]./2; [0,0,0] ];  % >>> USER CHOICE: k-line
  flags.bands.kline.Nk      = [ 200; 100; 200 ];                                % >>> USER CHOICE: points per segment

end

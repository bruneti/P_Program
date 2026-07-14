%% SIESTA_reader.m
%
% Orchestrates the reading of the SIESTA files and fills the 'siesta'
% structure with everything the dipole pipeline needs.
%
% Sub-structures produced:
%   siesta.lattice              -  Lattice constant, direct and reciprocal
%                                  vectors, super-cell.
%   siesta.H, .O, .DM           -  Rectangular (sparse) form.
%   siesta.Hsc, .Osc, .DMsc     -  Dense super-cell form.
%   siesta.spin                 -  Nspin and description of the calculation.
%   siesta.species              -  Species and PAO information
%                                  (Rlz, ordering, grid, ...).
%   siesta.atom                 -  Unit-cell and super-cell atoms with their
%                                  coordinates.
%   siesta.orb                  -  Each orbital with (atom, l, m, species,
%                                  zeta, ...).
%   siesta.Norbs.{ucell,scell}  -  Number of orbitals.
%   siesta.kpoints              -  K-points used by SIESTA.
%   siesta.EF, .eigen_energies  -  (legacy) taken from the .EIG. The updated
%                                  value lives in 'post'.
%
% Extracted from GROGU (J. Ferrer, 2020).

function siesta = SIESTA_reader( constants, flags )

%% ADD THE SIESTA DATA DIRECTORY TO THE PATH
  addpath( flags.input.directory );

%% READ THE MAIN SIESTA OUTPUT FILE
  logfile        = fileread( flags.input.file );
  siesta.logfile = logfile;

%% SIESTA VERSION
  foundat        = strfind( logfile, 'Siesta version' );
  dummy          = splitlines( logfile( foundat:foundat+100));
  dummy          = strsplit( dummy{1},{' ', ':','siesta-'});
  siesta.version = dummy(3);

%% SystemLabel
  foundat             = strfind( logfile, 'SystemLabel' );
  dummy               = strsplit( logfile( foundat:foundat+80 ) );
  label               = dummy(2);
  siesta.system_label = label{1};

%% ELECTRONIC TEMPERATURE
  foundat      = strfind( logfile, 'ElectronicTemperature');
  dummy        = strsplit( logfile( foundat:foundat+40 ) );
  siesta.T     = str2double(dummy(2));

%% NUMBER OF ELECTRONS
  foundat           = strfind( logfile, 'electrons:');
  dummy             = strsplit( logfile( foundat:foundat+40 ) );
  siesta.Nelectrons = str2double(dummy(2));

%% GEOMETRY AND SUPER-CELL
  siesta = read_lattice( siesta, constants, flags );

%% H, O, DM AND SPIN STRUCTURE
  siesta = read_hodm( siesta, constants, flags );

%% SPECIES AND PAOs
  siesta = read_species( logfile, siesta, flags );

%% UNIT-CELL AND SUPER-CELL ATOMS
  siesta = read_atom( siesta, constants, flags );

%% ORBITALS (indices, quantum numbers, species, zeta, ...)
  siesta = read_orb( siesta, flags );

%% K-POINTS OF THE ORIGINAL SIESTA SIMULATION
  siesta = read_kpoints( siesta, constants, flags );

%% EIGEN-ENERGIES READ FROM THE .EIG (initial reference; post recomputes them)
  siesta = read_eig( flags, siesta );

%% CLEAN-UP
  siesta = rmfield(siesta,'logfile');

end

%% SIESTA_reader.m
%
% Orquesta la lectura de los ficheros de SIESTA y rellena la estructura
% 'siesta' con todo lo necesario para el pipeline dipolar.
%
% Sub-estructuras producidas:
%   siesta.lattice              -  Constante de red, vectores directos y
%                                  recíprocos, supercelda.
%   siesta.H, .O, .DM           -  En forma rectangular (sparse).
%   siesta.Hsc, .Osc, .DMsc     -  En forma supercelda densa.
%   siesta.spin                 -  Nspin y descripción del tipo de cálculo.
%   siesta.species              -  Información de especies y PAOs
%                                  (Rlz, ordering, grid, ...).
%   siesta.atom                 -  Átomos de unit cell y supercelda con
%                                  sus coordenadas.
%   siesta.orb                  -  Cada orbital con (atom, l, m, species,
%                                  zeta, ...).
%   siesta.Norbs.{ucell,scell}  -  Número de orbitales.
%   siesta.kpoints              -  K-points usados por SIESTA.
%   siesta.EF, .eigen_energies  -  (legacy) tomados del .EIG. El valor
%                                  actualizado vive en 'post'.
%
% Extraído de GROGU (J. Ferrer, 2020).

function siesta = SIESTA_reader( constants, flags )

%% AÑADIR EL DIRECTORIO DE DATOS DE SIESTA AL PATH
  addpath( flags.input.directory );

%% LEER FICHERO PRINCIPAL DE SALIDA DE SIESTA
  logfile        = fileread( flags.input.file );
  siesta.logfile = logfile;

%% VERSIÓN DE SIESTA
  foundat        = strfind( logfile, 'Siesta version' );
  dummy          = splitlines( logfile( foundat:foundat+100));
  dummy          = strsplit( dummy{1},{' ', ':','siesta-'});
  siesta.version = dummy(3);

%% SystemLabel
  foundat             = strfind( logfile, 'SystemLabel' );
  dummy               = strsplit( logfile( foundat:foundat+80 ) );
  label               = dummy(2);
  siesta.system_label = label{1};

%% TEMPERATURA ELECTRÓNICA
  foundat      = strfind( logfile, 'ElectronicTemperature');
  dummy        = strsplit( logfile( foundat:foundat+40 ) );
  siesta.T     = str2double(dummy(2));

%% NÚMERO DE ELECTRONES
  foundat           = strfind( logfile, 'electrons:');
  dummy             = strsplit( logfile( foundat:foundat+40 ) );
  siesta.Nelectrons = str2double(dummy(2));

%% GEOMETRÍA Y SUPERCELDA
  siesta = read_lattice( siesta, constants, flags );

%% H, O, DM Y ESTRUCTURA DE ESPÍN
  siesta = read_hodm( siesta, constants, flags );

%% ESPECIES Y PAOs
  siesta = read_species( logfile, siesta, flags );

%% ÁTOMOS DE LA UNIT CELL Y LA SUPERCELDA
  siesta = read_atom( siesta, constants, flags );

%% ORBITALES (índices, números cuánticos, especies, zeta, ...)
  siesta = read_orb( siesta, flags );

%% K-POINTS DE LA SIMULACIÓN ORIGINAL DE SIESTA
  siesta = read_kpoints( siesta, constants, flags );

%% AUTOENERGÍAS LEÍDAS DEL .EIG (referencia inicial; post las recalcula)
  siesta = read_eig( flags, siesta );

%% LIMPIEZA
  siesta = rmfield(siesta,'logfile');

end

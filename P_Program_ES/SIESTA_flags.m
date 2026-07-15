%% SIESTA_flags.m
%
% Configuración mínima necesaria para cargar una simulación de SIESTA en
% el contexto del cálculo de la tasa de transición dipolar.
%
% Sólo se conservan los flags efectivamente consumidos por:
%   - los lectores de SIESTA (geometría, H/O/DM, k-points, .EIG, ...)
%   - el cálculo de la energía de Fermi
%   - el cálculo de autoenergías al k-point seleccionado
%   - el cálculo (opcional) de bandas para obtener post.k cuando se trabaja
%     con más de un k-point
%
% Extraído de GROGU (J. Ferrer, 2020).

function flags = SIESTA_flags

%% RUTAS DE ENTRADA   >>> ELECCIÓN DEL USUARIO: apuntar a tu simulación de SIESTA
  flags.input.directory = ['../Ejemplo_1/'];
  flags.input.file      = 'hBN.out';

%% FORMATO DE LECTURA DE O, H, DM:  NetCDF ('yes') vs BINARIO ('no')
  flags.netcdf          = 'yes';

%% MATRICES O, H, DM EN FORMA DENSA  (no recomendable para sistemas enormes)
  flags.full_OHDM       = 'yes';

%% PRECISIÓN DEL FORMATO BINARIO (sólo se utiliza si flags.netcdf = 'no')
%  SIESTA escribe DM en single precision.
%  Para H, S: single (sin spin-orbit) o double (con spin-orbit).
  if strcmp(flags.netcdf, 'no')
    flags.precision.HO  = 'single';
    flags.precision.DM  = 'double';
  end

%% TOLERANCIA PARA TRUNCAR ELEMENTOS PEQUEÑOS DE O, H, DM
  flags.tol             = 1.e-7;

%% AUTOENERGÍAS A UN K-POINT  (necesario para identificar HOMO/LUMO)
  flags.post.eigenenergy = 'yes';
  flags.kpoint.ik        = 1;   % >>> ELECCIÓN DEL USUARIO: índice del k-point para las autoenergías

%% BANDAS A LO LARGO DE UNA LÍNEA DE K-POINTS
%  Genera la lista 'post.k' que Reordenar.m consume cuando hay más de un
%  k-point en la simulación. Si el sistema es sólo en Gamma se puede dejar
%  en 'no' sin afectar al cálculo dipolar.
  flags.post.bands          = 'yes';
  flags.bands.kline.vectors = [ [ 0 0 0 ]; [1 2 0]./3; [1,1,0]./2; [0,0,0] ];  % >>> ELECCIÓN DEL USUARIO: línea de k
  flags.bands.kline.Nk      = [ 200; 100; 200 ];                                % >>> ELECCIÓN DEL USUARIO: puntos por segmento

end

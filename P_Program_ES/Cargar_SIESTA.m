%% Cargar_SIESTA.m
%
% Carga los datos relevantes de una simulación de SIESTA en la memoria de
% MATLAB y los deja listos para alimentar el cálculo de la tasa de
% transición dipolar radiativa.
%
% Pipeline completo:
%     Cargar_SIESTA  ->  Pain_Program  ->  Reordenar  ->  Calculo_Final
%
% Estructuras que quedan en el workspace al terminar:
%
%     siesta     -  Geometría, especies, orbitales (con sus números
%                   cuánticos), PAOs, Hamiltoniano y solapamiento en forma
%                   supercelda.
%     post       -  Autoenergías al k-point seleccionado (necesarias para
%                   identificar HOMO y LUMO en Calculo_Final), energía de
%                   Fermi y, opcionalmente, bandas a lo largo de una línea
%                   de k-points.
%     constants  -  Constantes físicas en SI utilizadas a lo largo del
%                   pipeline.
%
% Extraído y reducido del paquete original GROGU
% (J. Ferrer & L. Oroszlany, 2021-2022). Esta versión conserva sólo lo
% estrictamente necesario para el cálculo de la tasa de transición dipolar
% y elimina las ramas relativas a exchange, DOS/PDOS, fat-bands, Mulliken,
% populations, split-Hamiltonian, invariante Z2 y exportación a Xatu.

%% INITIALIZE
  clear all;
  runstart = num2cell(clock);
  disp( [ ' Ejecución iniciada el ', datestr(datenum(runstart{:})) ] );
  tic

%% CONSTANTES FÍSICAS
  constants = SIESTA_constants;

%% FLAGS Y RUTAS DE ENTRADA
  flags = SIESTA_flags;

%% LECTURA DE LOS FICHEROS DE SIESTA
  disp('Leyendo SIESTA')
  tic
  siesta = SIESTA_reader( constants, flags );
  toc

%% POST-PROCESADO MÍNIMO (Fermi, autoenergías al k-point, bandas opcionales)
  post = SIESTA_postprocess( siesta, constants, flags );

%% GUARDADO DEL WORKSPACE
%   siesta.mat   - compatible con Pain_Program, que espera la variable
%                  'siesta' al ejecutarse.
%   post.mat     - contiene 'post' y 'constants' para Calculo_Final.
  save('siesta.mat','siesta')
  save('post.mat','post','constants')

%% TIEMPO DE EJECUCIÓN
  runstop = num2cell(clock);
  pp      = toc;
  disp( [ ' Ejecución terminada el ', datestr(datenum(runstop{:})) ] );
  disp( [ ' La carga ha durado ', num2str(pp), ' segundos' ]);
  clear pp directory runstart runstop

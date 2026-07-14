%%
%% Lee los datos estructurales de SIESTA.
%%
%% Rellena la sub-estructura 'lattice' (constante de red, vectores directos
%% y recíprocos, información de supercelda).
%%
%% Escrito por J. Ferrer, 2020
%%

function [ siesta ] = read_lattice( siesta, constants, flags );

logfile = siesta.logfile;
fstr    = [ siesta.system_label '.STRUCT_OUT' ];

%% Constante de red
foundat = strfind(logfile,'LatticeConstant');
dummy   = strsplit(logfile(foundat:foundat+100));
a       = str2double(dummy{2});
if strcmp(lower(dummy{3}),'bohr'), a = a*constants.Bohr; end;
siesta.lattice.constant = a;

%% Lee 'label.STRUCT_OUT' y lo parte por líneas
file  = fileread(fstr);
split = strsplit(file,'\n');

%% Vectores primitivos
line = strsplit(strtrim(split{1}));
a1   = [str2double(line{1}) str2double(line{2}) str2double(line{3})]./a;
line = strsplit(strtrim(split{2}));
a2   = [str2double(line{1}) str2double(line{2}) str2double(line{3})]./a;
line = strsplit(strtrim(split{3}));
a3   = [str2double(line{1}) str2double(line{2}) str2double(line{3})]./a;
siesta.lattice.vectors.direct = [ a1; a2; a3; ];

%% Genera los vectores de la red recíproca
Vol = dot(a1,cross(a2,a3));
b1  = cross(a2,a3)./Vol;
b2  = cross(a3,a1)./Vol;
b3  = cross(a1,a2)./Vol;
siesta.lattice.vectors.reciprocal = [ b1; b2; b3; ];

%% Información de la supercelda
supercell.Ncells            = 1;
supercell.cells.coordinates = [ 1 1 1 ];
supercell.R                 = [ 0 0 0 ].';
if contains(logfile,'Internal auxiliary supercell:');
  found1   = strfind(logfile,'Internal auxiliary supercell:');
  found2   = strfind(logfile,'superc: Number of atoms');
  dummy    = strsplit(strtrim(logfile(found1+30:found2-1)));

  N1 = str2double(dummy{1});
  N2 = str2double(dummy{3});
  N3 = str2double(dummy{5});
  if (rem([N1 N2 N3],2) == 0), 
	disp('El número de celdas unidad en SIESTA es incorrecto');
	stop
  end
  supercell.Ncells             = N1*N2*N3;
  supercell.cells.coordinates  = [ N1 N2 N3 ];

  n3 = [ 0:1:(N3-1)/2 -(N3-1)/2:1:-1 ];
  n2 = [ 0:1:(N2-1)/2 -(N2-1)/2:1:-1 ];
  n1 = [ 0:1:(N1-1)/2 -(N1-1)/2:1:-1 ];

  pp = [];
  for i3 = 1:length(n3)
    for i2 = 1:length(n2)
      for i1 = 1:length(n1)
        pp = [ pp; n1(i1) n2(i2) n3(i3) ];
      end
    end
  end

  supercell.R = pp.';
end

siesta.lattice.supercell = supercell;

end
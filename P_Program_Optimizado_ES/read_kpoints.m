%%
%% Lee el fichero 'label.KP' y lo guarda en sub-estructuras.
%%
%% Nk  = nº de k-points
%%
%% kpoints( ik = 1:Nk ):  [ kx ky kz peso ]
%% Sub-estructuras:
%%   'cartesian' -> coordenadas cartesianas (en Angstrom^-1)
%%   'reciprocal_vectors' -> coordenadas en unidades de vectores recíprocos
%%
%% SIESTA escribe los k-points en coordenadas cartesianas y Bohr^(-1).
%% Escrito por J. Ferrer, octubre 2020
%%

function [ siesta ] = read_kpoints( siesta, constants, flags )

fk     = [ siesta.system_label '.KP' ];

%% Lee el fichero y lo parte por líneas
file = fileread(fk);
kp   = strsplit(file,'\n');

%% Lee la primera línea y extrae Nk
firstline         = strsplit(strtrim(kp{1}));
siesta.kpoints.Nk = str2double(firstline{1});

G = 2*pi*siesta.lattice.vectors.reciprocal'./siesta.lattice.constant;
G = inv(G');
Bohr = constants.Bohr;

%% Lee (kx, ky, kz) y el peso
for i = 1:siesta.kpoints.Nk
    k      = strsplit(strtrim(kp{i+1}));
    kx     = str2double(k{2})/Bohr;
    ky     = str2double(k{3})/Bohr;
    kz     = str2double(k{4})/Bohr;
    weight = str2double(k{5});
    k = [ kx ky kz ];
    siesta.kpoints.cartesian(i,1:4)  = [ k weight ];
    siesta.kpoints.reciprocal_vectors(i,1:4) = [ k*G weight ];
end

end
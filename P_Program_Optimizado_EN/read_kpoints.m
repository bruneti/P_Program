%%
%% This function reads file 'label.KP' and retrieves it in sub-structures
%%
%% Nk  = # of kpoints
%%
%% kpoints( ik = 1:Nk ):  [ kx ky kz weight ]
%% Sub-structure 
%%   'cartesian' contains cartesian coordinates in Angstrom^-1)
%%   'reciprocal_vectors' contains coordinates in units of reciprocal lattice vectors
%%
%% Siesta writes kpoints in cartesian coordinates and Bohr^(-1).
%% Written by J. Ferrer, October 2020
%%

function [ siesta ] = read_kpoints( siesta, constants, flags )

fk     = [ siesta.system_label '.KP' ];

%% Read file & split it by lines
file = fileread(fk);
kp   = strsplit(file,'\n');

%% Read first line and extract no_s & no_u
firstline         = strsplit(strtrim(kp{1}));
siesta.kpoints.Nk = str2double(firstline{1});

G = 2*pi*siesta.lattice.vectors.reciprocal'./siesta.lattice.constant;
G = inv(G');
Bohr = constants.Bohr;

%% Read (kx, ky, kz) and weight
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
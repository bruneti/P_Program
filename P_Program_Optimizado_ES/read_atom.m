%%
%% Lee el fichero 'label.STRUCT_OUT', que contiene la información atómica.
%%
%% Rellena las sub-estructuras:
%%   Natoms = nº de átomos en la celda unidad (ojo: NO la supercelda)
%%
%%   atom(ia = 1: natom), que contiene:
%%       species_number = número fdf de la especie en SIESTA
%%       Z              = número atómico
%%       coords         = coordenadas cartesianas en unidades de la constante de red
%%
%% Escrito por D. Visontai & J. Ferrer, ~ 2017-2020
%% Autodetección de offset y columnas añadida para compatibilidad
%% con distintas versiones de SIESTA.
%%

function [ siesta ] = read_atom( siesta, constants, flags )

fstr = [ siesta.system_label '.STRUCT_OUT' ];

%% Leer fichero y partir por líneas
file  = fileread(fstr);
split = strsplit(file, '\n');

%% === 1. Localizar la línea de Natoms =================================
%  La primera línea que contenga un único entero positivo es Natoms.
%  Los vectores de red tienen siempre 3 tokens, así que no hay confusión.

natoms_idx = [];
for k = 1 : min(10, numel(split))
    tok = strsplit(strtrim(split{k}));
    if numel(tok) == 1
        n = str2double(tok{1});
        if ~isnan(n) && n == floor(n) && n > 0
            natoms_idx = k;
            break;
        end
    end
end

if isempty(natoms_idx)
    error('read_atom: no se encontró la línea de Natoms en %s.\n  Revisa el formato del fichero.', fstr);
end

siesta.Natoms.ucell = str2double(strtrim(split{natoms_idx}));
siesta.Natoms.scell = siesta.lattice.supercell.Ncells * siesta.Natoms.ucell;

%% === 2. Detectar el orden de columnas (species, Z) o (Z, species) ===========
%  El STRUCT_OUT puede escribir las columnas en dos órdenes según la versión
%  de SIESTA:
%    Formato A (antiguo):  species_idx   Z   xf  yf  zf
%    Formato B (nuevo):    Z   species_idx   xf  yf  zf
%
%  Heurística: el índice de especie es un entero en [1, Nspecies];
%  el número atómico Z suele ser > Nspecies para sistemas reales.
%  Leemos la primera línea de átomo y comprobamos qué columna cae dentro
%  del rango válido [1, Nspecies].

first_atom_line = strsplit(strtrim(split{natoms_idx + 1}));
c1 = round(str2double(first_atom_line{1}));
c2 = round(str2double(first_atom_line{2}));

c1_ok = (c1 >= 1) && (c1 <= siesta.Nspecies);
c2_ok = (c2 >= 1) && (c2 <= siesta.Nspecies);

if c1_ok && ~c2_ok
    species_col = 1;  Z_col = 2;
elseif c2_ok && ~c1_ok
    species_col = 2;  Z_col = 1;
elseif c1_ok && c2_ok
    % Ambas dentro de rango: elegir la menor (más probable que sea especie)
    if c1 <= c2
        species_col = 1;  Z_col = 2;
    else
        species_col = 2;  Z_col = 1;
    end
    warning('read_atom: columnas de especie y Z ambiguas (col1=%d, col2=%d, Nspecies=%d). Se asume species_col=%d.', ...
            c1, c2, siesta.Nspecies, species_col);
else
    error(['read_atom: no se puede determinar la columna de especie en %s.\n' ...
           '  Primera línea de átomo: "%s"\n  col1=%d, col2=%d, Nspecies=%d.\n' ...
           '  Comprueba el formato del STRUCT_OUT.'], ...
          fstr, strjoin(first_atom_line), c1, c2, siesta.Nspecies);
end

%% === 3. Leer átomos de la unit cell =====================================
R = siesta.lattice.vectors.direct.';

for i = 1 : siesta.Natoms.ucell
    tok = strsplit(strtrim(split{natoms_idx + i}));
    atom(i).species_number = str2double(tok{species_col});
    atom(i).Z              = str2double(tok{Z_col});
    atom(i).species_symbol = siesta.species(atom(i).species_number).symbol;
    atom(i).cell           = siesta.lattice.supercell.R(:,1);
    coords_frac(1)         = str2double(tok{3});
    coords_frac(2)         = str2double(tok{4});
    coords_frac(3)         = str2double(tok{5});
    atom(i).coords_ucell   = R * coords_frac.';
    atom(i).coords_scell   = R * coords_frac.';
    atom(i).PAO.PAOs       = siesta.species(atom(i).species_number).PAO.PAOs;
    atom(i).PAO.nPAO       = siesta.species(atom(i).species_number).PAO.nPAO;
    atom(i).PAO.ordering   = siesta.species(atom(i).species_number).PAO.ordering;
end

%% === 4. Expandir a la supercelda ==========================================
siesta.lattice.supercell.cells.atoms(1,:) = 1 : siesta.Natoms.ucell;

for ic = 1 : siesta.lattice.supercell.Ncells - 1
    siesta.lattice.supercell.cells.atoms(ic+1,:) = ic*siesta.Natoms.ucell + (1:siesta.Natoms.ucell)';
    for i = 1 : siesta.Natoms.ucell
        iatom                      = siesta.Natoms.ucell * ic + i;
        atom(iatom).species_number = atom(i).species_number;
        atom(iatom).species_symbol = atom(i).species_symbol;
        atom(iatom).Z              = atom(i).Z;
        atom(iatom).cell           = siesta.lattice.supercell.R(:, ic+1);
        atom(iatom).coords_ucell   = atom(i).coords_ucell;
        atom(iatom).coords_scell   = R * atom(iatom).cell + atom(i).coords_ucell;
        atom(iatom).PAO.PAOs       = atom(i).PAO.PAOs;
        atom(iatom).PAO.nPAO       = atom(i).PAO.nPAO;
        atom(iatom).PAO.ordering   = atom(i).PAO.ordering;
    end
end

siesta.atom = atom;

end

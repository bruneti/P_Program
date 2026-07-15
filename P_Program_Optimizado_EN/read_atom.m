%%
%% This function reads file 'label.STRUCT_OUT'.
%% The file contains the relevant atomic info.
%%
%% Output data are the sub-structures:
%%   Natoms = # of atoms in the unit cell (mind: NOT the supercell)
%%
%%   atom(ia = 1: natom), that contains the sub-sub-structures:
%%       species_number = siesta fdf number for the species
%%       Z              = Atomic number
%%       coords         = cartesian coordinates in units of the lattice constant
%%
%% Written by D. Visontai & J. Ferrer, ~ 2017-2020
%% Auto-detection of offset and columns added for compatibility with
%% different SIESTA versions.
%%

function [ siesta ] = read_atom( siesta, constants, flags )

fstr = [ siesta.system_label '.STRUCT_OUT' ];

%% Read the file and split it into lines
file  = fileread(fstr);
split = strsplit(file, '\n');

%% ==1. Locate the Natoms line ======================================
%  The first line containing a single positive integer is Natoms.
%  The lattice vectors always have 3 tokens, so there is no ambiguity.

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
    error('read_atom: could not find the Natoms line in %s.\n  Check the file format.', fstr);
end

siesta.Natoms.ucell = str2double(strtrim(split{natoms_idx}));
siesta.Natoms.scell = siesta.lattice.supercell.Ncells * siesta.Natoms.ucell;

%% == 2. Detect the column order (species, Z) or (Z, species) ==============
%  The STRUCT_OUT may write the columns in two orders depending on the
%  SIESTA version:
%    Format A (old):  species_idx   Z   xf  yf  zf
%    Format B (new):  Z   species_idx   xf  yf  zf
%
%  Heuristic: the species index is an integer in [1, Nspecies]; the atomic
%  number Z is usually > Nspecies for real systems. We read the first atom
%  line and check which column falls in the valid range [1, Nspecies].

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
    % Both in range: pick the smaller one (more likely to be the species)
    if c1 <= c2
        species_col = 1;  Z_col = 2;
    else
        species_col = 2;  Z_col = 1;
    end
    warning('read_atom: ambiguous species and Z columns (col1=%d, col2=%d, Nspecies=%d). Assuming species_col=%d.', ...
            c1, c2, siesta.Nspecies, species_col);
else
    error(['read_atom: cannot determine the species column in %s.\n' ...
           '  First atom line: "%s"\n  col1=%d, col2=%d, Nspecies=%d.\n' ...
           '  Check the STRUCT_OUT format.'], ...
          fstr, strjoin(first_atom_line), c1, c2, siesta.Nspecies);
end

%% == 3. Read the unit-cell atoms ========================================
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

%% ==== 4. Expand to the super-cell ==========================================
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

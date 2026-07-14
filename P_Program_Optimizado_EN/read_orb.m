%%
%% This function reads file label.ORB_INDX and retrieves the sub-structures:
%%
%%
%% Written by J. Ferrer, October 2020
%%

function [ siesta ] = read_orb( siesta, flags )

forb    = [ siesta.system_label '.ORB_INDX' ];
logfile = siesta.logfile;
Nspin   = siesta.spin.Nspin;

%% Read file & split it by lines
file    = fileread(forb);
orbindx = strsplit(file,'\n');

%% Read first line and extract nos & nou
firstline   = strsplit(strtrim(orbindx{1}));
siesta.Norbs.ucell = str2double(firstline{1});
siesta.Norbs.scell = str2double(firstline{2});
nos                = siesta.Norbs.scell;
nou                = siesta.Norbs.ucell;

if Nspin > 2
  siesta.Norbs.ucell = 2*siesta.Norbs.ucell;
  siesta.Norbs.scell = 2*siesta.Norbs.scell;
end

%% Read orbital info
for i = 1:nos
  o = strsplit(strtrim(orbindx{i+2}));
  if Nspin < 4 
    siesta.orb(i).index_scell    = str2double(o{1});
    siesta.orb(i).index_ucell    = str2double(o{16});
    siesta.orb(i).R(1:3)         = [ str2double(o{13}) str2double(o{14}) str2double(o{15}) ];
    siesta.orb(i).index_in_atom  = str2double(o{5});
    siesta.orb(i).atom           = str2double(o{2});
    siesta.orb(i).species        = str2double(o{3});
    siesta.orb(i).species_symbol = string(o{4});
    siesta.orb(i).n_number       = str2double(o{6});
    siesta.orb(i).l_number       = str2double(o{7});
    siesta.orb(i).M_number       = str2double(o{8});
    siesta.orb(i).symmetry       = string(o{11});
    siesta.orb(i).zeta_number    = str2double(o{9});
    siesta.orb(i).polarization   = string(o{10});
    siesta.orb(i).Rcutoff        = str2double(o{12});
  else
    up = 2*i-1;
    dn = 2*i;
    siesta.orb(up).index_scell    = 2*str2double(o{1})-1;
    siesta.orb(dn).index_scell    = 2*str2double(o{1});
    siesta.orb(up).index_ucell    = 2*str2double(o{16})-1;
    siesta.orb(dn).index_ucell    = 2*str2double(o{16});
    siesta.orb(up).R(1:3)         = [ str2double(o{13}) str2double(o{14}) str2double(o{15}) ];
    siesta.orb(dn).R(1:3)         = [ str2double(o{13}) str2double(o{14}) str2double(o{15}) ];
    siesta.orb(up).index_in_atom  = str2double(o{5});
    siesta.orb(dn).index_in_atom  = str2double(o{5});
    siesta.orb(up).atom           = str2double(o{2});
    siesta.orb(dn).atom           = str2double(o{2});
    siesta.orb(up).species        = str2double(o{3});
    siesta.orb(dn).species        = str2double(o{3});
    siesta.orb(up).species_symbol = string(o{4});
    siesta.orb(dn).species_symbol = string(o{4});
    siesta.orb(up).n_number       = str2double(o{6});
    siesta.orb(dn).n_number       = str2double(o{6});
    siesta.orb(up).l_number       = str2double(o{7});
    siesta.orb(dn).l_number       = str2double(o{7});
    siesta.orb(up).M_number       = str2double(o{8});
    siesta.orb(dn).M_number       = str2double(o{8});
    siesta.orb(up).symmetry       = string(o{11});
    siesta.orb(dn).symmetry       = string(o{11});
    siesta.orb(up).zeta_number    = str2double(o{9});
    siesta.orb(dn).zeta_number    = str2double(o{9});
    siesta.orb(up).polarization   = string(o{10});
    siesta.orb(dn).polarization   = string(o{10});
    siesta.orb(up).Rcutoff        = str2double(o{12});
    siesta.orb(dn).Rcutoff        = str2double(o{12});
    siesta.orb(up).spin           = string('up');
    siesta.orb(dn).spin           = string('down');
  end
end

for io = 1:siesta.Norbs.scell
  pp(io) = siesta.orb(io).atom;
  ppp(io) = siesta.orb(io).l_number;
end

for ia = 1:siesta.Natoms.scell
  siesta.atom(ia).orbs = find(pp(:) == ia);
  shells = siesta.atom(ia).PAO.PAOs(:,1);
  for ishell = 1:size(shells,1)
    siesta.atom(ia).l_number(ishell).orbs = intersect(find(pp(:) == ia),find(ppp(:) == shells(ishell)));
  end
end

for ic = 1:siesta.lattice.supercell.Ncells
  siesta.lattice.supercell.cells.orbs(ic,:) = (ic-1)*siesta.Norbs.ucell+(1:siesta.Norbs.ucell)';
end
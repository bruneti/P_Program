function[Delta_r_complete,LM]=goDelta(siesta)

%Getting every R-R' comparing also L, M, S, Z 
cuenta=1;
tamano=siesta.Norbs.ucell*siesta.Norbs.scell;
Delta_r1=zeros(tamano,1);
Delta_r2=zeros(tamano,1);
Delta_r3=zeros(tamano,1);
L1=zeros(tamano,1);
L2=zeros(tamano,1);
M1=zeros(tamano,1);
M2=zeros(tamano,1);
S1=zeros(tamano,1);
S2=zeros(tamano,1);
Z1=zeros(tamano,1);
Z2=zeros(tamano,1);

for    bb=1:siesta.Norbs.ucell           %bb for orbs on unit cell
   for cc=1:siesta.Norbs.scell           %cc for orbs on super cell
      
       aux1=siesta.orb(bb).atom;
       aux2=siesta.orb(cc).atom;
L1(cuenta,1)=siesta.orb(bb).l_number;
L2(cuenta,1)=siesta.orb(cc).l_number;
M1(cuenta,1)=siesta.orb(bb).M_number;
M2(cuenta,1)=siesta.orb(cc).M_number;
S1(cuenta,1)=siesta.orb(bb).species;
S2(cuenta,1)=siesta.orb(cc).species;
Z1(cuenta,1)=siesta.orb(bb).zeta_number;
Z2(cuenta,1)=siesta.orb(cc).zeta_number;

%l,m,species and z associated to each orb 1 and 2, compared on each case.
Delta_r1(cuenta,1)=(siesta.atom(aux2).coords_scell(1)-siesta.atom(aux1).coords_scell(1))*siesta.lattice.constant; %Delta_r component x
Delta_r2(cuenta,1)=(siesta.atom(aux2).coords_scell(2)-siesta.atom(aux1).coords_scell(2))*siesta.lattice.constant; %Delta_r component y
Delta_r3(cuenta,1)=(siesta.atom(aux2).coords_scell(3)-siesta.atom(aux1).coords_scell(3))*siesta.lattice.constant; %Delta_r component z

cuenta=cuenta+1;
      %end
      
   end
end
Bohr=0.529177210903;
Delta_r_complete=[Delta_r1,Delta_r2,Delta_r3]./Bohr;
LM=[L1,L2,M1,M2,S1,S2,Z1,Z2];

%We got all the R-R' ordered with its l, m, s, z, associated for each pair
%of orbitals.
end
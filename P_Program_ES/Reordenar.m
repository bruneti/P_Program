% Reordenar.m
% Pasa los elementos de la matriz P de forma de cadena (PLM_Delta) a forma
% matricial (PLMi_prueba). Con más de un k-point hace además la suma de Bloch
% sobre la supercelda usando post.k.
%
% NOTA: antes de ejecutar, carga el fichero de la dirección deseada y renómbralo:
%   load('PLM_Definitivo.mat');  % o PLM_Definitivo_z.mat en la versión optimizada
%   PLM_Delta = PLM_Delta;       % (renombra a PLM_Delta si tuviera otro nombre)

Nu=siesta.Norbs.ucell;
Ns=siesta.Norbs.scell;
p=Ns/Nu;
PLMi_prueba=reshape(PLM_Delta,[Ns,Nu]).'; % cadena -> forma matricial. Suficiente para un solo k-point.
PLM_k=zeros(Nu,Nu,p);

% Para más de un k-point.
for n=1:p

PLM_k(:,:,n)=PLMi_prueba(1:Nu,Nu*n-(Nu-1):Nu*n);



end

if p>1

 R=siesta.lattice.supercell.R;
 nsc = siesta.lattice.supercell.Ncells;
 ks=post.k;
 
 PLM_eikr=zeros(Nu,Nu,length(ks),p);
 for Rik=1:p
     for ik=1:length(ks)
         PLM_eikr(:,:,ik,Rik)=PLM_k(:,:,Rik)*exp(-1i*ks(ik,:)*R(:,Rik));
     end
 end
   PLM=sum(PLM_eikr,4);
   pk=1; % >>> ELECCIÓN DEL USUARIO: índice del k-point a inspeccionar. 1 = primer punto de la banda (Gamma).
   PLMi_prueba=PLM(:,:,pk);
   else
    ks=1;
end


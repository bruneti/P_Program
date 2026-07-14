% Reordenar.m
% Reshapes the P-matrix elements from string form (PLM_Delta) into matrix
% form (PLMi_prueba). With more than one k-point it also performs the Bloch
% sum over the super-cell using post.k.
%
% NOTE: before running, load the desired direction file and rename it, e.g.
%   load('PLM_Definitivo.mat');  % or PLM_Definitivo_z.mat in the optimized version
%   PLM_Delta = PLM_Delta;       % (rename to PLM_Delta if it has another name)

Nu=siesta.Norbs.ucell;
Ns=siesta.Norbs.scell;
p=Ns/Nu;
PLMi_prueba=reshape(PLM_Delta,[Ns,Nu]).'; % string -> matrix form. Enough for a single k-point.
PLM_k=zeros(Nu,Nu,p);

% For more than one k-point.
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
   pk=1; % >>> USER CHOICE: k-point index to inspect. 1 = first band point (Gamma).
   PLMi_prueba=PLM(:,:,pk);
   else
    ks=1;
end

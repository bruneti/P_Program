%Inferno_cero
Inferno_zero=zeros(3,4,4,7,7); % il,l1,l2,m1,m2
CL=Full_CL;
auxiliar=zeros(1,siesta.Norbs.ucell);
[v,k,~,~,SLZ]=goVl(siesta);
for i=1:siesta.Norbs.ucell

auxiliar(i)=siesta.orb(i).zeta_number;

end
zmax=max(auxiliar);




for il=-1:1
    for l1=0:3
        for l2=0:3
           
                 for m1=-l1:l1
                    for m2=-l2:l2
                       
                        
                        Inferno_zero(il+2,l1+1,l2+1,m1+l1+1,m2+l2+1)=(-1)^(il)*sqrt(4*pi/3)*(-1)^(m1+m2)*Extra_calc_2(il,l1,l2,m1,m2);      %*(-1)^(m1+m2)
                              
                        
                    end
                end
           
        end
    end
end
%save Inferno.mat Inferno_JiLM -mat



Inferno_extra=zeros(4,4,siesta.Nspecies,siesta.Nspecies,zmax,zmax); 

for l1=0:3
        for l2=0:3
           
                 for s1=1:siesta.Nspecies
                    for s2=1:siesta.Nspecies
                       for z1=1:zmax
                           for z2=1:zmax
                        
                        Inferno_extra(l1+1,l2+1,s1,s2,z1,z2)=Extra_calc(v,k,SLZ,l1,l2,s1,s2,z1,z2);         
                           end       
                       end
                    end
                end
           
        end
end
save('Inferno_000','Inferno_zero','Inferno_extra')

% Cálculo de PLM cuando R-R' = 0. Usamos Extra_calc y Extra_calc_2 para
% calcular ILLL y JiLM en este caso particular, y los tabulamos en
% Inferno_zero e Inferno_extra (guardados en Inferno_000.mat), reduciendo las
% iteraciones necesarias para calcular PLM.
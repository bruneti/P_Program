function [Illl]=goILLL(v,k,l1,l2,s1,s2,z1,z2,l3,SLZ,Delta_r)


slz1=[s1,l1,z1];
slz2=[s2,l2,z2];

% Para seleccionar el v correcto hay que fijar su especie, l y z. Esta parte
% simplemente obtiene el índice asociado al v deseado.

vl1=ismember(SLZ',slz1,'rows');

for i=1:length(vl1)
    if vl1(i)

        VL1=i;

    end
end
%VL1=2;
vl2=ismember(SLZ',slz2,'rows');

for i=1:length(vl2)
    if vl2(i)
       
        VL2=i;

    end
end



          Illl=trapz(k,k.^3.*conj(v(VL1,:)).*v(VL2,:).* sphbesselj(l3,k*norm(Delta_r)));  % Implementación de Illl del desarrollo teórico.
          % Como depende de Delta_r, es uno de los cálculos más costosos en
          % tiempo, al tener que iterar muchísimas veces.
          


end
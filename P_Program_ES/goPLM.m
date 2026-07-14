
function [PiLM]=goPLM(n,v,k,Delta_r,il,l1,l2,m1,m2,s1,s2,z1,z2,SLZ,Inferno_JiLM)
% Construcción general de PLM: el valor del elemento de la matriz P según los
% orbitales, la distancia (distinta de 0) entre ellos, etc.


x=Delta_r(1);
y=Delta_r(2);
z=Delta_r(3);


if l1>n
     fprintf('l1 demasiado grande')
     return
 end
 
 
 if l2>n
     fprintf('l2 demasiado grande')
     return
 end
 % Salvaguarda por si l1 o l2 exceden los límites esperados.


%%% Construcción principal de PLM %%%

lmax=l1+l2+1;

% lmax=l1+l2+1 es el máximo alcanzable en los cálculos siguientes, así que se
% fija en el valor máximo que puede dar contribución.


PiLM=0;
for lr=0:lmax 

    m=-lr:lr;
    for mr=m
        Yrml=Spherical_armonic_realB(x,y,z,lr,mr); % armónicos esféricos reales tabulados; esta función devuelve el valor correspondiente.
 

    PiLM=PiLM+(-1i)^(lr).*Inferno_JiLM(il+2,l1+1,l2+1,lr+1,m1+l1+1,m2+l2+1,mr+lr+1).*goILLL(v,k,l1,l2,s1,s2,z1,z2,lr,SLZ,Delta_r)*Yrml;  
  
    end
    
end

PiLM=(-1)^(il)*(-1)^(m1+m2)*4*pi*sqrt(4*pi/3)*PiLM; 

% Valor de la matriz P para los orbitales considerados, en el caso R-R' != 0.
end
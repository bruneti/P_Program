

function [JILM,CL1,CL2,CL3,CL4]=JiLM(CL,il,l1,l2,l3,m1,m2,m3)
% Calcula el JiLM de la teoría, donde se integran 4 armónicos esféricos.


if abs(m1)>l1
   fprintf('m1 y l1 no concuerdan')
   return
end
if il~=-1 && il~=0  && il~=1
fprintf('il solo puede valer -1, 0 o 1')
end
if abs(m2)>l2
    fprintf('m2 y l2 no concuerdan')
    return
end
if abs(m3)>l3
    fprintf('m3 y l3 no concuerdan')
return
end
% reglas para evitar errores

lmax=l1+l2; % con l en {0,1,2,3} basta; en este caso se toma Lmax=l1+l2
if l1==0
    CL1=1;
else 
    CL1=CL(1:2*l1+1,1:2*l1+1,l1);
end

if l2==0
    CL2=1;
else 
    CL2=CL(1:2*l2+1,1:2*l2+1,l2);
end
if l3==0
    CL3=1;
else 
    CL3=CL(1:2*l3+1,1:2*l3+1,l3);
end
l4=1;
m4=il;

CL4=CL(1:2*l4+1,1:2*l4+1,1);




P2=sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)*(2*l4+1))/(4*pi); 



m1c=m1+l1+1;  % pasa m a un índice válido para la matriz CL.
m2c=m2+l2+1;
m3c=m3+l3+1;
m4c=m4+l4+1;

P1S=0;
for mu1=-l1:l1
    for mu2=-l2:l2
        for mu3=-l3:l3
            for mu4=-l4:l4
                if mu1+mu2+mu3+mu4==0 
              for L=0:lmax
                  if (abs(l1-l2)<=L) && L<=l1+l2 && (abs(l3-l4)<=L) && L<=l3+l4 % acelera el método (regla de selección).
                P1S=P1S+(-1)^(mu1+mu2)*P2*CL1(m1c,mu1+l1+1)*CL2(m2c,mu2+l2+1)*CL3(m3c,mu3+l3+1)*CL4(m4c,mu4+l4+1)*(2*L+1)*Wigner3j([l1,l2,L],[0,0,0])*Wigner3j([l3,l4,L],[0,0,0])*Wigner3j([l1,l2,L],[mu1,mu2,mu3+mu4])*Wigner3j([l3,l4,L],[mu3,mu4,mu1+mu2]);% %mu1+mu2+mu1+mu2
               % Se incluye el factor (-1)^(mu1+mu2) por la solución de la
               % integral de 4 armónicos esféricos.
                  end
              end
               % En P1S incluimos P2 y todo el JiLM. CL transforma de armónicos
               % esféricos reales a complejos, para poder usar los símbolos 3j al
               % resolver la integral; el desarrollo usa armónicos reales porque
               % es la base que emplea SIESTA.
                end
            end
        end
    end
end


JILM=P1S;


end

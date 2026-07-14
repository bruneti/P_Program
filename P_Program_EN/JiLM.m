

function [JILM,CL1,CL2,CL3,CL4]=JiLM(CL,il,l1,l2,l3,m1,m2,m3)
%Here we will compute the JiLM from theory, where 4 spherical harmonics are
%integrated.


if abs(m1)>l1
   fprintf('m1 and l1 are mismatching')
   return
end
if il~=-1 && il~=0  && il~=1
fprintf('il only can have the values -1,0 and 1')
end
if abs(m2)>l2
    fprintf('m2 and l2 are mismatching')
    return
end
if abs(m3)>l3
    fprintf('m3 and l3 are mismatching')
return
end
%rules to avoid bugs

lmax=l1+l2; %We fix lmax to 3, as 0,1,2,3 should be enough values of l. In this scenario, Lmax=l1+l2 should be selected
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



m1c=m1+l1+1;  %This change the value of m to one understandable for the matrix CL.
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
                  if (abs(l1-l2)<=L) && L<=l1+l2 && (abs(l3-l4)<=L) && L<=l3+l4 %This should do the method even faster. Uncomment if needed. 
                P1S=P1S+(-1)^(mu1+mu2)*P2*CL1(m1c,mu1+l1+1)*CL2(m2c,mu2+l2+1)*CL3(m3c,mu3+l3+1)*CL4(m4c,mu4+l4+1)*(2*L+1)*Wigner3j([l1,l2,L],[0,0,0])*Wigner3j([l3,l4,L],[0,0,0])*Wigner3j([l1,l2,L],[mu1,mu2,mu3+mu4])*Wigner3j([l3,l4,L],[mu3,mu4,mu1+mu2]);% %mu1+mu2+mu1+mu2
               %I am including the product (-1)^(mu1+mu2) due to the
               % 4 spherical harmonics integral solution.
                  end
              end
               %On P1S we include P2, and all the JiLM. CL is used to
               %transform from real spherical armonics to complex spherical
               %armonics, in order to use 3jSymbols to solve the integral
               %and our developement includes this real spherical harmonics
               %just because SIESTA used that basis.
                end
            end
        end
    end
end


JILM=P1S;


end

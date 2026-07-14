
function [PiLM]=goPLM(n,v,k,Delta_r,il,l1,l2,m1,m2,s1,s2,z1,z2,SLZ,Inferno_JiLM)
%The general construction of PLM, the value of the element of the P matrix
%depending on the orbitals, the distance (different from 0) between them and so on. 


x=Delta_r(1);
y=Delta_r(2);
z=Delta_r(3);


if l1>n
     fprintf('l1 too large')
     return
 end
 
 
 if l2>n
     fprintf('l2 too large')
     return
 end
 %Just to avoid bugs, if l1 or l2 exceeds the limitations expected.
 

%%% The main PLM construction %%%

lmax=l1+l2+1;

%This lmax lmax=l1+l2+1 is the maximum reachable in the next calculations,
%so it is fixed at the maximum that could give us value.


PiLM=0;
for lr=0:lmax 

    m=-lr:lr;
    for mr=m
        Yrml=Spherical_armonic_realB(x,y,z,lr,mr); %I tabulated the values of the Real Spherical harmonics. With this function I simply take the corresponding value. 
 

    PiLM=PiLM+(-1i)^(lr).*Inferno_JiLM(il+2,l1+1,l2+1,lr+1,m1+l1+1,m2+l2+1,mr+lr+1).*goILLL(v,k,l1,l2,s1,s2,z1,z2,lr,SLZ,Delta_r)*Yrml;  
  
    end
    
end

PiLM=(-1)^(il)*(-1)^(m1+m2)*4*pi*sqrt(4*pi/3)*PiLM; 

%Getting the corresponding value of the P matrix to the orbitals we are working with, in case of R-R' not equal to 0.
end
function [CL]=Full_CL

lmax= 7; % lmax=7 es el máximo posible con el que trabajamos.

mtot=2*lmax+1;
CL=zeros(mtot,mtot,lmax+1);
for k=1:lmax
    mmax=2*k+1;
for ii=1:mmax
    for jj=1:mmax
        if ii==jj && ii>k+1 && k/2==floor(k/2)
        CL(ii,jj,k)=(-1)^(ii+1)/sqrt(2); 
        end
        if ii==jj && ii>k+1 && k/2~=floor(k/2)
        CL(ii,jj,k)=(-1)^(ii)/sqrt(2); 
        end

       if ii==jj  && ii<k+1
           CL(ii,jj,k)=1i/sqrt(2);
       end
       if ii>k+1 && jj<k+1 && ii-(k+1)==(k+1)-jj
           CL(ii,jj,k)=1/sqrt(2);
       end
       if ii<k+1 && jj>k+1 && jj-(k+1)==(k+1)-ii && k/2==floor(k/2)
           CL(ii,jj,k)=-1i*(-1)^(jj+1)/sqrt(2); 
       end
       if ii<k+1 && jj>k+1 && jj-(k+1)==(k+1)-ii && k/2~=floor(k/2)
           CL(ii,jj,k)=-1i*(-1)^(jj)/sqrt(2); 
       end
    end
end
CL(k+1,k+1,k)=1;
end
CL(1,1,lmax+1)=1;

% Devuelve CL: para l=0,1,2,3, todas las transformaciones de armónicos
% esféricos reales a armónicos esféricos complejos.
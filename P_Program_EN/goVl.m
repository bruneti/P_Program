function [v,k,n,R,SLZ]=goVl(siesta)

dr=zeros(1,siesta.Nspecies);
for ss=1:siesta.Nspecies

dr(ss)=min(siesta.species(ss).PAO.grid(:,2));

end   
[Dr,s]=min(dr); %fixing the min step.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SLZ recover in order %%%%%%%%
L_i=[];
z=[];
Species=[];
contador=1;
     for ss=1:siesta.Nspecies
         mm=length(siesta.species(ss).PAO.ordering);
         for ii=1:mm
             L_i(contador)=siesta.species(ss).PAO.ordering(ii);
             z(contador)=1;
             if ii>1 && L_i(contador)==L_i(contador-1) 
             z(contador)=z(contador-1)+1;
             end
             Species(contador)=ss;
             contador=contador+1;

         end
     end




SLZ=[Species;L_i;z];

%R depends on Species, L and Z coming from SIESTA. And we have to adjust to
%each orbital





n=0; %reseting any other value of n
for i=1:siesta.Nspecies
n=size(siesta.species(i).PAO.grid,1)+n;

end

Nrgrid=n;     %Number of PAOs under SIESTA notation and results. Needed for the next steps

  Ps=size(siesta.species(s).PAO.Rlz(1).Rlz,1); %The longest size of the Rlz from PAOs 
  R=zeros(Ps,Nrgrid);
  r=zeros(Ps,Nrgrid);

 %trapz |R*r|^2dr=1
 %trapz(R.^2.*r.^2.*dr)  must be =1 (normalization) Just for checking
 %afterwards
contador=1;
for s=1:siesta.Nspecies
 for ii=1:size(siesta.species(s).PAO.grid,1)
     R(:,contador)=siesta.species(s).PAO.Rlz(ii).Rlz(:,2);
     r(:,contador)=siesta.species(s).PAO.Rlz(ii).Rlz(:,1);
     contador=contador+1;
 end
end

% Building the new mesh equal for each PAO.

 g=zeros(1,Nrgrid);
 h=g;
 contador=1;
for ss=1:siesta.Nspecies
     for ii=1:size(siesta.species(ss).PAO.grid,1)
 g(contador)=min(min(siesta.species(ss).PAO.Rlz(ii).Rlz(:,1))); %Min of the mesh
 h(contador)=max(max(siesta.species(ss).PAO.Rlz(ii).Rlz(:,1))); %Max of the mesh
 contador=contador+1;
     end
end
 Malla=min(g):Dr:max(h); %fixing max and min of the final mesh and the step with Dr

 Rb=zeros(length(Malla),Nrgrid);
 for ii=1:Nrgrid
     F=griddedInterpolant(r(:,ii),R(:,ii));
     Rb(:,ii)=F(Malla);
 
 for i=1:length(Malla)
    for j=1:Nrgrid
    if Malla(i)> h(j)
     Rb(i,j)=0;
     end

    end
end
end
    %Rb=[Rb1;Rb2;Rb3;Rb4;Rb5];
R=Rb;
r=zeros(length(Malla),Nrgrid);
for i=1:Nrgrid

r(:,i)=Malla';
dr(:,i)=Dr;
end

%New mesh accomplished. Now every PAO is working with the same length and
%step.



 kmax=pi/(Dr); %Based on SIESTA manual
 k=0:0.1:kmax; %dk fixed as 0.1. It could be changed if needed, but it should be more than enough.
% k=linspace(0,kmax,length(r));
 v=zeros(Nrgrid,length(k));

 for kk=0:Nrgrid-1
     for jj=1:length(k)
         v(kk+1,jj)=sqrt(2/pi)*(-1i)^(L_i(kk+1))*trapz(r(:,kk+1),r(:,kk+1).^2.*sphbesselj(L_i(kk+1),k(jj)*r(:,kk+1)).*R(:,kk+1));
     end
 end    
% Here we have finally the v_l terms of the construction of P, associated to each possible orbital.

 end



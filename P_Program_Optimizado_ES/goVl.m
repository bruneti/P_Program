function [v,k,n,R,SLZ]=goVl(siesta)
% Construye los términos radiales v_l(k) de la matriz P, uno por PAO, sobre
% una malla radial común. Devuelve también la tabla de índices
% (especie,l,zeta) en SLZ.

dr=zeros(1,siesta.Nspecies);
for ss=1:siesta.Nspecies

dr(ss)=min(siesta.species(ss).PAO.grid(:,2));

end   
[Dr,s]=min(dr); % paso radial mínimo entre especies (fija el paso de malla)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SLZ recuperado en orden  %%%%%%%%
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




SLZ=[Species;L_i;z]; % fila 1: especie, fila 2: l, fila 3: zeta (por orbital)

% R depende de especie, l y zeta de SIESTA, ajustado a cada orbital.





n=0; % reinicio
for i=1:siesta.Nspecies
n=size(siesta.species(i).PAO.grid,1)+n;

end

Nrgrid=n;     % número de PAOs (notación SIESTA). Necesario más abajo.

  Ps=size(siesta.species(s).PAO.Rlz(1).Rlz,1); % mayor tamaño de Rlz entre los PAOs
  R=zeros(Ps,Nrgrid);
  r=zeros(Ps,Nrgrid);

 % Comprobación de normalización (opcional): trapz(R.^2.*r.^2.*dr) debe ser 1.
contador=1;
for s=1:siesta.Nspecies
 for ii=1:size(siesta.species(s).PAO.grid,1)
     R(:,contador)=siesta.species(s).PAO.Rlz(ii).Rlz(:,2);
     r(:,contador)=siesta.species(s).PAO.Rlz(ii).Rlz(:,1);
     contador=contador+1;
 end
end

% Construye una malla común compartida por todos los PAOs.

 g=zeros(1,Nrgrid);
 h=g;
 contador=1;
for ss=1:siesta.Nspecies
     for ii=1:size(siesta.species(ss).PAO.grid,1)
 g(contador)=min(min(siesta.species(ss).PAO.Rlz(ii).Rlz(:,1))); % mínimo de la malla
 h(contador)=max(max(siesta.species(ss).PAO.Rlz(ii).Rlz(:,1))); % máximo de la malla
 contador=contador+1;
     end
end
 Malla=min(g):Dr:max(h); % malla final: min..max con paso Dr

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

% Malla común lista: todos los PAOs comparten longitud y paso.



 kmax=pi/(Dr); % del manual de SIESTA
 k=0:0.1:kmax; % >>> ELECCIÓN DEL USUARIO (rara vez): paso en k, dk = 0.1 (suele sobrar)
% k=linspace(0,kmax,length(r));
 v=zeros(Nrgrid,length(k));

 for kk=0:Nrgrid-1
     for jj=1:length(k)
         v(kk+1,jj)=sqrt(2/pi)*(-1i)^(L_i(kk+1))*trapz(r(:,kk+1),r(:,kk+1).^2.*sphbesselj(L_i(kk+1),k(jj)*r(:,kk+1)).*R(:,kk+1));
     end
 end    
% Términos v_l(k) de la construcción de P, uno por cada orbital posible.

 end



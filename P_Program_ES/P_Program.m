% ========================= EJECUTABLE PRINCIPAL =========================
% Construye los elementos de la matriz P (operador momento) para UNA
% dirección de emisión, elegida mediante 'il' más abajo.
%
% Pipeline:  Cargar_SIESTA -> P_Program -> Reordenar -> Calculo_Final
% Tras este script, ejecuta Reordenar (para poner P en forma matricial) y
% luego Calculo_Final, que devuelve <|P|>, Mu, Gamma_R y la vida media.
% ------------------------------------------------------------------------

%load('siesta.mat')
[v,k,n,R,SLZ]=goVl(siesta);           % términos v_l y malla radial
[Delta_r_complete,LM]=goDelta(siesta);% vectores R-R' + (l,m,s,z) por par de orbitales
Inferno     % construye Inferno_JiLM (JiLM tabulado, tensor 7-D)
Inferno_cero   % construye Inferno_zero e Inferno_extra (necesarios si R-R' = 0)

l1=LM(:,1);
l2=LM(:,2); 
m1=LM(:,3);
m2=LM(:,4);
s1=LM(:,5);
s2=LM(:,6);
z1=LM(:,7);
z2=LM(:,8);
% Números cuánticos tomados de cada par de orbitales para el cálculo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
il= 0; % >>> ELECCIÓN DEL USUARIO: dirección de emisión.  -1 -> y,  0 -> z,  1 -> x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PLM_Delta=zeros(length(Delta_r_complete),1);
tic % Cronometra este bucle.
for i=1:length(Delta_r_complete)
    if Delta_r_complete(i,1)==0 && Delta_r_complete(i,2)==0 && Delta_r_complete(i,3)==0
        PLM_Delta(i)=Inferno_zero(il+2,l1(i)+1,l2(i)+1,m1(i)+l1(i)+1,m2(i)+l2(i)+1)*Inferno_extra(l1(i)+1,l2(i)+1,s1(i),s2(i),z1(i),z2(i));  % elemento de matriz para R-R' = 0
    else
PLM_Delta(i)=goPLM(n,v,k,Delta_r_complete(i,:),il,l1(i),l2(i),m1(i),m2(i),s1(i),s2(i),z1(i),z2(i),SLZ,Inferno_JiLM); % elemento de matriz para R-R' != 0
    end
end

%save('PLM_x','PLM_Delta')
save('PLM_Definitivo.mat','PLM_Delta') % resultado final: elementos de P (forma de cadena)
toc
% PLM_Definitivo.mat contiene los elementos de la matriz P para la dirección
% elegida (il) en forma de cadena. Usa Reordenar para pasarlos a forma
% matricial.

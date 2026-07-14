% ========================= MAIN EXECUTABLE ==============================
% Builds the P-matrix elements (momentum operator) for ONE emission
% direction, selected through 'il' below.
%
% Pipeline:  Cargar_SIESTA -> Pain_Program -> Reordenar -> Calculo_Final
% After this script, run Reordenar (to reshape P into matrix form) and then
% Calculo_Final, which returns <|P|>, Mu, Gamma_R and the lifetime.
% ------------------------------------------------------------------------

%load('siesta.mat')
[v,k,n,R,SLZ]=goVl(siesta);            % v_l terms and radial mesh
[Delta_r_complete,LM]=goDelta(siesta); % R-R' vectors + (l,m,s,z) per orbital pair
Inferno        % builds Inferno_JiLM (tabulated JiLM, 7-D tensor)
Inferno_cero   % builds Inferno_zero and Inferno_extra (needed when R-R' = 0)

l1=LM(:,1);
l2=LM(:,2);
m1=LM(:,3);
m2=LM(:,4);
s1=LM(:,5);
s2=LM(:,6);
z1=LM(:,7);
z2=LM(:,8);
% Quantum numbers taken from each orbital pair for the calculation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
il= 0; % >>> USER CHOICE: emission direction.  -1 -> y,  0 -> z,  1 -> x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PLM_Delta=zeros(length(Delta_r_complete),1);
tic % Time this loop.
for i=1:length(Delta_r_complete)
    if Delta_r_complete(i,1)==0 && Delta_r_complete(i,2)==0 && Delta_r_complete(i,3)==0
        PLM_Delta(i)=Inferno_zero(il+2,l1(i)+1,l2(i)+1,m1(i)+l1(i)+1,m2(i)+l2(i)+1)*Inferno_extra(l1(i)+1,l2(i)+1,s1(i),s2(i),z1(i),z2(i));  % matrix element for R-R' = 0
    else
PLM_Delta(i)=goPLM(n,v,k,Delta_r_complete(i,:),il,l1(i),l2(i),m1(i),m2(i),s1(i),s2(i),z1(i),z2(i),SLZ,Inferno_JiLM); % matrix element for R-R' != 0
    end
end

%save('PLM_x','PLM_Delta')
save('PLM_Definitivo.mat','PLM_Delta') % final result: P-matrix elements (string form)
toc
% PLM_Definitivo.mat holds the P-matrix elements for the chosen direction
% (il) in string form. Use Reordenar to turn them into matrix form.

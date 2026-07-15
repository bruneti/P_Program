% Calculo_Final.m
% Paso final: a partir de la matriz P (PLMi_prueba) y los autoestados calcula
% el elemento de matriz de transición <psi|P|psi>, el dipolo de transición Mu,
% la tasa radiativa Gamma_R y la vida media asociada.
%
% Queremos calcular: Sum C_n,mu * <P> C_n,mu


% Cuando trabajamos con más de un k-point aparece esto:
R   = siesta.lattice.supercell.R;
nsc = siesta.lattice.supercell.Ncells;

k=[0,0,0];% >>> ELECCIÓN DEL USUARIO: k-point. Por defecto Gamma = [0 0 0]. %k=ks(pk,:);

eikr = reshape( exp(-1i*k*R), [ 1 1 nsc ] );



if siesta.spin.Nspin == 2 % Caso spin-polarized
H=siesta.Hsc;
Hk(:,:,1) = sum(bsxfun(@times,H(1).sp,eikr), 3);
Hk(:,:,2) = sum(bsxfun(@times,H(2).sp,eikr), 3);
Hk(:,:,2) = (Hk(:,:,2)+Hk(:,:,2)')/2;
Hk(:,:,1) = (Hk(:,:,1)+Hk(:,:,1)')/2;
O=siesta.Osc;
Ok = sum( bsxfun( @times, O, eikr), 3);
Ok=(Ok+Ok')/2;
[eigen_vector1,~]=eig(Hk(:,:,1),Ok);
[eigen_vector2,~]=eig(Hk(:,:,2),Ok);
% >>> ELECCIÓN DEL USUARIO: estados inicial y final más abajo.
% eigen_vector1 = spin up, eigen_vector2 = spin down.
 C_mu=eigen_vector1(:,142);   % >>> índice del estado inicial
 C_mup=eigen_vector1(:,143);  % >>> índice del estado final
else
 
H_eikr=zeros(Nu,Nu,length(ks),p);
    H_k=siesta.Hsc;
    Hk=H_k;
    if p>1 || length(ks)>1
 for Rik=1:p
     for ik=1:length(ks)
         H_eikr(:,:,ik,Rik)=H_k(:,:,Rik)*exp(-1i*ks(ik,:)*R(:,Rik));%*exp(-1i*ks(ik,:)*R(:,Rik));
     end
 end
Hsum= sum(H_eikr,4);
Hk=Hsum(:,:,pk); % pk se selecciona en Reordenar
    else      
    Hk = (Hk+Hk')/2;
    %end
O=siesta.Osc;
Ok = sum( bsxfun( @times, O, eikr), 3);
Ok=(Ok+Ok')/2;
   
[eigen_vector,~]=eig(Hk,Ok);

 C_mu=eigen_vector(:,142);   % >>> ELECCIÓN DEL USUARIO: índice del estado inicial
 C_mup=eigen_vector(:,143);  % >>> ELECCIÓN DEL USUARIO: índice del estado final
    end
end


hbar=4.135667696*10^(-15)*1.6*10^(-19)/(2*pi);
Bohr=0.529177210903; % 1 Bohr = 0.529... Angstrom
P=PLMi_prueba;
%P=(P+P')/2;
Prob=C_mu'*P*C_mup; % <psi|P|psi>, todavía no en unidades SI
Prob=Prob*10^(10)*hbar*Bohr; % convertido a unidades SI
me=9.1*10^(-31);
if min(size(post.eigen_energies))==1
d=post.eigen_energies;
else
d=diag(post.eigen_energies);
end
DE=(d(143)-d(142))*(1.6*10^(-19));% >>> ELECCIÓN DEL USUARIO: gap de energía. Los índices DEBEN coincidir con los estados de arriba.
% Aquí no hay ..1 ..2 como en los eigen_vectors. Con más de un k-point puede
% ser necesario usar eig(Hk(:,:...)) en lugar de post.eigen_energies.


Mu=1i*hbar*Prob/(me*DE); % Dipolo de transición (Cholsuk et al., hBN Defects Database, JPCC 2024, DOI:10.1021/acs.jpcc.4c03404)
eps_0=8.85*10^(-12);
c=3*10^8;
%1.85 si hBN 4 si WSe2   % >>> factor de índice de refracción del medio (informativo)
Gamma_R=1.85*constants.q^2*DE^3*(abs(Mu))^2/(3*pi*eps_0*hbar^4*c^3); % Tasa de transición radiativa (misma referencia que Mu)
vidamedia=1/Gamma_R; % Vida media asociada a Gamma_R

% ===================== RESULTADOS =====================
fprintf('\n');
fprintf('============================================\n');
fprintf('       RESULTADOS - Calculo_Final.m\n');
fprintf('============================================\n');
%fprintf('  <psi|P|psi> (SI)    = %e + %ei\n', real(Prob), imag(Prob));
fprintf('  Delta E             = %e J  (%.4f eV)\n', DE, DE/(1.6e-19));
fprintf('  Mu (dipolo trans.)  = %e + %ei m\n', real(Mu), imag(Mu));
fprintf('  |Mu|                = %e m\n', abs(Mu));
fprintf('  Gamma_R             = %e s^-1\n', Gamma_R);
fprintf('  Vida media          = %e s  (%.4f ns)\n', vidamedia, vidamedia*1e9);
fprintf('============================================\n');

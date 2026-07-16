% Calculo_Final.m
% Final step: from the P matrix (PLMi_prueba) and the eigenstates it computes
% the transition matrix element <psi|P|psi>, the transition dipole Mu, the
% radiative rate Gamma_R and the associated lifetime.
%
% We want to calculate: Sum C_n,mu * <P> C_n,mu

% When we work with more than one k-point we get this:
R   = siesta.lattice.supercell.R;
nsc = siesta.lattice.supercell.Ncells;

k=[0,0,0]; % >>> USER CHOICE: k-point. Default Gamma = [0 0 0]. %k=ks(pk,:);

eikr = reshape( exp(-1i*k*R), [ 1 1 nsc ] );



if siesta.spin.Nspin == 2
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
% >>> USER CHOICE: initial and final states below.
% eigen_vector1 = spin up, eigen_vector2 = spin down.
 C_mu=eigen_vector1(:,4);   % >>> initial state index
 C_mup=eigen_vector1(:,5);  % >>> final state index
else

H_eikr=zeros(Nu,Nu,length(ks),p);
    H_k=siesta.Hsc;
    Hk=H_k;
    if p>1 || length(ks)>1
 for Rik=1:p
     for ik=1:length(ks)
         H_eikr(:,:,ik,Rik)=H_k(:,:,Rik)*exp(-1i*ks(ik,:)*R(:,Rik));
     end
 end
Hsum= sum(H_eikr,4);
Hk=Hsum(:,:,pk); % pk is selected in Reordenar
    else
    Hk = (Hk+Hk')/2;
    %end
O=siesta.Osc;
Ok = sum( bsxfun( @times, O, eikr), 3);
Ok=(Ok+Ok')/2;

[eigen_vector,~]=eig(Hk,Ok);

% >>> USER CHOICE: initial and final states.
 C_mu=eigen_vector(:,4);   % >>> initial state index
 C_mup=eigen_vector(:,5);  % >>> final state index
    end
end


hbar=4.135667696*10^(-15)*1.6*10^(-19)/(2*pi);
Bohr=0.529177210903; % 1 Bohr = 0.529... Angstrom
P=PLMi_prueba;
%P=(P+P')/2;
Prob=C_mu'*P*C_mup; % <psi|P|psi>, not yet in SI units
Prob=Prob*10^(10)*hbar*Bohr; % converted to SI units
me=9.1*10^(-31);
if min(size(post.eigen_energies))==1
d=post.eigen_energies;
else
d=diag(post.eigen_energies);
end
DE=(d(5)-d(4))*(1.6*10^(-19)); % >>> USER CHOICE: energy gap. Indices MUST match the states above.
% Note there are no ..1 ..2 here as there were for the eigen_vectors. With
% more than one k-point it may be necessary to use eig(Hk(:,:...)) instead
% of post.eigen_energies.


Mu=1i*hbar*Prob/(me*DE); % Transition dipole (Cholsuk et al., hBN Defects Database, JPCC 2024, DOI:10.1021/acs.jpcc.4c03404)
eps_0=8.85*10^(-12);
c=3*10^8;
%1.85 if hBN 4 if WSe2   % >>> refractive-index factor of the medium (informational)
Gamma_R=1.85*constants.q^2*DE^3*(abs(Mu))^2/(3*pi*eps_0*hbar^4*c^3); % Radiative transition rate (same reference as Mu)
vidamedia=1/Gamma_R; % Mean lifetime associated with Gamma_R
% ===================== RESULTS =====================
fprintf('\n');
fprintf('============================================\n');
fprintf('       RESULTS - Calculo_Final.m\n');
fprintf('============================================\n');
%fprintf('  <psi|P|psi> (SI)    = %e + %ei\n', real(Prob), imag(Prob));
fprintf('  Delta E             = %e J  (%.4f eV)\n', DE, DE/(1.6e-19));
fprintf('  Mu (dipolo trans.)  = %e + %ei m\n', real(Mu), imag(Mu));
fprintf('  |Mu|                = %e m\n', abs(Mu));
fprintf('  Gamma_R             = %e s^-1\n', Gamma_R);
fprintf(' Lifetime             = %e s  (%.4f ns)\n', vidamedia, vidamedia*1e9);
fprintf('============================================\n');

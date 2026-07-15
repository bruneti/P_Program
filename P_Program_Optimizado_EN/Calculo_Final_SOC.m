% =========================================================================
% Calculo_Final_SOC.m  -  FINAL VERSION
% Version adapted for Nspin = 8 (Spin-Orbit Coupling, non-collinear SIESTA)
%
% ==== Differences with respect to Calculo_Final.m (no-SOC / spin-polarised) ==
%
%   1. siesta.Hsc is a SINGLE spinor matrix of dimension (2*Nu x 2*Nu x nsc).
%      read_hodm.m assembles it via kron products, so no spin index needed.
%
%   2. Eigenvectors are spinors of length 2*Nu with INTERLEAVED ordering:
%        C(1,3,5,...) = spin-up   components  (odd  indices)
%        C(2,4,6,...) = spin-down components  (even indices)
%      This ordering comes from  H = kron(H_spatial, tau_spin)  in read_hodm.
%
%   3. PLMi_prueba (from Pain_Program_3il + Reordenar) has size (2*Nu x 2*Nu),
%      but its structure is NOT  P_spatial x I_2 . Because goDelta assigns
%      the same spatial integral to ALL four spin combinations of the same
%      orbital pair (up-up, up-dn, dn-up, dn-dn), the full matrix is:
%
%            PLMi_prueba  =  P_spatial  x  [[1,1],[1,1]]
%
%      Using  C' * PLMi_prueba * C  directly would include spurious up-dn
%      cross terms that the momentum operator (a spin scalar) must NOT have.
%      The correct spinor contraction is:
%
%         <f|P|i>  =  C_f_up' * P_sp * C_i_up  +  C_f_dn' * P_sp * C_i_dn
%
%      where  P_sp = PLMi_prueba(1:2:end, 1:2:end)  is the up-up block
%      (which equals the dn-dn block and equals the spatial matrix).
%
%   4. Kramers degeneracy: with time-reversal symmetry, states come in
%      degenerate pairs. The radiative rate FROM a final state must SUM
%      over BOTH degenerate initial states:
%
%         Gamma_R(emp1) = A * ( |<occ1|P|emp1>|^2 + |<occ2|P|emp1>|^2 )
%
%      By time-reversal symmetry,  Gamma_R(emp2) = Gamma_R(emp1), so only
%      ONE final state of the pair needs to be computed.
%
% === Workflow ============================================================
%   Pain_Program_3il  ->  Reordenar  ->  this script
%   (identical to the no-SOC workflow; just swap the final script)
%
% ====== User input =======================================================
%   Set idx_occ_1, idx_occ_2, idx_emp_1, idx_emp_2 to the Kramers-pair
%   indices you want to evaluate. The script will verify the pairs are
%   actually degenerate before proceeding.
% =========================================================================

%% -- Bloch phase (Gamma point; change k if needed) -----------------------
R    = siesta.lattice.supercell.R;
nsc  = siesta.lattice.supercell.Ncells;
k    = [0, 0, 0];
eikr = reshape( exp(-1i*k*R), [1 1 nsc] );

%% -- Build H(k) and O(k) from the spinor supercell matrices --------------
Hk = sum( bsxfun(@times, siesta.Hsc, eikr), 3 );
Hk = (Hk + Hk') / 2;

Ok = sum( bsxfun(@times, siesta.Osc, eikr), 3 );
Ok = (Ok + Ok') / 2;

%% -- Diagonalise (generalised Hermitian problem, eigs come sorted asc) ---
[eigen_vector, eigen_vals_mat] = eig(Hk, Ok);
eigen_energies = diag(eigen_vals_mat);

%% =========================================================================
%  USER INPUT
% =========================================================================
idx_occ_1 = 647;   % >>> USER CHOICE: Kramers partner 1 of the occupied pair
idx_occ_2 = 648;   % >>> USER CHOICE: Kramers partner 2 of the occupied pair
idx_emp_1 = 649;   % >>> USER CHOICE: final state for which the rate is computed
idx_emp_2 = 650;   % >>> (only used for the degeneracy check; rate is the same)

%% -- Sanity check: verify Kramers degeneracy ----------------------------
tol_deg = 1e-4;   % eV - tolerance for "degenerate"
E_occ1  = eigen_energies(idx_occ_1);
E_occ2  = eigen_energies(idx_occ_2);
E_emp1  = eigen_energies(idx_emp_1);
E_emp2  = eigen_energies(idx_emp_2);

fprintf('Occupied pair: E(%d)=%.6f eV, E(%d)=%.6f eV, split=%.2e eV\n', ...
    idx_occ_1, E_occ1, idx_occ_2, E_occ2, abs(E_occ1-E_occ2));
fprintf('Empty    pair: E(%d)=%.6f eV, E(%d)=%.6f eV, split=%.2e eV\n', ...
    idx_emp_1, E_emp1, idx_emp_2, E_emp2, abs(E_emp1-E_emp2));

if abs(E_occ1 - E_occ2) > tol_deg
    warning('Occupied states are not degenerate within %.0e eV. Check indices.', tol_deg);
end
if abs(E_emp1 - E_emp2) > tol_deg
    warning('Empty states are not degenerate within %.0e eV. Check indices.', tol_deg);
end

%% -- Extract spin components from the spinor eigenvectors ---------------
% Interleaved ordering: odd rows = spin up, even rows = spin down
idx_up = 1:2:size(eigen_vector,1);
idx_dn = 2:2:size(eigen_vector,1);

C_occ1_up = eigen_vector(idx_up, idx_occ_1);
C_occ1_dn = eigen_vector(idx_dn, idx_occ_1);

C_occ2_up = eigen_vector(idx_up, idx_occ_2);
C_occ2_dn = eigen_vector(idx_dn, idx_occ_2);

C_emp1_up = eigen_vector(idx_up, idx_emp_1);
C_emp1_dn = eigen_vector(idx_dn, idx_emp_1);

%% -- Extract spatial P matrix from PLMi_prueba --------------------------
% The up-up block is the true spatial P. Any off-diagonal spin block in
% PLMi_prueba is a spurious artefact of goDelta not distinguishing spin.
P_sp = PLMi_prueba(idx_up, idx_up);

%% -- Spinor matrix elements (sum over spin components) -------------------
% <f|P|i> = <f_up|P_sp|i_up> + <f_dn|P_sp|i_dn>
Prob_occ1_emp1 = C_occ1_up' * P_sp * C_emp1_up + C_occ1_dn' * P_sp * C_emp1_dn;
Prob_occ2_emp1 = C_occ2_up' * P_sp * C_emp1_up + C_occ2_dn' * P_sp * C_emp1_dn;

%% -- Physical constants (same as original Calculo_Final.m) --------------
hbar  = 4.135667696e-15 * 1.6e-19 / (2*pi);   % J*s
Bohr  = 0.529177210903;                        % 1 Bohr in Angstrom
me    = 9.1e-31;                               % kg
eps_0 = 8.85e-12;                              % F/m
c_luz = 3e8;                                   % m/s
q_e   = 1.6e-19;                               % C

%% -- Energy difference (use emp1 and mean of occupied pair) --------------
DE = ( E_emp1 - (E_occ1 + E_occ2)/2 ) * q_e;    % Joules
fprintf('\nEnergy gap DE = %.6f eV\n', DE / q_e);

%% -- Unit conversion (same factor as original Calculo_Final.m) ----------
unit_factor       = 1e10 * hbar * Bohr;
Prob_occ1_emp1_SI = Prob_occ1_emp1 * unit_factor;
Prob_occ2_emp1_SI = Prob_occ2_emp1 * unit_factor;

%% -- Transition dipole moments (Cholsuk et al., JPCC 2024) --------------
Mu_occ1_emp1 = 1i * hbar * Prob_occ1_emp1_SI / (me * DE);
Mu_occ2_emp1 = 1i * hbar * Prob_occ2_emp1_SI / (me * DE);

fprintf('|Mu(occ1->emp1)| = %.4e  C*m\n', abs(Mu_occ1_emp1));
fprintf('|Mu(occ2->emp1)| = %.4e  C*m\n', abs(Mu_occ2_emp1));

%% -- Radiative rate - SAME FORMULA as Calculo_Final.m, with sum over ----
%%    degenerate initial states ------------------------------------------
% Gamma_R(emp1) = A * ( |Mu(occ1->emp1)|^2 + |Mu(occ2->emp1)|^2 )
% Gamma_R(emp2) = Gamma_R(emp1)  (Kramers symmetry)
%
% NOTE: to include a refractive-index (medium) correction, MULTIPLY the
%       rate by n (n ~ 1.85 for hBN, ~ 4 for WSe2). Do NOT divide.  >>> USER CHOICE
%       It is NOT included below to stay faithful to the original formula.

prefactor = 4 * q_e^2 * DE^3 / (3 * pi * eps_0 * hbar^4 * c_luz^3);

Gamma_R  = prefactor * ( abs(Mu_occ1_emp1)^2 + abs(Mu_occ2_emp1)^2 );
lifetime = 1 / Gamma_R;

fprintf('\nGamma_R  = %.6e  s^-1\n', Gamma_R);
fprintf('Lifetime = %.6e  s\n',  lifetime);

%% -- Store results -------------------------------------------------------
results_SOC.idx_occ       = [idx_occ_1, idx_occ_2];
results_SOC.idx_emp       = [idx_emp_1, idx_emp_2];
results_SOC.DE_eV         = DE / q_e;
results_SOC.Mu_occ1_emp1  = Mu_occ1_emp1;
results_SOC.Mu_occ2_emp1  = Mu_occ2_emp1;
results_SOC.Gamma_R       = Gamma_R;
results_SOC.lifetime_s    = lifetime;

fprintf('\nResults stored in  results_SOC\n');

% -------------------------------------------------------------------------
% Run this script three times (one per spatial direction x, y, z) by
% loading the corresponding PLM_Definitivo_?.mat and passing through
% Reordenar first. Compare the three lifetimes to identify the emission
% polarisation, exactly as in the original hBN study.
% -------------------------------------------------------------------------

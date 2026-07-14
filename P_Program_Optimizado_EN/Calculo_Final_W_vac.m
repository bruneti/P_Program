% =========================================================================
% Calculo_Final_W_vac.m
%
% Adaptation of Calculo_Final_SOC.m for the particular W_vac case:
%
%   - The defect induces a net magnetic moment.
%   - Time-reversal symmetry is broken.
%   - There is NO Kramers degeneracy: every state is non-degenerate.
%
% DIFFERENCES with respect to Calculo_Final_SOC.m:
%
%   1. No sum over the Kramers pair. A single <f|P|i> is computed for the
%      state pair specified by the user.
%
%   2. The prior degeneracy check is removed. Instead, <S_z> of the states
%      involved is computed and the transition is verified to be
%      spin-conserving (optical selection rule).
%
%   3. <S_z>_i and <S_z>_f are reported explicitly, and a warning is issued
%      if the transition is spin-flip (suppressed) rather than computing
%      blindly.
%
% The rest of the logic (extraction of the spatial block of PLMi_prueba,
% spinor contraction, prefactor with n=4) stays identical.
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

%% -- Diagonalise --------------------------------------------------------
[eigen_vector, eigen_vals_mat] = eig(Hk, Ok);
eigen_energies = diag(eigen_vals_mat);

%% =========================================================================
%  USER INPUT
% =========================================================================
idx_occ = 642;   % >>> USER CHOICE: initial (occupied) state index
idx_emp = 643;   % >>> USER CHOICE: final (empty) state index

%% -- Compute <S_z> for the two states involved -------------------------
N_total = size(eigen_vector, 1);
N_orbs  = N_total / 2;
Sz      = kron(eye(N_orbs), [1 0; 0 -1]/2);

sz_occ = real( eigen_vector(:,idx_occ)' * Sz * eigen_vector(:,idx_occ) );
sz_emp = real( eigen_vector(:,idx_emp)' * Sz * eigen_vector(:,idx_emp) );

E_occ = eigen_energies(idx_occ);
E_emp = eigen_energies(idx_emp);

fprintf('========================================================\n');
fprintf('  W_vac calculation without Kramers pairs\n');
fprintf('========================================================\n');
fprintf('Initial state (i): index %d\n', idx_occ);
fprintf('  E = %.6f eV,  <S_z> = %+.4f hbar\n', E_occ, sz_occ);
fprintf('Final state   (f): index %d\n', idx_emp);
fprintf('  E = %.6f eV,  <S_z> = %+.4f hbar\n', E_emp, sz_emp);
fprintf('Energy difference: DE = %.6f eV\n\n', E_emp - E_occ);

%% -- Spin selection-rule check -----------------------------------------
threshold_sz = 0.1;   % below this we consider it "mixed"

if abs(sz_occ) < threshold_sz || abs(sz_emp) < threshold_sz
    fprintf('WARNING: at least one state has |<S_z>| < %.2f.\n', threshold_sz);
    fprintf('       The transition is of mixed character. The result can\n');
    fprintf('       be computed but must be interpreted with caution.\n\n');
elseif sign(sz_occ) ~= sign(sz_emp)
    fprintf('WARNING: the transition is SPIN-FLIP.\n');
    fprintf('       <S_z>_i and <S_z>_f have opposite signs.\n');
    fprintf('       The momentum operator does not induce spin flip, so\n');
    fprintf('       the result should come out close to zero.\n');
    fprintf('       Proceeding with the calculation anyway.\n\n');
else
    fprintf('OK: SPIN-CONSERVING transition (both in channel %s).\n\n', ...
            ternary(sz_occ > 0, 'UP', 'DOWN'));
end

%% -- Extract the spinor components --------------------------------------
idx_up = 1:2:N_total;
idx_dn = 2:2:N_total;

C_occ_up = eigen_vector(idx_up, idx_occ);
C_occ_dn = eigen_vector(idx_dn, idx_occ);

C_emp_up = eigen_vector(idx_up, idx_emp);
C_emp_dn = eigen_vector(idx_dn, idx_emp);

%% -- Extract the spatial block of PLMi_prueba --------------------------
% Same logic as in Calculo_Final_SOC.m: the four spin blocks are identical,
% so it is enough to extract the up-up block.
P_sp = PLMi_prueba(idx_up, idx_up);

%% -- Spinor matrix element ---------------------------------------------
% <f|P|i> = (C_f^up)' * P_sp * C_i^up + (C_f^dn)' * P_sp * C_i^dn
% The momentum operator is a spin scalar, so there are NO cross terms.
Prob = C_occ_up' * P_sp * C_emp_up + C_occ_dn' * P_sp * C_emp_dn;
Prob_up = C_occ_up' * P_sp * C_emp_up;
Prob_dn = C_occ_dn' * P_sp * C_emp_dn;

fprintf('\n--- Spinor decomposition ---\n');
fprintf('UP contribution:    %.4e + %.4ei\n', real(Prob_up), imag(Prob_up));
fprintf('DN contribution:    %.4e + %.4ei\n', real(Prob_dn), imag(Prob_dn));
fprintf('Sum:                %.4e + %.4ei\n', real(Prob),    imag(Prob));
fprintf('|UP|/|DN|:          %.4f\n', abs(Prob_up)/abs(Prob_dn));
fprintf('Cancellation?  |sum|/(|UP|+|DN|) = %.4f\n', ...
        abs(Prob)/(abs(Prob_up)+abs(Prob_dn)));
%% -- Physical constants -------------------------------------------------
hbar  = 4.135667696e-15 * 1.6e-19 / (2*pi);   % J*s
Bohr  = 0.529177210903;                        % 1 Bohr in Angstrom
me    = 9.1e-31;                               % kg
eps_0 = 8.85e-12;                              % F/m
c_luz = 3e8;                                   % m/s
q_e   = 1.6e-19;                               % C

%% -- Transition energy --------------------------------------------------
DE = (E_emp - E_occ) * q_e;                    % in Joules

%% -- Conversion to SI units --------------------------------------------
unit_factor = 1e10 * hbar * Bohr;
Prob_SI     = Prob * unit_factor;

%% -- Transition dipole moment ------------------------------------------
% Mu = i*hbar / ((E_f - E_i)*m) * <f|p|i>
Mu = 1i * hbar * Prob_SI / (me * DE);
fprintf('Matrix element:\n');
fprintf('  |<f|P|i>|^2  = %.4e\n', abs(Prob)^2);
fprintf('  |Mu|         = %.4e  C*m\n', abs(Mu));

%% -- Radiative rate WITHOUT the Kramers sum ----------------------------
% Gamma_R = (n_D * e^2 * DE^3) / (3*pi*eps_0*hbar^4*c^3) * |Mu|^2
% n_D = 4 for WSe2   >>> USER CHOICE: refractive-index factor of the medium
prefactor = 4 * q_e^2 * DE^3 / (3 * pi * eps_0 * hbar^4 * c_luz^3);
Gamma_R   = prefactor * abs(Mu)^2;
lifetime  = 1 / Gamma_R;

fprintf('\nResults:\n');
fprintf('  Gamma_R   = %.4e  s^-1\n', Gamma_R);
fprintf('  Lifetime  = %.4e  s\n',    lifetime);
fprintf('  log10(Gamma_R) = %.3f\n',  log10(Gamma_R));

%% -- Store results ------------------------------------------------------
results_Wvac.idx_i        = idx_occ;
results_Wvac.idx_f        = idx_emp;
results_Wvac.E_i_eV       = E_occ;
results_Wvac.E_f_eV       = E_emp;
results_Wvac.Sz_i         = sz_occ;
results_Wvac.Sz_f         = sz_emp;
results_Wvac.DE_eV        = E_emp - E_occ;
results_Wvac.Mu           = Mu;
results_Wvac.Gamma_R      = Gamma_R;
results_Wvac.lifetime_s   = lifetime;

fprintf('\nResults stored in  results_Wvac\n');

% =========================================================================
%  Inline helper (MATLAB has no native "ternary")
% =========================================================================
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
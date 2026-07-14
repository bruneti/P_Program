% =========================================================================
% Calculo_Final_SOC.m  -  VERSIÓN FINAL
% Versión adaptada para Nspin = 8 (acoplo espín-órbita, SIESTA no colineal)
%
% ─── Diferencias respecto a Calculo_Final.m (sin SOC / polarizado en espín) ─
%
%   1. siesta.Hsc es UNA SOLA matriz espinorial de dimensión (2*Nu x 2*Nu x nsc).
%      read_hodm.m la ensambla vía productos kron, sin índice de espín.
%
%   2. Los autovectores son espinores de longitud 2*Nu con orden ENTRELAZADO:
%        C(1,3,5,...) = componentes spin-up   (índices impares)
%        C(2,4,6,...) = componentes spin-down (índices pares)
%      Este orden viene de  H = kron(H_espacial, tau_espin)  en read_hodm.
%
%   3. PLMi_prueba (de Pain_Program_3il + Reordenar) tiene tamaño (2*Nu x 2*Nu),
%      pero su estructura NO es  P_espacial x I_2 . Como goDelta asigna la misma
%      integral espacial a las CUATRO combinaciones de espín del mismo par de
%      orbitales (up-up, up-dn, dn-up, dn-dn), la matriz completa es:
%
%            PLMi_prueba  =  P_espacial  x  [[1,1],[1,1]]
%
%      Usar  C' * PLMi_prueba * C  directamente incluiría términos cruzados
%      up-dn espurios que el operador momento (escalar en espín) NO debe tener.
%      La contracción espinorial correcta es:
%
%         <f|P|i>  =  C_f_up' * P_sp * C_i_up  +  C_f_dn' * P_sp * C_i_dn
%
%      donde  P_sp = PLMi_prueba(1:2:end, 1:2:end)  es el bloque up-up
%      (igual al bloque dn-dn e igual a la matriz espacial).
%
%   4. Degeneración de Kramers: con simetría de inversión temporal los estados
%      vienen en pares degenerados. La tasa radiativa DESDE un estado final debe
%      SUMAR sobre AMBOS estados iniciales degenerados:
%
%         Gamma_R(emp1) = A * ( |<occ1|P|emp1>|^2 + |<occ2|P|emp1>|^2 )
%
%      Por simetría de inversión temporal,  Gamma_R(emp2) = Gamma_R(emp1), así
%      que basta calcular UN estado final del par.
%
% ─── Flujo de trabajo ──────────────────────────────────────────────────────
%   Pain_Program_3il  ->  Reordenar  ->  este script
%   (idéntico al flujo sin SOC; solo se cambia el script final)
%
% ─── Entrada del usuario ───────────────────────────────────────────────────
%   Fija idx_occ_1, idx_occ_2, idx_emp_1, idx_emp_2 con los índices del par de
%   Kramers a evaluar. El script verifica que los pares sean realmente
%   degenerados antes de continuar.
% =========================================================================

%% -- Fase de Bloch (punto Gamma; cambia k si hace falta) -----------------
R    = siesta.lattice.supercell.R;
nsc  = siesta.lattice.supercell.Ncells;
k    = [0, 0, 0];
eikr = reshape( exp(-1i*k*R), [1 1 nsc] );

%% -- Construye H(k) y O(k) a partir de las matrices espinoriales ---------
Hk = sum( bsxfun(@times, siesta.Hsc, eikr), 3 );
Hk = (Hk + Hk') / 2;

Ok = sum( bsxfun(@times, siesta.Osc, eikr), 3 );
Ok = (Ok + Ok') / 2;

%% -- Diagonaliza (problema hermítico generalizado; autovalores ascendentes) --
[eigen_vector, eigen_vals_mat] = eig(Hk, Ok);
eigen_energies = diag(eigen_vals_mat);

%% =========================================================================
%  USER INPUT
% =========================================================================
idx_occ_1 = 647;   % >>> ELECCIÓN DEL USUARIO: pareja de Kramers 1 del par ocupado
idx_occ_2 = 648;   % >>> ELECCIÓN DEL USUARIO: pareja de Kramers 2 del par ocupado
idx_emp_1 = 649;   % >>> ELECCIÓN DEL USUARIO: estado final para el que se calcula la tasa
idx_emp_2 = 650;   % >>> (solo para comprobar la degeneración; la tasa es la misma)

%% -- Comprobación: verifica la degeneración de Kramers ------------------
tol_deg = 1e-4;   % eV — tolerancia para "degenerado"
E_occ1  = eigen_energies(idx_occ_1);
E_occ2  = eigen_energies(idx_occ_2);
E_emp1  = eigen_energies(idx_emp_1);
E_emp2  = eigen_energies(idx_emp_2);

fprintf('Par ocupado:    E(%d)=%.6f eV, E(%d)=%.6f eV, split=%.2e eV\n', ...
    idx_occ_1, E_occ1, idx_occ_2, E_occ2, abs(E_occ1-E_occ2));
fprintf('Par vacío:      E(%d)=%.6f eV, E(%d)=%.6f eV, split=%.2e eV\n', ...
    idx_emp_1, E_emp1, idx_emp_2, E_emp2, abs(E_emp1-E_emp2));

if abs(E_occ1 - E_occ2) > tol_deg
    warning('Los estados ocupados no son degenerados dentro de %.0e eV. Revisa los índices.', tol_deg);
end
if abs(E_emp1 - E_emp2) > tol_deg
    warning('Los estados vacíos no son degenerados dentro de %.0e eV. Revisa los índices.', tol_deg);
end

%% -- Extrae las componentes de espín de los autovectores espinoriales ---
% Orden entrelazado: filas impares = spin up, filas pares = spin down
idx_up = 1:2:size(eigen_vector,1);
idx_dn = 2:2:size(eigen_vector,1);

C_occ1_up = eigen_vector(idx_up, idx_occ_1);
C_occ1_dn = eigen_vector(idx_dn, idx_occ_1);

C_occ2_up = eigen_vector(idx_up, idx_occ_2);
C_occ2_dn = eigen_vector(idx_dn, idx_occ_2);

C_emp1_up = eigen_vector(idx_up, idx_emp_1);
C_emp1_dn = eigen_vector(idx_dn, idx_emp_1);

%% -- Extrae la matriz P espacial de PLMi_prueba -------------------------
% El bloque up-up es la P espacial real. Cualquier bloque de espín fuera de
% la diagonal en PLMi_prueba es un artefacto de que goDelta no distingue espín.
P_sp = PLMi_prueba(idx_up, idx_up);

%% -- Elementos de matriz espinoriales (suma sobre componentes de espín) --
% <f|P|i> = <f_up|P_sp|i_up> + <f_dn|P_sp|i_dn>
Prob_occ1_emp1 = C_occ1_up' * P_sp * C_emp1_up + C_occ1_dn' * P_sp * C_emp1_dn;
Prob_occ2_emp1 = C_occ2_up' * P_sp * C_emp1_up + C_occ2_dn' * P_sp * C_emp1_dn;

%% -- Constantes físicas (iguales que en Calculo_Final.m) ----------------
hbar  = 4.135667696e-15 * 1.6e-19 / (2*pi);   % J*s
Bohr  = 0.529177210903;                        % 1 Bohr en Angstrom
me    = 9.1e-31;                               % kg
eps_0 = 8.85e-12;                              % F/m
c_luz = 3e8;                                   % m/s
q_e   = 1.6e-19;                               % C

%% -- Diferencia de energía (usa emp1 y la media del par ocupado) --------
DE = ( E_emp1 - (E_occ1 + E_occ2)/2 ) * q_e;    % Joules
fprintf('\nGap de energía DE = %.6f eV\n', DE / q_e);

%% -- Conversión de unidades (mismo factor que en Calculo_Final.m) -------
unit_factor       = 1e10 * hbar * Bohr;
Prob_occ1_emp1_SI = Prob_occ1_emp1 * unit_factor;
Prob_occ2_emp1_SI = Prob_occ2_emp1 * unit_factor;

%% -- Momentos dipolares de transición (Cholsuk et al., JPCC 2024) -------
Mu_occ1_emp1 = 1i * hbar * Prob_occ1_emp1_SI / (me * DE);
Mu_occ2_emp1 = 1i * hbar * Prob_occ2_emp1_SI / (me * DE);

fprintf('|Mu(occ1->emp1)| = %.4e  C*m\n', abs(Mu_occ1_emp1));
fprintf('|Mu(occ2->emp1)| = %.4e  C*m\n', abs(Mu_occ2_emp1));

%% -- Tasa radiativa - MISMA FÓRMULA que Calculo_Final.m, con suma sobre --
%%    los estados iniciales degenerados ------------------------------------
% Gamma_R(emp1) = A * ( |Mu(occ1->emp1)|^2 + |Mu(occ2->emp1)|^2 )
% Gamma_R(emp2) = Gamma_R(emp1)  (simetría de Kramers)
%
% NOTA: para incluir la corrección por índice de refracción del medio,
%       MULTIPLICA la tasa por n (n ~ 1.85 para hBN, ~ 4 para WSe2). NO dividas.
%       Abajo NO se incluye, para ser fiel a la fórmula original.

prefactor = 4 * q_e^2 * DE^3 / (3 * pi * eps_0 * hbar^4 * c_luz^3);

Gamma_R  = prefactor * ( abs(Mu_occ1_emp1)^2 + abs(Mu_occ2_emp1)^2 );
lifetime = 1 / Gamma_R;

fprintf('\nGamma_R  = %.6e  s^-1\n', Gamma_R);
fprintf('Vida media = %.6e  s\n',  lifetime);

%% -- Almacena los resultados --------------------------------------------
results_SOC.idx_occ       = [idx_occ_1, idx_occ_2];
results_SOC.idx_emp       = [idx_emp_1, idx_emp_2];
results_SOC.DE_eV         = DE / q_e;
results_SOC.Mu_occ1_emp1  = Mu_occ1_emp1;
results_SOC.Mu_occ2_emp1  = Mu_occ2_emp1;
results_SOC.Gamma_R       = Gamma_R;
results_SOC.lifetime_s    = lifetime;

fprintf('\nResultados almacenados en  results_SOC\n');

% -------------------------------------------------------------------------
% Ejecuta este script tres veces (una por dirección espacial x, y, z),
% cargando el PLM_Definitivo_?.mat correspondiente y pasando antes por
% Reordenar. Compara las tres vidas medias para identificar la polarización
% de emisión, igual que en el estudio original de hBN.
% -------------------------------------------------------------------------

% =========================================================================
% Calculo_Final_W_vac.m
%
% Adaptación de Calculo_Final_SOC.m para el caso particular de W_vac:
%
%   - El defecto induce momento magnético neto.
%   - La simetría de inversión temporal está rota.
%   - NO hay degeneración de Kramers: cada estado es no degenerado.
%
% DIFERENCIAS respecto a Calculo_Final_SOC.m:
%
%   1. No se suma sobre el par de Kramers. Se calcula un único <f|P|i>
%      para la pareja de estados que el usuario especifique.
%
%   2. Se elimina la verificación de degeneración previa. En su lugar,
%      se calcula <S_z> de los estados involucrados y se verifica que
%      la transición sea spin-conservante (regla de selección óptica).
%
%   3. Se reportan explícitamente <S_z>_i, <S_z>_f, y se advierte si la
%      transición es spin-flip (suprimida) en lugar de calcular a ciegas.
%
% El resto de la lógica (extracción del bloque espacial de PLMi_prueba,
% contracción espinorial, prefactor con n=4) permanece idéntica.
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

%% -- Diagonaliza --------------------------------------------------------
[eigen_vector, eigen_vals_mat] = eig(Hk, Ok);
eigen_energies = diag(eigen_vals_mat);

%% =========================================================================
%  USER INPUT
% =========================================================================
idx_occ = 642;   % >>> ELECCIÓN DEL USUARIO: estado inicial (ocupado)
idx_emp = 643;   % >>> ELECCIÓN DEL USUARIO: estado final (desocupado)

%% -- Cálculo de <S_z> para los dos estados involucrados ----------------
N_total = size(eigen_vector, 1);
N_orbs  = N_total / 2;
Sz      = kron(eye(N_orbs), [1 0; 0 -1]/2);

sz_occ = real( eigen_vector(:,idx_occ)' * Sz * eigen_vector(:,idx_occ) );
sz_emp = real( eigen_vector(:,idx_emp)' * Sz * eigen_vector(:,idx_emp) );

E_occ = eigen_energies(idx_occ);
E_emp = eigen_energies(idx_emp);

fprintf('========================================================\n');
fprintf('  Cálculo W_vac sin pares de Kramers\n');
fprintf('========================================================\n');
fprintf('Estado inicial (i): índice %d\n', idx_occ);
fprintf('  E = %.6f eV,  <S_z> = %+.4f hbar\n', E_occ, sz_occ);
fprintf('Estado final   (f): índice %d\n', idx_emp);
fprintf('  E = %.6f eV,  <S_z> = %+.4f hbar\n', E_emp, sz_emp);
fprintf('Diferencia de energía: DE = %.6f eV\n\n', E_emp - E_occ);

%% -- Comprobación de regla de selección de spin ------------------------
threshold_sz = 0.1;   % por debajo de esto consideramos "mezcla"

if abs(sz_occ) < threshold_sz || abs(sz_emp) < threshold_sz
    fprintf('AVISO: al menos uno de los estados tiene |<S_z>| < %.2f.\n', threshold_sz);
    fprintf('       La transición es de carácter mixto. El resultado se\n');
    fprintf('       puede calcular pero su interpretación requiere cautela.\n\n');
elseif sign(sz_occ) ~= sign(sz_emp)
    fprintf('AVISO: la transición tiene carácter SPIN-FLIP.\n');
    fprintf('       <S_z>_i y <S_z>_f tienen signos opuestos.\n');
    fprintf('       El operador momento no induce flip de spin, por lo que\n');
    fprintf('       el resultado debería salir cercano a cero.\n');
    fprintf('       Procedo con el cálculo de todos modos.\n\n');
else
    fprintf('OK: transición SPIN-CONSERVANTE (ambos en canal %s).\n\n', ...
            ternary(sz_occ > 0, 'UP', 'DOWN'));
end

%% -- Extracción de las componentes espinoriales ------------------------
idx_up = 1:2:N_total;
idx_dn = 2:2:N_total;

C_occ_up = eigen_vector(idx_up, idx_occ);
C_occ_dn = eigen_vector(idx_dn, idx_occ);

C_emp_up = eigen_vector(idx_up, idx_emp);
C_emp_dn = eigen_vector(idx_dn, idx_emp);

%% -- Extracción del bloque espacial de PLMi_prueba ---------------------
% Misma lógica que en Calculo_Final_SOC.m: los cuatro bloques de spin son
% idénticos, así que basta extraer el bloque up-up.
P_sp = PLMi_prueba(idx_up, idx_up);

%% -- Elemento matricial espinorial -------------------------------------
% <f|P|i> = (C_f^up)' * P_sp * C_i^up + (C_f^dn)' * P_sp * C_i^dn
% El operador momento es escalar en spin, por lo que NO hay términos cruzados.
Prob = C_occ_up' * P_sp * C_emp_up + C_occ_dn' * P_sp * C_emp_dn;
Prob_up = C_occ_up' * P_sp * C_emp_up;
Prob_dn = C_occ_dn' * P_sp * C_emp_dn;

fprintf('\n--- Descomposición espinorial ---\n');
fprintf('Contribución UP:    %.4e + %.4ei\n', real(Prob_up), imag(Prob_up));
fprintf('Contribución DN:    %.4e + %.4ei\n', real(Prob_dn), imag(Prob_dn));
fprintf('Suma:               %.4e + %.4ei\n', real(Prob),    imag(Prob));
fprintf('|UP|/|DN|:          %.4f\n', abs(Prob_up)/abs(Prob_dn));
fprintf('Cancelación?  |suma|/(|UP|+|DN|) = %.4f\n', ...
        abs(Prob)/(abs(Prob_up)+abs(Prob_dn)));
%% -- Constantes físicas ------------------------------------------------
hbar  = 4.135667696e-15 * 1.6e-19 / (2*pi);   % J*s
Bohr  = 0.529177210903;                        % 1 Bohr en Å
me    = 9.1e-31;                               % kg
eps_0 = 8.85e-12;                              % F/m
c_luz = 3e8;                                   % m/s
q_e   = 1.6e-19;                               % C

%% -- Energía de transición --------------------------------------------
DE = (E_emp - E_occ) * q_e;                    % en Joules

%% -- Conversión a unidades SI -----------------------------------------
unit_factor = 1e10 * hbar * Bohr;
Prob_SI     = Prob * unit_factor;

%% -- Momento dipolar de transición -----------------------------------
% Mu = i*hbar / ((E_f - E_i)*m) * <f|p|i>
Mu = 1i * hbar * Prob_SI / (me * DE);
fprintf('Elemento matricial:\n');
fprintf('  |<f|P|i>|^2  = %.4e\n', abs(Prob)^2);
fprintf('  |Mu|         = %.4e  C*m\n', abs(Mu));

%% -- Tasa radiativa SIN suma sobre Kramers ----------------------------
% Gamma_R = (n_D * e^2 * DE^3) / (3*pi*eps_0*hbar^4*c^3) * |Mu|^2
% n_D = 4 para WSe2   >>> ELECCIÓN DEL USUARIO: factor de índice de refracción del medio
prefactor = 4 * q_e^2 * DE^3 / (3 * pi * eps_0 * hbar^4 * c_luz^3);
Gamma_R   = prefactor * abs(Mu)^2;
lifetime  = 1 / Gamma_R;

fprintf('\nResultados:\n');
fprintf('  Gamma_R   = %.4e  s^-1\n', Gamma_R);
fprintf('  Vida media = %.4e  s\n',    lifetime);
fprintf('  log10(Gamma_R) = %.3f\n',  log10(Gamma_R));

%% -- Almacenar resultados ---------------------------------------------
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

fprintf('\nResultados almacenados en  results_Wvac\n');

% =========================================================================
%  Helper inline (MATLAB no tiene "ternary" nativo)
% =========================================================================
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
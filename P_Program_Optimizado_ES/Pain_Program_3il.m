%% Pain_Program_3il.m
% Versión que calcula las 3 direcciones (il = -1, 0, 1) en una sola ejecución.
% Las tablas ILLL y Yrml se construyen una sola vez y se reutilizan para las 3.
% Salida: PLM_Definitivo_y.mat, PLM_Definitivo_z.mat, PLM_Definitivo_x.mat

%load('siesta.mat')
[v, k, n, R, SLZ] = goVl(siesta);
[Delta_r_complete, LM] = goDelta(siesta);
Inferno        % genera Inferno_JiLM  (ya incluye las 3 direcciones en dim 1)
Inferno_cero   % genera Inferno_zero e Inferno_extra

l1 = LM(:,1);  l2 = LM(:,2);
m1 = LM(:,3);  m2 = LM(:,4);
s1 = LM(:,5);  s2 = LM(:,6);
z1 = LM(:,7);  z2 = LM(:,8);

N         = length(Delta_r_complete);
lmax_max  = 7;

% ── Máscaras ─────────────────────────────────────────────────────────────
zero_mask   = Delta_r_complete(:,1)==0 & Delta_r_complete(:,2)==0 & Delta_r_complete(:,3)==0;
nonzero_idx = find(~zero_mask);
M_nz        = length(nonzero_idx);

% ── Resultado: 3 columnas, una por dirección (il=-1, il=0, il=1) ─────────
PLM_Delta = zeros(N, 3);

% ═════════════════════════════════════════════════════════════════════════
%  Caso Delta_r = 0
% ═════════════════════════════════════════════════════════════════════════
fprintf('Calculando elementos con Delta_r = 0...\n');
idx0 = find(zero_mask);
for i = idx0(:)'
    for il = -1:1
        PLM_Delta(i, il+2) = ...
            Inferno_zero(il+2, l1(i)+1, l2(i)+1, m1(i)+l1(i)+1, m2(i)+l2(i)+1) * ...
            Inferno_extra(l1(i)+1, l2(i)+1, s1(i), s2(i), z1(i), z2(i));
    end
end

% ═════════════════════════════════════════════════════════════════════════
%  PRECÁLCULO 1 — Tabla ILLL
%  Independiente de il -> se calcula una sola vez para las 3 direcciones
% ═════════════════════════════════════════════════════════════════════════
fprintf('Construyendo tabla ILLL...\n');
tic

% Mapa (species, l, zeta) -> índice de fila en v
slz_to_VL = containers.Map('KeyType','char','ValueType','int32');
for idx = 1:size(SLZ,2)
    key = sprintf('%d_%d_%d', SLZ(1,idx), SLZ(2,idx), SLZ(3,idx));
    if ~isKey(slz_to_VL, key)
        slz_to_VL(key) = int32(idx);
    end
end

% Índices VL para los pares no-cero
VL1_v = zeros(M_nz,1,'int32');
VL2_v = zeros(M_nz,1,'int32');
for ii = 1:M_nz
    i = nonzero_idx(ii);
    key1 = sprintf('%d_%d_%d', s1(i), l1(i), z1(i));
    key2 = sprintf('%d_%d_%d', s2(i), l2(i), z2(i));
    if isKey(slz_to_VL, key1), VL1_v(ii) = slz_to_VL(key1); end
    if isKey(slz_to_VL, key2), VL2_v(ii) = slz_to_VL(key2); end
end

% Normas únicas de Delta_r
norms_nz = vecnorm(Delta_r_complete(nonzero_idx,:), 2, 2);
[unique_norms, ~, norm_uid] = unique(norms_nz);
N_norms = length(unique_norms);

% Pares (VL1,VL2) únicos
unique_VL = unique([VL1_v, VL2_v], 'rows');
N_pairs   = size(unique_VL, 1);

VL_pair_map = containers.Map('KeyType','char','ValueType','int32');
for pp = 1:N_pairs
    key = sprintf('%d_%d', unique_VL(pp,1), unique_VL(pp,2));
    VL_pair_map(key) = int32(pp);
end
pair_uid = zeros(M_nz, 1, 'int32');
for ii = 1:M_nz
    key = sprintf('%d_%d', VL1_v(ii), VL2_v(ii));
    pair_uid(ii) = VL_pair_map(key);
end

% Tabla ILLL(par, lr+1, norm_idx)
ILLL_table = zeros(N_pairs, lmax_max+1, N_norms);
k_col = k(:);
for pp = 1:N_pairs
    vl1 = unique_VL(pp,1);
    vl2 = unique_VL(pp,2);
    if vl1==0 || vl2==0, continue; end
    base = k_col.^3 .* conj(v(vl1,:)).' .* v(vl2,:).';  % Nk x 1
    for lr = 0:lmax_max
        for nr = 1:N_norms
            jl = sphbesselj(lr, k_col * unique_norms(nr));
            ILLL_table(pp, lr+1, nr) = trapz(k, base .* jl);
        end
    end
end
fprintf('  Tabla ILLL lista en %.1f s\n', toc);

% ═════════════════════════════════════════════════════════════════════════
%  PRECÁLCULO 2 — Armónicos esféricos reales
%  Independiente de il -> se calcula una sola vez para las 3 direcciones
% ═════════════════════════════════════════════════════════════════════════
fprintf('Precalculando Yrml...\n');
tic

[unique_dr, ~, dr_uid] = unique(Delta_r_complete(nonzero_idx,:), 'rows');
N_udelta = size(unique_dr, 1);

Yrml_cache = zeros(N_udelta, lmax_max+1, 2*lmax_max+1);
for ud = 1:N_udelta
    for lr = 0:lmax_max
        for mr = -lr:lr
            Yrml_cache(ud, lr+1, mr+lmax_max+1) = ...
                Spherical_armonic_realB(unique_dr(ud,1), unique_dr(ud,2), unique_dr(ud,3), lr, mr);
        end
    end
end
fprintf('  Yrml listo en %.1f s\n', toc);

% ═════════════════════════════════════════════════════════════════════════
%  BUCLE PRINCIPAL — las 3 direcciones en cada iteración
% ═════════════════════════════════════════════════════════════════════════
fprintf('Calculando elementos PLM para il = -1, 0, +1...\n');
tic

l1_nz = l1(nonzero_idx);
l2_nz = l2(nonzero_idx);
m1_nz = m1(nonzero_idx);
m2_nz = m2(nonzero_idx);

for ii = 1:M_nz
    ll1  = l1_nz(ii);
    ll2  = l2_nz(ii);
    mm1  = m1_nz(ii);
    mm2  = m2_nz(ii);
    lmax = ll1 + ll2 + 1;
    pp   = pair_uid(ii);
    nr   = norm_uid(ii);
    ud   = dr_uid(ii);

    val = zeros(1,3);   % val(1)=il-1, val(2)=il0, val(3)=il+1

    for lr = 0:lmax
        ILLLv = ILLL_table(pp, lr+1, nr);
        if ILLLv == 0, continue; end
        for mr = -lr:lr
            Yrml = Yrml_cache(ud, lr+1, mr+lmax_max+1);
            if Yrml == 0, continue; end
            factor = (-1i)^lr * ILLLv * Yrml;
            for il = -1:1
                val(il+2) = val(il+2) + factor * ...
                    Inferno_JiLM(il+2, ll1+1, ll2+1, lr+1, mm1+ll1+1, mm2+ll2+1, mr+lr+1);
            end
        end
    end

    prefactor = (-1)^(mm1+mm2) * 4*pi * sqrt(4*pi/3);
    for il = -1:1
        PLM_Delta(nonzero_idx(ii), il+2) = (-1)^il * prefactor * val(il+2);
    end
end
fprintf('  Bucle principal listo en %.1f s\n', toc);

% ═════════════════════════════════════════════════════════════════════════
%  Guardar resultados — mismo formato que el programa original
% ═════════════════════════════════════════════════════════════════════════
PLM_Delta_y = PLM_Delta(:,1);   save('PLM_Definitivo_y.mat', 'PLM_Delta_y');
PLM_Delta_z = PLM_Delta(:,2);   save('PLM_Definitivo_z.mat', 'PLM_Delta_z');
PLM_Delta_x = PLM_Delta(:,3);   save('PLM_Definitivo_x.mat', 'PLM_Delta_x');
fprintf('Guardados PLM_Definitivo_x/y/z.mat\n');

% Nota para Reordenar.m y Calculo_Final.m:
%   Cargar el archivo de la dirección deseada y renombrar la variable:
%   load('PLM_Definitivo_z.mat'); PLM_Delta = PLM_Delta_z;
%   Luego ejecutar Reordenar y Calculo_Final como siempre.

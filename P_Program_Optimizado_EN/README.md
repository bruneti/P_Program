# P_Program_Optimizado (optimized version) — English comments

Optimized version of the radiative dipole transition-rate calculation. It gives
the same physics as the original but runs much faster: the P-matrix is built by
**precomputing** the ILLL table (unique Δr norms and unique (species,l,zeta)
pairs) and caching the spherical harmonics per unique Δr, and the three emission
directions are done in a single pass.

Reading routines extracted and trimmed from the **GROGU** package
(J. Ferrer & L. Oroszlany, 2021-2022).

## Pipeline

```
Cargar_SIESTA  ->  Pain_Program_3il  ->  Reordenar  ->  Calculo_Final[_SOC | _W_vac]
```

1. **`Cargar_SIESTA.m`** — reads the SIESTA files (`siesta`, `post`, `constants`).
2. **`Pain_Program_3il.m`** — computes the P-matrix elements for the THREE
   directions at once and saves `PLM_Definitivo_x/y/z.mat`.
3. **`Reordenar.m`** — load the desired direction, rename to `PLM_Delta`, then
   run it to get the matrix form (and the Bloch sum with more than one k-point).
4. **Final step** — pick the script that matches the calculation:
   - **`Calculo_Final.m`** — no SOC / spin-polarised (Nspin = 1 or 2).
   - **`Calculo_Final_SOC.m`** — spin-orbit (Nspin = 8), sums over the Kramers
     pair.
   - **`Calculo_Final_W_vac.m`** — broken time-reversal (no Kramers), checks the
     spin selection rule via ⟨S_z⟩.

## Configuration

Edit `SIESTA_flags.m` (input paths, NetCDF, k-point). Points where the user must
choose a value (input paths, state indices, k-point, refractive-index factor)
are marked in the code with `>>> USER CHOICE`.


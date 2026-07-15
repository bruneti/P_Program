# P_Program (original version) - English comments

Original version of the radiative dipole transition-rate calculation, the one
explained step by step in the thesis. Reads a SIESTA simulation and computes
the P-matrix (momentum operator), the transition dipole μ, the radiative rate
Γ_R and the associated lifetime.

Reading routines extracted and trimmed from the **GROGU** package
(J. Ferrer & L. Oroszlany, 2021-2022).

## Pipeline

```
Cargar_SIESTA  ->  Pain_Program  ->  Reordenar  ->  Calculo_Final
```

1. **`Cargar_SIESTA.m`** - reads the SIESTA files and leaves the `siesta`,
   `post` and `constants` structures in the workspace. Saves `siesta.mat`
   and `post.mat`.
2. **`Pain_Program.m`** - computes the P-matrix elements for ONE emission
   direction. Choose it with `il` (-1 -> y, 0 -> z, 1 -> x). Run it three times
   (once per direction) to obtain x, y and z.
3. **`Reordenar.m`** - turns the string-form output into matrix form and, with
   more than one k-point, performs the Bloch sum with `post.k`.
4. **`Calculo_Final.m`** - computes ⟨|P|⟩, μ, Γ_R and the lifetime.

## Configuration

Edit `SIESTA_flags.m`:
- `flags.input.directory` / `flags.input.file` - the SIESTA simulation to read.
- `flags.netcdf` - H/O/DM in NetCDF (`'yes'`) or old binary (`'no'`).
- `flags.kpoint.ik` - k-point for the eigen-energies (HOMO/LUMO).

Points where the user must choose a value are marked in the code with
`>>> USER CHOICE`.



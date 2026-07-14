# P_Program (versión original) — comentarios en español

Versión original del cálculo de la tasa de transición dipolar radiativa, la que
se explica paso a paso en la tesis. Lee una simulación de SIESTA y calcula la
matriz P (operador momento), el dipolo de transición μ, la tasa radiativa Γ_R y
la vida media asociada.

Rutinas de lectura extraídas y recortadas del paquete **GROGU**
(J. Ferrer & L. Oroszlany, 2021-2022).

## Pipeline

```
Cargar_SIESTA  ->  Pain_Program  ->  Reordenar  ->  Calculo_Final
```

1. **`Cargar_SIESTA.m`** — lee los ficheros de SIESTA y deja en el workspace las
   estructuras `siesta`, `post` y `constants`. Guarda `siesta.mat` y `post.mat`.
2. **`Pain_Program.m`** — calcula los elementos de la matriz P para UNA dirección
   de emisión. Se elige con `il` (−1 -> y, 0 -> z, 1 -> x). Ejecútalo tres veces
   (una por dirección) para obtener x, y y z.
3. **`Reordenar.m`** — pasa la salida en forma de cadena a forma matricial y, con
   más de un k-point, hace la suma de Bloch con `post.k`.
4. **`Calculo_Final.m`** — calcula ⟨|P|⟩, μ, Γ_R y la vida media.

## Configuración

Editar `SIESTA_flags.m`:
- `flags.input.directory` / `flags.input.file` — la simulación de SIESTA a leer.
- `flags.netcdf` — H/O/DM en NetCDF (`'yes'`) o binario antiguo (`'no'`).
- `flags.kpoint.ik` — k-point para las autoenergías (HOMO/LUMO).

Los puntos donde el usuario debe escoger un valor están marcados en el código
con `>>> ELECCIÓN DEL USUARIO`.


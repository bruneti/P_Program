# P_Program_Optimizado (versión optimizada) — comentarios en español

Versión optimizada del cálculo de la tasa de transición dipolar radiativa. Da la
misma física que la original pero es mucho más rápida: la matriz P se construye
**precalculando** la tabla ILLL (normas de Δr únicas y pares (especie,l,zeta)
únicos) y cacheando los armónicos esféricos por cada Δr único, y las tres
direcciones de emisión se hacen en una sola pasada.

Rutinas de lectura extraídas y recortadas del paquete **GROGU**
(J. Ferrer & L. Oroszlany, 2021-2022).

## Pipeline

```
Cargar_SIESTA  ->  Pain_Program_3il  ->  Reordenar  ->  Calculo_Final[_SOC | _W_vac]
```

1. **`Cargar_SIESTA.m`** — lee los ficheros de SIESTA (`siesta`, `post`, `constants`).
2. **`Pain_Program_3il.m`** — calcula los elementos de la matriz P para las TRES
   direcciones a la vez y guarda `PLM_Definitivo_x/y/z.mat`.
3. **`Reordenar.m`** — carga la dirección deseada, renómbrala a `PLM_Delta` y
   ejecútalo para obtener la forma matricial (y la suma de Bloch con más de un
   k-point).
4. **Paso final** — elige el script según el cálculo:
   - **`Calculo_Final.m`** — sin SOC / polarizado en espín (Nspin = 1 o 2).
   - **`Calculo_Final_SOC.m`** — espín-órbita (Nspin = 8), suma sobre el par de
     Kramers.
   - **`Calculo_Final_W_vac.m`** — inversión temporal rota (sin Kramers),
     comprueba la regla de selección de espín vía ⟨S_z⟩.

## Configuración

Editar `SIESTA_flags.m` (rutas de entrada, NetCDF, k-point). Los puntos donde el
usuario debe escoger un valor (rutas, índices de estados, k-point, factor de
índice de refracción) están marcados en el código con `>>> ELECCIÓN DEL USUARIO`.

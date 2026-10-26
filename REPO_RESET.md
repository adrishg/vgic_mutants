# Reinicio del repositorio — 2026-10-26

Este repositorio reemplaza al anterior, que quedó archivado en GitHub.

- **Este repo (código, activo):** `~/Repositories/vgci_mutants`
- **El anterior (completo, congelado):** `~/Repositories/vgci_mutants_archive`

## Por qué

1. **El masking approach tenía un bug.** Todos los resultados derivados del
   repositorio anterior (RMSD, RMSF, distancias, estadísticas y figuras) se
   consideran inválidos y se van a regenerar desde cero.
2. El repositorio anterior acumuló ~13 GB en Git LFS (905 objetos), en su
   mayoría versiones superadas de tablas intermedias. Mantener esos datos en
   la historia no aportaba nada y arriesgaba que alguien volviera a leer
   números inválidos.

## Qué se trajo

- Todo el código: scripts de análisis, paquetes de cluster (`scripts/*_cluster_package/`),
  configuraciones de Slurm y YAML de máscaras.
- Los 47 notebooks, **con los outputs borrados**. No hay un solo número ni
  figura del análisis anterior incrustado en ellos.
- Los `run_config.json` de cada corrida de convergencia (configuración, no resultados).
- Los registries de `docs/` y toda la documentación en Markdown, incluida la
  que describe los problemas conocidos.
- Las estructuras experimentales de referencia en `*/experimental/` (vienen del
  PDB, no las afecta el bug).

## Qué NO se trajo

Ningún dato generado: CSV, NPY, figuras, logs, resúmenes de corrida ni salidas
de MolProbity. El árbol de carpetas sí se conservó, con `.gitkeep`, para que el
pipeline escriba donde corresponde desde la primera corrida.

## Dónde está lo anterior

| Qué | Dónde |
|---|---|
| Historia completa (34 commits, todo el código) | `~/Archive/vgci_mutants_history_2026-10-26.bundle` |
| `.git` completo, incluidos los objetos LFS | `~/Archive/vgci_mutants_git_2026-10-26.tar.zst` |
| Archivos sin commitear al momento del corte | `~/Archive/vgci_mutants_untracked_2026-10-26.tar.zst` |
| Copia local completa del anterior | `~/Repositories/vgci_mutants_archive` |
| Repositorio original en GitHub | `adrishg/vgci_mutants_archive` (archivado) |

Para restaurar la historia:

```bash
GIT_LFS_SKIP_SMUDGE=1 git clone vgci_mutants_history_2026-10-26.bundle restaurado
```

Para restaurar además los datos LFS:

```bash
mkdir restaurado_full && cd restaurado_full
zstd -dc ../vgci_mutants_git_2026-10-26.tar.zst | tar -xf -
git checkout .
zstd -dc ../vgci_mutants_untracked_2026-10-26.tar.zst | tar -xf -
```

## Reglas para este repositorio

- **Nada de Git LFS.** No hay `.gitattributes` y no se va a agregar.
- Los datos viven en Hugging Face (versionados, con revision pinneado) y en el
  almacenamiento del cluster. El repositorio solo guarda el código que los produce.
- Los notebooks se commitean sin outputs.
- Cada figura del manuscrito debería salir de un script ejecutable, no de una
  corrida manual. Ver `docs/FIGURE_INDEX.md`.

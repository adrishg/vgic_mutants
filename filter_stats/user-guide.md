# User guide: inputs and settings

The pipeline runs on one dataset described by a config file:

1. Step 1 (`pipeline/build_dv.py`) superposes every ColabFold model and computes how far each lies from
   its group's other models (`spread`) and, with an experimental structure, from the structure (`offset`).
2. The topology screen (`pipeline/screen.py`, optional) removes misfolded models.
3. Step 2 (`pipeline/region_dv.py`) computes the same measures for each region named in the config.
4. The statistics (`pipeline/stats.py`) compare each condition with the control.

Paths are relative to the folder.

---

## 1. Requirements

- Python 3.11 or later, with numpy, scipy, pandas, statsmodels, gemmi and biopython
  (`requirements.txt`). On the lab machine use `~/anaconda3/bin/python3`.
- About 0.9 GB of disk and 2.5 GB of memory for a Kv2.1-sized dataset (four archives).

---

## 2. Inputs

1. **One ColabFold archive per condition**: a `.tar.gz` of one ColabFold job.
   - Every condition must reach the same highest recycle.
2. **The same construct in every condition.** The only differences allowed are at listed mutation
   sites (markers).
3. **The protein's full UniProt sequence**, as a FASTA file holding one sequence.
4. **Optionally, an experimental structure** (PDB or mmCIF, plain or `.gz`), with one chain per model
   chain, the author's chain IDs and residue numbers, and every fit residue resolved in every chain.
   Without a structure, only spread is computed.
5. **A config file**, `datasets/<dataset>/<dataset>.toml`, with the FASTA beside it.

---

## 3. The config

Paths are relative to the config file. Residue numbers are UniProt numbers; ranges are inclusive
`[start, end]` pairs.

| Key | Required | Meaning |
|---|---|---|
| `name` | yes | Dataset name; also the default output folder `results/<name>/` |
| `uniprot_fasta` | yes | The UniProt FASTA file |
| `output_dir` | no | Another output folder (not inside `data/`) |
| `model_chains` | yes | Chain IDs in the model files, e.g. `["A", "B", "C", "D"]` |
| `[reference] path` | with a structure | The structure file |
| `[reference] chains` | with a structure | Its chain IDs, as many as `model_chains` |
| `[reference] mutations` | no | Known differences from UniProt: `{ "50" = "W" }`; `X` for an unknown residue |
| `[reference] ignore_chains` | no | Structure chains to skip, e.g. an antibody |
| `[ranges] fit` | yes | The core all models are superposed on |
| `[ranges] measure` | yes | Where variation is measured |
| `[ranges] treated` | no | Residues the treatment changes; must not overlap `fit` |
| `[regions] <name>` | no | A region for step 2, as ranges inside `measure`: `VSD = [[187, 316]]`. Regions may overlap and have gaps |
| `[design] factors` | yes | Factor names, e.g. `["seq", "protocol"]` |
| `[design] control` | yes | The control condition's name |
| `[frame]` | no | Only without a structure: the model to superpose on (`model`, `seed`; optional `condition`, `recycle`) |
| `[[conditions]] name` | yes | Letters, digits, `_`, `-`, `.` |
| `[[conditions]] archive` | yes | The condition's archive |
| `[[conditions]] levels` | yes | One level per factor, as text: `{ seq = "WT", protocol = "plain" }` |
| `[[conditions]] markers` | no | Expected residue at each mutation site: `{ "416" = "F" }` |
| `[screen] domain_swapped` | with a screen | `true`: each voltage sensor sits on the next pore domain clockwise seen from outside (Kv1–9, KCNQ, Nav, Cav); `false`: on its own (KCNH, HCN, BK) |
| `[screen] wiring`, `own_place`, `box_corners`, `s6_on_own_s5` | no (default `true`) | Switch a check off with `false` (below) |
| `[screen] filter_quarter` | no (default `false`) | Switch this check on with `true` |
| `[screen] omit` | no | A CSV of predictions removed by hand: header `condition,model,seed,reason` |
| `[segments]` | with a screen | One subunit's segments, in order, inside `measure`: `S1`, `S2`, `S3`, `S4`, `linker`, `S5`, `pore_helix`, `filter`, `S6_upper`, `S6_lower`, each `[start, end]` |
| `[[segments]]` | instead, for one chain | Four such tables, one per repeat (DI–DIV of a Nav or Cav channel) |

**The screen's checks:**

- `wiring`: each voltage sensor touches most the pore domain that `domain_swapped` says it should.
- `own_place`: each helix stays nearer its own place in the structure than any other helix's (needs a
  structure).
- `box_corners`: the four lower S6 helices, seen down the pore, keep a box shape.
- `s6_on_own_s5`: each S6 has more contacts with its own subunit's S5 than with any other.
- `filter_quarter`: each selectivity filter stays nearest its own place around the pore (needs a
  structure).

**Naming rules:** config.py checks these and stops with a message if one is broken.
**Region names:** may contain letters, digits and - only (no _).
**Condition anmes, region names, factor levels must not look like a number, a true/false value, or an empty value**

**Allowed:**. WT. F412L_masked, example10
**Do Not Use:** 10, 01, True, NA, None, inf
---

## 4. Run

```bash
PYTHON=~/anaconda3/bin/python3 ./run_all.sh datasets/<dataset>/<dataset>.toml
```

This runs step 1, the screen (only if the config has `[screen]`), step 2 and the statistics, in order.
Each can also be run alone:

```bash
~/anaconda3/bin/python3 -B pipeline/build_dv.py  --config datasets/<dataset>/<dataset>.toml
~/anaconda3/bin/python3 -B pipeline/screen.py    --config datasets/<dataset>/<dataset>.toml
~/anaconda3/bin/python3 -B pipeline/region_dv.py --config datasets/<dataset>/<dataset>.toml
~/anaconda3/bin/python3 -B pipeline/stats.py     --config datasets/<dataset>/<dataset>.toml
```

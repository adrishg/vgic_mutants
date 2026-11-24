# The filter chain, and why its steps run in this order

The chain is `filter_stats/filter_chain.py`. It screens AlphaFold/ColabFold predictions for
structural plausibility before any statistics are run. This note records the order the steps run
in, what fixes that order, and what the order does *not* change — the last point matters, because
it is what keeps the chain from being an argument about which check is "most important".

The unit the chain judges is a **prediction**: one (model, seed) trajectory, its eleven recycles
taken together. Every step but the last asks a question of the whole trajectory.

## The order

| # | Step | What it asks | Level |
|---|------|--------------|-------|
| 1 | all-OK-3 convergence | from which recycle on does every later successive-recycle transition stay under 3 Å? That recycle is the trajectory's **representative**, and it and everything after it are **retained** | trajectory |
| 2 | selectivity-filter ring | does the ring open past 20 Å in any retained recycle? (Kv2.1 G377, Cav1.2 EEEE, Nav1.5 DEKA) | trajectory |
| 3a | pore–VSD interface, generalized | does the lower S6 let go of its own S4–S5 linker (past 12 Å) in any retained recycle? | trajectory |
| 3b | pore–VSD interface, published | does S6's cytoplasmic end leave the residues before S1 (past 27 Å) in any retained recycle? **Kv2.1 only** | trajectory |
| 4 | topology screen | does the representative pass the reference-free checks `wiring`, `box_corners`, `s6_on_own_s5`? | model |

## What fixes the order

**Convergence has to be first, because the later steps are defined on its output.** It is also the
one step whose own definition is fragile: it walks a *contiguous* recycle series, comparing each
recycle with the previous available one. Any step that removes recycles out of the middle of a
trajectory changes what "the previous available recycle" means. Running such a step first does not
reorder two filters — it silently redefines one of them, and "converged" stops meaning what the
published number means.

**The topology screen has to be last, because it is the only per-model step.** Two reasons, both
measured (26-11/filterComparison/`FILTER_ORDER_FINDINGS.md`, section 9):

- The trajectory-level checks exist to catch a trajectory that goes wrong in *one* recycle. A
  per-model screen run first deletes exactly those recycles, and the trajectory-level checks then
  find almost nothing: 8 Kv2.1 L403A trajectories rejected instead of 107, and 11 WT instead of
  297. The ensemble is not cleaner; the QC has been blinded.
- The headline metric then rises in the WT control as well as in the mutant. That is
  permissiveness, not signal.

The screen judges the representative only, not all eleven recycles. That is about eleven times
cheaper, but the reason is that the representative is the unit the rest of the analysis reports.

**Recycle 0 is judged.** The convergence step frequently selects recycle 0 as the representative —
for over half of Kv2.1 L403A's trajectories. A topology screen that starts at recycle 1 therefore
does not clean the ensemble; it lets exactly those representatives through unjudged. Both our chain
and the statistics configs (`[screen] first_recycle = 0`) judge it.

## What the order does **not** change

Steps 2, 3a and 3b are **absolute thresholds applied per model and then taken to the trajectory**.
None of them computes a statistic over the set another step has already thinned: whether a given
trajectory fails one of them does not depend on which other trajectories are still present. So
their order among themselves changes **no** survivor and **no** count — it changes only which name
a doomed trajectory is filed under.

That is still worth fixing, because the per-prediction verdict table reports one reason per
prediction. The convention is: **3a before 3b**, the generalized check before the Kv2.1-only
published one, so that one structural failure — the gate off its linker — carries one name in all
three channels rather than appearing as a Kv2.1-specific reason in Kv2.1 and a different reason
elsewhere.

On Kv2.1 the two overlap completely: 3a is a **strict subset** of 3b. Running both removes exactly
what 3b removes alone (38+259 = 297 WT, 19+88 = 107 L403A, 24+199 = 223 F412L), which is only
possible if nested. 3a earns its place for Nav1.5 and Cav1.2, where 3b does not exist, not for
Kv2.1.

The topology verdicts are also order-free in practice: the expensive `judge` step is run once over
*all* convergence representatives, before any of steps 2–3b, and the chain then reads the table.
Changing a cutoff is therefore a question asked of a CSV, not a reason to walk the models again.

## What each step is worth, measured

- **The selectivity-filter ring has never rejected anything** — 0 trajectories in all 110 rows of
  the setup sweep, with the widest ring 10.8 Å in Cav1.2, 12.0 Å in Nav1.5 and 18.4 Å in Kv2.1
  F412L against a 20 Å cutoff. Keep it as a tripwire and report it as one; do not report it as a
  filter that cleaned anything.
- **The pore–VSD interface is where the discrimination is**, in Kv2.1. It is also the step with no
  published equivalent outside Kv2.1, which is what 3a was calibrated for: fourteen experimental
  structures (six Kv2.1, two Cav1.2, six Nav1.5), closest S6-to-own-linker approach 4.8–9.2 Å, the
  subunit's own linker always nearest, cutoff set at 12 Å.
- **Nav1.5 and Cav1.2 models sit inside the experimental envelope**: widest S6-to-linker 8.24 Å
  over roughly 42,000 models, against 9.17 Å in the experimental set. The chain is not finding
  these ensembles clean by failing to look.

## Open question, recorded rather than resolved

Kv2.1 WT has a much looser gate–linker joint than either mutant (median 8.72 Å against 5.46 and
5.51 Å, 95th percentile above the 12 Å cutoff). Sequence or masking artifact is unresolved, and it
is a confound for any WT-versus-mutant comparison that leans on step 3a.

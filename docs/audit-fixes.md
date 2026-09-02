# Plan: closing audit findings A1, A2, A4 and A6

## Context

The audit found six class-A findings — defects that change scientific results.
Two are now closed: A5 (out-of-box penalty measured from the box centre) and A3
(swap-pool draws that never terminate), in `7ea27e1` and `c6dcf24`.

The remaining four are **not** fixed in either upstream branch. `master` and
`cyclic` are byte-identical in all four regions; the hunks near them are
comments and whitespace. So unlike A5, there is no upstream implementation to
port — each needs a decision and a fix written here.

They are also not the same kind of problem, and lumping them together is how
this goes wrong. Two are unambiguous defects with one obviously correct answer.
Two are force-field changes that no amount of code review can settle, because
the parameters were fitted against the current behaviour.

## The four, sized against what we can actually measure

| | what | targets that exercise it | can the benchmark adjudicate it? |
|---|---|---|---|
| **A6** | `ramabias` indexes `ramaprob[]` with no range check | measured: **0 out-of-range in 200k steps** | n/a — latent |
| **A2** | `hasCYS` tests only `'C'`, so methionine's Sδ scores against `C.map` | **6 of 49** (M without C) | barely |
| **A4** | `lower_gridenergy`'s log compression is applied to `e.map` and `d.map` too | every target with a charged residue | no — see below |
| **A1** | `hydrophobic_low` returns a ramp that **rises** with distance | **48 of 50** | no — see below |

### A6 — latent, guard it and move on

`segphi`/`segpsi` come from `atan2`, whose range is `(-π, π]`. At φ = π exactly,
`segphi` reaches 180 and `ind` reaches 32400, one past the end of the array
`ramaprob_initialise()` allocates. At ψ = π the index silently wraps into the
next φ row and returns the wrong bin.

Instrumented and run for 200,000 steps on 1SFI: **zero** out-of-range indices.
It needs geometry idealised enough to hit φ = π exactly, which the MC does not
produce. So this is a latent read, not an active one.

**Fix**: clamp both segment indices to `[0, 179]`. Three lines.
**Validation**: must be *observationally inert* — the clamp only fires where
the old code read out of bounds, which measurably never happens. Prove it the
way A3 was proven: byte-identical trajectories across the benchmark.

### A2 — a real defect, but the benchmark barely covers it

`AD_init()` loads `SA.map` into the sulfur slot only when the peptide contains
cysteine. Methionine's Sδ is also atom type 4 (`canonicalAA.cpp`), so a peptide
with Met and no Cys silently scores its sulfur against the **carbon** map.

Proven by deleting `rigidReceptor.SA.map`: a Met peptide runs and produces a
**bit-identical** trajectory, so the map was never read; the Cys control dies
with `Missing gridmap_file.map file`. Magnitude, comparing the two maps over the
3Q47 pocket outside the steric wall: median −0.003, **p95 +3.04 kcal/mol per
atom**, and the best pocket point is −1.075 (SA) against −0.838 (C) — 28% of the
sulfur's best affinity lost.

**Fix**: `if (id == 'C' || id == 'M') hasCYS = 1;` — and rename the variable,
since it no longer means "has cysteine". One line plus a rename.

**Coverage problem**: only **6 of 49** manifest targets have Met without Cys, and
none of them is in a set with a per-target published reference that isolates
sulfur. The fix is right on inspection, but the benchmark cannot demonstrate it
is an improvement. Either accept it on the argument (a sulfur scored on a carbon
map is wrong regardless of outcome) or add Met-rich targets first.

### A4 — needs a decision about what "correct" means

`lower_gridenergy()` compresses map values above +2.718 to `ln(E) + 1.718`. It
is applied inside `gridmap_initialise()`, so it hits all nine maps — including
`e.map` (electrostatic potential) and `d.map` (desolvation), loaded at
`main.cpp:916-917`.

Because it only fires above +2.718, the **positive** lobe of the Coulomb field
is squashed while the **negative** lobe is untouched. The field is no longer
antisymmetric in the sign of the charge: a positive ligand atom near a positive
receptor patch is penalised far less than a negative atom near a negative patch
is. That is a physics error, not a tuning choice.

But "don't compress e/d" is not obviously the whole answer:

- the compression exists to soften the receptor's steric wall, which is a real
  need for a Monte Carlo search that must be able to climb out of clashes
- `e.map` has no steric wall — its large values are genuine electrostatics
- `d.map` is desolvation, always ≥ 0 in AutoDock4, so compression there is a
  different question again

**Fix**: give `gridmap_initialise()` a flag, and pass it false for slots 7 and 8.
Small. The decision is whether the desolvation map should also be exempt, and
that wants a measurement, not an opinion.

**Validation gap**: no target in the manifest was chosen to stress electrostatics.
Adjudicating this needs targets with charged binding sites, which is a manifest
question before it is a code question.

### A1 — the one that is not a bug fix

```c
if (distance > contact_cutoff + range) return 0.0;
if (distance < contact_cutoff)         return 1.0;
return (distance - contact_cutoff) / range;   /* rises 0 -> 1 */
```

Its sibling `linear_decay()`, 300 lines up with the same signature and doc-comment
style, returns `1.0 - (distance - cutoff) / width`. Both commented-out
alternatives in the same file (`hydrophobic_low_recip`, `hydrophobic_low_spline`)
also decay. Three of four candidate forms decay; this one rises. The header
comment above it describes `f(d) = d^-1`, a decaying form.

As written the well is non-monotonic with two discontinuities: a pair at contact
weighs 1, just past contact ~0, then climbs back to 1 at the far edge of the
2.8 Å range and drops to 0 beyond it.

**Measured**: rebuilding with `1.0 -` restored and re-running the whole
validation reshuffles the ranking completely (top-5 becomes 15, 14, 8, 3, 12
instead of 15, 12, 13, 4, 10), makes energies systematically more negative
(−29.06 vs −28.55 RT), and moves the top-1 RMSD from **0.66 Å to 0.84 Å**.

That last number is why this is not a simple fix. `kauzmann_param = 0.122` is a
contrastive-divergence fit, and the code does not record which functional form
it was fitted against. Correcting the form without refitting the weight trades a
known bias for an unknown one — and on the one system we can check, it makes the
published result *worse*.

**This is a recalibration, not a patch.** Doing it properly means refitting
`kauzmann_param` against a validation set, which is a research task, not a code
change.

## Order

1. **A6** — clamp. Inert by construction, provable, no decision needed.
2. **A2** — one condition. Right on the argument; benchmark coverage is thin and
   the commit message should say so rather than claim a measured improvement.
3. **A4** — exempt `e.map`; decide separately about `d.map`. Land the flag, and
   record that the manifest cannot yet adjudicate it.
4. **A1** — **do not patch.** Open it as a recalibration item with the measured
   evidence attached. If it is ever changed, `kauzmann_param` moves with it.

## Validation

Per fix, using `tests/validation/`:

- `ctest` (22 tests) and `func_adcp_fold_determinism` 12×, per MIGRATION.md's
  standing rule for `energy`/`main` changes
- `val_3q47_redock` must still reach 0.66 Å, **except** where a fix is expected
  to move it, in which case the new value is recorded
- A6 gate: byte-identical trajectories over the 49-target smoke tier. Anything
  else means the clamp is firing where it should not
- A2 gate: the 43 targets without Met-without-Cys must be byte-identical; the 6
  affected targets are reported as a delta, not as an improvement
- A4 gate: results will move; record the per-target diff in `docs/compares/` style

## What this plan does not solve

The recurring problem is the same for A2 and A4: **the manifest was assembled to
reproduce published tables, not to isolate individual force-field terms.** Six
Met targets and no charged-site selection cannot adjudicate a sulfur-typing bug
or an electrostatic-map bug.

Closing A2 and A4 with confidence needs targets chosen for those terms — a
different kind of set from A/B/C, and one no paper publishes. That is a
prerequisite worth naming before promising these can be "validated".

# Plan: closing audit findings A1, A2, A4 and A6

## Context

The audit found six class-A findings — defects that change scientific results.
Four are now closed: A5 (out-of-box penalty measured from the box centre), A3
(swap-pool draws that never terminate), A6 (unchecked `ramaprob[]` index) and A2
(methionine's sulfur scored on the carbon map), in `7ea27e1`, `c6dcf24`,
`797eb72` and `7e38b98`.

The two that remain — A1 and A4 — are both **force-field recalibrations, not
patches**, and both were written, measured and reverted rather than argued
about. Neither is closed and neither should be patched in isolation.

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
| **A4** | `lower_gridenergy`'s log compression is applied to `e.map` and `d.map` too | measured: **32 of 50**, 0.0097% of map points | yes, and it did — the fix loses |
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

### A4 — implemented, measured, and reverted. Not a patch either.

`lower_gridenergy()` compresses map values above +2.718 to `ln(E) + 1.718`. It
is applied inside `gridmap_initialise()`, so it hits all nine maps — including
`e.map` (electrostatic potential) and `d.map` (desolvation), loaded at
`main.cpp:950-951`.

Because it only fires above +2.718, the **positive** lobe of the Coulomb field
is squashed while the **negative** lobe is untouched. The field is no longer
antisymmetric in the sign of the charge: a positive ligand atom near a positive
receptor patch is penalised far less than a negative atom near a negative patch
is. That reasoning is unchanged and still correct.

**The fix was written and it does not survive measurement.** A `soften` flag on
`gridmap_initialise()`, passed false for slots 7 and 8, builds clean and passes
all 22 tests. Then the 3Q47 redock, four independent blocks of 16 seeds each:

| seeds | 1-16 | 101-116 | 201-216 | 301-316 |
|---|---|---|---|---|
| compressed (HEAD) | 0.70 | 0.84 | 0.78 | 0.74 |
| `e.map` exempt | 0.85 | **2.58** | **9.47** | 0.90 |

The search is not what breaks. Ranking every pose in each ensemble shows the
near-native pose is still **found** every time — it just stops **winning**:

| | top-1 | best in ensemble | its rank |
|---|---|---|---|
| compressed, seeds 201-216 | 0.78 | 0.78 | 1 of 16 |
| exempt, seeds 201-216 | 9.47 | 0.94 | 2 of 16 |
| exempt, seeds 101-116 | 2.58 | 1.30 | 3 of 16 |

A 0.94 A pose loses to a 9.47 A pose under the uncompressed map. Two of four
seed blocks fail, against zero of four before.

**Why a 0.01% change in the map does this.** Measured over all 50 prepared
targets: 1,467 of 15,192,578 `e.map` points are above the threshold (0.0097%),
in 32 of the 50 targets, largest delta 1.043 kcal/mol/e — about 0.36 kcal/mol
once multiplied by ADCP's backbone charges. That is small in absolute terms and
comparable to the gaps that separate adjacent poses in the ranking, which is a
0.25/0.75 blend of total and external energy. The compression was in place when
those weights and `kauzmann_param` were fitted.

So A4 lands in exactly the same category as A1: **the physics argument is right
and the recalibration it implies has not been done.** Removing the compression
without refitting the ranking weights trades a known bias for an unknown one,
and on the only system with a crystallographic answer it measurably loses.

`d.map` is a separate and simpler story: **0 of 15,192,578 points across all 50
targets exceed the threshold.** The compression provably never fires on the
desolvation map, so exempting it is a literal no-op and there is nothing to
decide.

**Status**: not patched. Reopened as a recalibration item alongside A1. If the
ranking weights are ever refitted, both move together — they are the same
experiment, not two.

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
3. **A4** — **do not patch.** Implemented and reverted: the exemption breaks the
   3Q47 ranking in 2 of 4 seed blocks. Recalibration item, paired with A1.
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
- A4 gate: **failed**. Four seed blocks of the 3Q47 redock, top-1 backbone RMSD
  0.70/0.84/0.78/0.74 compressed against 0.85/2.58/9.47/0.90 exempt

## Revisiting A1 and A4: how the refit would actually work

Parked, not abandoned. Recording what was established so it does not have to be
re-derived.

### The refit does not require re-docking

Both failures are **ranking** failures — the near-native pose is found every
time and loses. That decouples the fit from the search, and three properties of
the code make the fit arithmetic over a table rather than a compute campaign:

1. `kauzmann_param` is a **pure linear prefactor** on the hydrophobic term
   (`energy.cpp:956`, `return -k_h * energy * intensity`). So
   `E(k_h) = E_rest + k_h * H`, and any `k_h` is a reweighting, not a re-run.
2. A1 changes the shape of `hydrophobic_low`, but `H` follows from the pose
   geometry. That is two numbers per pose (`H_rise`, `H_decay`), not two runs.
3. A4 is which map file was loaded. Also two numbers per pose.

Dump roughly six components per pose **once**
(`E_rest, H_rise, H_decay, Eelec_compressed, Eelec_raw, totalE, extE`) and the
parameter sweep is seconds.

### The fitting rig is already in the tree

`energy.cpp:2916` `energy_probe_1()` is CRANKITE's contrastive-divergence
gradient machinery, intact: 36 parameter slots, `energy_probe_1_last`/`_this`
holding the CD data and model expectations, and `energy_probe_1_calc[]`
selecting what to fit. **Slot 9 is already `kauzmann_param`**, slot 10 is
`hydrophobic_cutoff_range`, slot 11 the dielectric. This is not a from-scratch
project; it is the rig that produced the 0.122.

`adcp -f pose.pdb -r 1x0` was tried and does load a structure and print
`totalE`/`extE` in 0 seconds. In that probe the pose was judged out of the box
(extE 699993, the out-of-box penalty) because the AutoSite translation points
were not staged alongside it. Wiring those up is setup, not architecture.

### Three things are missing, and only one is code

1. **A training set disjoint from A/B/C.** The 49 targets are the test set;
   fitting and validating on them is how a benchmark stops meaning anything.
   Use **peptiDB (set F) minus set E** — the plan already marks it `ref = NA`
   precisely because no ADCP number is published for it, so it contaminates no
   comparison.
2. **An objective.** CD maximises the likelihood of native *conformations* — a
   generative fit, and the one that produced these values. What broke here is
   **discrimination**: the right pose exists and loses. Fit against the *rank*
   of the lowest-RMSD pose, which is literally the metric that failed.
3. **The component dump.** The only code: a mode that reads a pose and prints
   the component vector above.

Fit A1 and A4 **jointly** — two continuous knobs (`kauzmann_param`,
`hydrophobic_cutoff_range`) and two discrete ones (hydrophobic form, `e.map`
compression), so four grid scans over two parameters. The 0.25/0.75 ranking
blend is a fifth knob if the four are not enough.

### The trap

The published ADCP numbers were produced **with** A1's inverted form. So
"recovers the paper's numbers" is not the same as "correct" — some of the
published performance may depend on the bias.

The honest tiebreaker is crystallographic RMSD, not published fnc. A
parameterisation that ranks the near-native pose first more often **and** matches
the paper settles it. One that improves RMSD while worsening published fnc is a
finding about the paper, not about this fork.

### Cost, and the cheap first move

| step | work | compute |
|---|---|---|
| feasibility probe on existing ensembles | ~1 day | minutes |
| real fit against a training set | +2 days | ~1 night of docking |
| validation: re-dock A/B/C | — | ~40 min (smoke tier) |

The feasibility probe pays for itself alone: if **no** combination of
(`k_h`, form, compression) ranks better than the current one even on the test
set, the recalibration idea dies there. Only if something wins is it worth
assembling the training set properly.

Note the ensembles it would run against (49 targets x 8 replicas x 30 models,
across five run sets) lived in a scratch directory and are not preserved.
Regenerating them is the ~40 min smoke tier, so this is not a blocker.

**First concrete step**: stage the translation points into the rescore probe and
dump the component vector for the 3Q47 ensemble. Everything else depends on it.

## What this plan does not solve

**The manifest was assembled to reproduce published tables, not to isolate
individual force-field terms.** Six Met targets and no charged-site selection
cannot adjudicate a sulfur-typing bug, so A2 landed on the argument rather than
on a measurement, and the commit says so.

A4 is the sharper lesson. The benchmark *did* adjudicate it — but only
negatively: it showed the fix loses on the one target with a crystallographic
answer. It cannot show what a refitted `kauzmann_param` and ranking blend would
do, because refitting needs a training set, and the 49 targets here are the
test set. Using them for both is how a benchmark stops meaning anything.

So A1 and A4 are blocked on the same missing thing, and it is not more code:
**a separate fitting set, disjoint from A/B/C, with charged binding sites.**
Naming that is more useful than promising these can be "validated" as they
stand.

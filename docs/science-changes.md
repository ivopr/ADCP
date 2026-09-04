# Science changes in this fork

Every change in this repository that **moves a scientific result** — an energy,
a trajectory, a pose, a ranking — with why it was necessary and how it was
done. Build fixes, memory hygiene, leaks and tooling are out of scope unless
they changed a number; those live in [MIGRATION.md](../MIGRATION.md).

Companion documents: [docs/audit-fixes.md](audit-fixes.md) (the class-A audit
and what was deliberately left alone), [docs/cyclic-port.md](cyclic-port.md)
(what upstream's `cyclic` branch has and why most of it is not ported), and
[docs/compares/5.md](compares/5.md) (the four-way measurement this table's
"status" column is checked against).

## Summary

| # | change | commit | kind | result moved? |
|---|---|---|---|---|
| 1 | `sqrt` on a float expression lost precision | `55b58b6` | port regression | yes — restored pristine values |
| 2 | denormal `currTargetEnergy` hung the fold path | `68dd3ae` | port regression | no — removed a hang |
| 3 | `swapChains` off-by-one crashed every production run | `7e2192d` | upstream defect | yes — runs complete at all |
| 4 | `transmutate`/`transmove`/`transopt` dropped `chi1`/`chi2` | `5dd2a23` | upstream defect | yes — side-chain dihedrals no longer garbage |
| 5 | swap-pool draws that can never terminate (A3) | `c6dcf24` | upstream defect | no — 375/375 healthy runs byte-identical |
| 6 | nested sampling seeded from the discarded point | `35fb3fb` | upstream defect | yes — NS sequence changed by design |
| 7 | out-of-box penalty, cysteine scoring, gamma centroid (A5) | `7ea27e1` | force field | yes — all 16 probe seeds |
| 8 | Ramachandran table index clamp (A6) | `797eb72` | latent read | no — inert by construction |
| 9 | sulfur map loaded for methionine, not only cysteine (A2) | `7e38b98` | force field | yes — 6 of 49 targets |
| 10 | `MaxRotamers`, opt-in rotamer cap | `58337c5` | opt-in knob | only when enabled |
| 11 | optimizing-strategy weights in `model_params` | `9cad378` | plumbing | no |
| 12 | `cdlearn` never set per-protein sequence | `ae0bfc6` | tool science | yes — CD learning ran at all |
| 13 | `cdlearn` missing `ramaprob_initialise()` | `3bf4003` | tool science | yes |
| 14 | three bugs across `probe.cpp`'s 32 diagnostic bits | `6db4341` | diagnostics | no |
| 15 | disulfide model: soft chi3, two-pair cap, `0.25*SSloss` | `dff5df3` | force field | yes — any peptide with 2+ cysteines |
| 16 | CA-CA harmonic on adjacent pairs + repaired ring closure | `240afae` | force field | yes — every peptide; fold path now matches `cyclic` bit for bit |
| 17 | `transmutate` returns `int`, FIXED guards kept | `0b1fce1` | signature | no |
| 18 | no-improvement stop raised 10M → 30M steps | `74185b4` | search budget | only on runs that would have stopped early |

Two further class-A findings — **A1** (inverted hydrophobic ramp) and **A4**
(log compression on `e.map`) — were implemented, measured, and **reverted**.
They are recalibrations, not patches. See §6.

---

## 1. Restoring what the C→C++ port broke

### `sqrt` on a float expression lost precision — `55b58b6`

**Necessity.** Docking diverged from pre-migration ADCP on every run past
200,000 steps. Not a cosmetic difference: the Monte Carlo search amplifies it
into an entirely different trajectory.

**Mechanism.** Five sites in `energy.cpp` of the form

```c
float v[3];
n = 1. / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
```

In C, `sqrt` is `double sqrt(double)` — the only one there is — so the float
sum widens and the root is taken in double precision. In C++ the overload set
is visible, the argument binds to `float sqrt(float)`, and the root is rounded
to single precision *before* the divide. One extra rounding, 2⁻²³ ≈ 1.2e-7
relative.

**How it was isolated.** After the two obvious suspects were tested and
cleared: reverting the VLA→`std::vector` change alone still diverged;
reverting the `View3` flattening still diverged; rebuilding both sides at `-O0`
gave the same values, ruling out codegen; tracing intermediates at `-O0` showed
`v1/v2/v3` differing on call #1, before any coordinate is read. Casting the
five `sqrt` arguments to double, changing nothing else, produced an exact
pristine match.

**Validation.** Bit-level equality with pristine across the 16-seed sweep — a
property still held today and re-measured in every comparison since.

### Denormal `currTargetEnergy` hung the fold path — `68dd3ae`

**Necessity.** `func_adcp_fold_determinism` hung roughly 1 run in 3 — same
binary, same seed, same command.

**Mechanism.** `currTargetEnergy` picked up a denormal (~1e-308, a different
value each run), which makes the retry loop at `main.cpp:306` spin forever,
since every `swapEnergy[i] >= currTargetEnergy` comparison succeeds.

**How.** Bisected by running the identical fold 12× per build: pre-migration C
0/12 hangs, `adf422e` 0/12, `5660a33` 3/12, `90c63a8` 4/12, and `90c63a8` with
only `energy` reverted to C 0/12 — which localised it to the ported energy
code rather than to the loop itself.

---

## 2. Defects inherited from upstream `master`

These were present, byte-identical, in pristine `1c1a330`. Comparison 5
established that upstream's own `cyclic` branch is free of the crash and the
seed-10 divergence — so these are `master`-lineage defects, and our fixes
closed a gap `cyclic` had already closed by another route (except A3, which
both upstream branches still carry).

### `swapChains` off-by-one crashed every production run — `7e2192d`

**Necessity.** ADCP segfaulted on its own documented production settings. The
wrapper's defaults (`runADCP.py: -N 50 -n 2500000`) crashed **100% of the
time**. Worse, the crashed runs were plausible: the process still wrote a
50-model PDB with sensible energies, and only the missing `best target energy`
line and the exit status revealed the failure.

**Mechanism.** `main.c:213` declared the pool one element short:

```c
double swapEnergy[swapLength + 1];   /* 11 elements, correct */
Chain*  swapChains[swapLength];      /* 10 elements, off by one */
for (int i = 0; i < swapLength + 1; i++)   /* writes [10] */
```

The intent is unambiguous — `swapEnergy` is sized `+1`, the comment says "last
element is with the best energy", and index `swapLength` is used directly at
three sites. Only the declaration was short.

**How the threshold was pinned.** Every use of the corrupted slot is gated
behind `swapMutateSteps = 200000`, so shorter runs never reach it. That made
the boundary exactly reproducible: 20k and 50k steps exit 0 before and after;
200k, 1M and 2.5M steps SIGSEGV before and complete after.

**Standing evidence.** Reproduced again in comparison 5 under
`-fsanitize=address,undefined` in earlier sessions (`main.c:237`, index 10 out
of bounds) and, this session, by the crash map: pristine loses 12/16, 14/16,
14/16 seeds at 50k/150k/250k steps; HEAD, `cyclic` and the shipped binary lose
none at any budget.

### `transmutate`/`transmove`/`transopt` never copied `chi1`/`chi2` — `5dd2a23`

**Necessity.** MSan caught `chain->aa[i].chi1` read uninitialized in
`crankshaft()` on **every docking run using `external=5` + `Opt=1`** — that is,
essentially every production docking configuration. Present byte-identical in
pristine C, 12 commits before the migration.

**Mechanism.** The three `transX()` functions each copy a hand-listed subset of
AA fields into a scratch `chaint->aat[j]`, omitting `chi1`/`chi2` (a
commented-out full-struct copy sits next to the incomplete list in all three),
then commit the move with a full struct assignment — silently overwriting
valid side-chain dihedrals with whatever sat in the scratch slot.
`crankshaft()` is innocent; it is only where the corruption is detected.

**How.** Add the two omitted fields to each copy loop, plus `aat_init()` as
defence in depth. Five lines, no refactor. `rotate_cyclic()` was checked and
ruled out (it rebuilds `chi1`/`chi2` via `acidate()`); `flex.cpp`/`nested.cpp`
have the same shape but are backfilled before first use and were deliberately
left alone rather than bundling an untested NS/MPI change into the fix.

### Swap-pool draws that can never terminate — A3, `c6dcf24`

**Necessity.** `simulate()`'s stagnation branch draws a pool member *strictly
better* than the current chain:

```c
while (swapEnergy[swapInd] >= currTargetEnergy) swapInd = rand() % 11;
```

The branch is reached **because** the search stagnated, and stagnation is
exactly the state in which the current chain is already the pool's best — so
the condition can be unsatisfiable and the draw never returns. Two sibling
loops have the same defect.

**Evidence.** Reproduced on the production cyclic path: 1SFI seeds 2 and 3,
3P8F seed 7, 4KEL seeds 6 and 8 spin at 100% CPU with no output and no
timeout — measured to 1h53m. Also reachable from the documented
`-A AMPLITUDE,FIX_AMP` flag. The trigger is the search *state*, not the seed,
which is why pristine and HEAD hang on different seeds from identical code,
and why any instrumentation made it "disappear".

**How.** Bounded to 10 draws per slot; with nothing better in the pool the
chain stays put and annealing handles the stagnation.

**Validation.** 375 of 375 healthy replica trajectories stay byte-identical,
and the 17 that previously produced nothing now complete. **Both upstream
branches still carry this loop unchanged** — this fix is not redundant with
`cyclic`.

### Nested sampling seeded new points from the point being discarded — `35fb3fb`

**Necessity.** Nested sampling's validity rests on drawing from the
**surviving** prior mass. The serial live-point draw was

```c
do  copies = 1 + ((int)(rand()/(double)RAND_MAX * N)) % N;   /* 1..N */
while (copies == 1 && N > 1);
```

but `find_worst()` has already swapped the heap root out to
`chainhash[heaplength]` and decremented, so with `P == 1` the discarded worst
point sits at index `N`. Drawing 1..N therefore re-seeded from the discarded
point with probability 1/N, biasing the evidence estimate.

**How it was confirmed** — by instrumentation, not inference: a run drew
`copies == 20` with `N == 20`, and that entry's `ll` equalled `logLstar`
exactly, which is by definition the discarded point's likelihood. The
`PARALLEL` branch a few hundred lines up already gets this right with
`% (N-P)`; the fix is the same expression at `P == 1`.

**Consequence, by design.** The NS energy sequence changed and stays changed —
`40.058829 31.225769 10.764044 7.463047 6.393763` before,
`40.058829 688.128472 81.995422 10.222738 7.086147` after. Comparison 5 shows
pristine, `cyclic` and the shipped ADFRsuite binary all still produce the old
sequence; this fork is alone in producing the corrected one.

---

## 3. Force-field corrections

### Out-of-box penalty, cysteine scoring, gamma centroid — A5, `7ea27e1`

Three changes taken from upstream's `cyclic` branch, each independently
flagged by this repo's own audit. **All three change docking results by
construction** — this is the commit that moved HEAD off pristine's frozen
table (comparison 5 attributes 100% of the movement to it, by per-commit
bisect).

1. **`gridenergy()` out-of-box penalty.** It measured distance from the box
   **centre** — `(exactGrid - N/2)² / 20` — which put a ~38 RT step at the
   wall, failed to point back into the box, used integer division for `N/2`,
   and let an atom flung past the 1e9 escape pay only 10, *less than one that
   barely left*. It now measures from the nearest face, damped and capped.
   **Beyond upstream**: the index computation is moved *below* the out-of-box
   return, because `getindex()` takes ints and converting an out-of-range
   double to int is undefined behaviour. Past that return every axis satisfies
   `0 < exactGrid < N-1`, so all eight stencil indices are in range by
   construction.
2. **`ADenergyNoClash()` cysteine.** Cysteine was the only residue scored by a
   lone `gridenergy()` call on the pseudo-gamma position — no rotamer library,
   no clash check, no SH hydrogen. It now goes through `scoreSideChainNoClash`
   like every other side chain.
3. **`scoreSideChain{,NoClash}()` gamma position.** The gamma pseudo-atom is
   the side chain's **centroid**, not its first rotamer atom. `aadict.cpp`
   calibrates every gamma radius — hydrophobic contact, vdW, H-bond — against
   a pseudo-atom near the charged tip (LYS 4.700 Å from CB, ARG 4.900).
   Writing `tc[SCRot][0]` put it at CG, ~1.5 Å out, so those radii were being
   applied in the wrong place. `bestSideChainCenter` was already computed and
   thrown away.

### Sulfur map loaded for methionine, not only cysteine — A2, `7e38b98`

**Necessity.** `AD_init()` loaded `rigidReceptor.SA.map` into grid slot 4 only
when the peptide contained cysteine. Slot 4 is the AutoDock **sulfur** map, and
it is needed by every residue whose rotamer library uses atom type 4 — CYS
*and* MET (`canonicalAA.cpp`: CYS `{4,3}`, MET `{4,0,0}`). A peptide with
methionine and no cysteine silently scored its Sδ against the **carbon** map.

**Proof it was never read.** Delete `SA.map` from the working directory:
before the fix a Met peptide runs and produces a **bit-identical** trajectory
(so the map was never loaded), while the Cys control dies with
`Missing gridmap_file.map file`. After the fix the Met peptide dies too, and
its energy moves (-7.538 vs -7.646, same seed and step count).

**Magnitude.** Comparing the two maps over the 3Q47 pocket outside the steric
wall: median −0.003, **p95 +3.04 kcal/mol per atom**, best pocket point −1.075
(SA) against −0.838 (C) — a methionine sulfur was losing **28%** of its best
affinity.

**How.** `hasCYS` becomes `hasS` and tests `'C' || 'M'`. One condition plus a
rename.

**Honest limits.** The benchmark cannot show this is an *improvement*: only 6
of 49 manifest targets have Met without Cys, and none sits in a set whose
published reference isolates sulfur. What it can confirm is the mechanism, and
the prediction — exactly 1CM1, 2F31, 3AVI, 3AVJ, 5UWI and 6CIT move; the three
targets with both Met and Cys stay byte-identical because they already loaded
the map; the other 40 are untouched — was written down **before** the run and
came back exact: 43 byte-identical, 6 changed.

### Ramachandran table index clamp — A6, `797eb72`

**Necessity.** `ramaprob`/`alaprob`/`glyprob` are 180×180 tables.
`dihedral_rama()` returns `atan2`'s range, (−π, π], so φ = π exactly yields
`segphi == 180` — one past the end of the allocation — while ψ = π yields
`segpsi == 180`, which does not overrun but silently wraps into the next φ row
and scores the **wrong bin**.

**How the scope was established.** Instrumented and run for 200,000 steps on
1SFI: **zero** out-of-range indices. Reaching φ = π needs geometry more
idealised than the MC produces, so this closes a *latent* read, not an active
one — which is why the fix is a clamp and not a redesign.

**Validation.** Required to be observationally inert, and was: across the
49-target smoke tier every target this could touch stays byte-identical. The 6
that move do so because of the methionine fix in the next commit, and which 6
was predicted before the run.

---

## 4. Opt-in knobs, off by default

### `MaxRotamers` — `58337c5`

Upstream's `cyclic` branch hardcodes a cap of 20 rotamers for residues that
have more (LYS/ARG 81, GLN 36, MET/GLU 27), drawn at random. That is where its
3–4× speedup on cyclic peptides comes from — the only hunk in the whole branch
diff that is neither a comment nor a science change.

**It is a speed-for-accuracy trade, not a free optimisation**: the subset is
redrawn on every call, so the energy becomes an even noisier function of the
coordinates. Measured over the 49-target benchmark it buys **2.25× wall clock**
(4825 s → 2140 s) and recovers some average fnc, but **halves the correlation**
with the published per-target values (0.27 → 0.17).

**Default is 0** — try every rotamer, which is what this code has always done.
`-p MaxRotamers=20` reproduces upstream, including its oddities (it overwrites
the loop counter with a fresh draw, so it samples *with replacement*, and tests
the counter after drawing) — reproducing that ordering is the point, since
otherwise the flag would consume `rand()` differently and not reproduce
upstream at all. `run_set.sh` records the setting in every `run_info.txt` so no
report can silently mix capped and uncapped runs.

### Optimizing-strategy weights in `model_params` — `9cad378`

Adds `opt`, `opt_totE_weight`, `opt_firstlastE_weight`, `opt_extE_weight`,
default-initialised. Plumbing for the `Opt=` option set; no result moves.

---

## 5. Tool-level science fixes

- **`cdlearn` never set per-protein sequence — `ae0bfc6`.** The per-protein
  setup loop cloned the shared `sim_params` into every slot, but nothing had
  ever set that struct's `seq`/`sequence`/`NAA` from an actually-read protein
  (`main.cpp` always calls `update_sim_params_from_chain()`; `cdlearn.cpp`
  never did). Every real CD-learning iteration therefore either died on the
  first `move()` call or — with a stale value — silently used the **wrong
  protein's sequence** for every slot but one. This is the rig that fits the
  force-field parameters, so it matters beyond the tool.
- **`cdlearn` missing `ramaprob_initialise()` — `3bf4003`.** Same file; the
  Ramachandran tables were never initialised on that path.
- **`probe.cpp` diagnostic bits — `6db4341`.** `tests()` walks a 32-entry
  mask-gated table; the default masks exercise only 5 of 32 bits, so 27 had
  never run under a sanitizer. Sweeping all of them found three bugs
  (`number_of_contacts` leaks, a `CA_geometry` heap-overflow read, one more);
  29 bits were clean. Diagnostics only — no production number moves.

---

## 6. Deliberately not changed

Two class-A findings were written, measured and **reverted**. Both are force-field
**recalibrations**: the physics argument is right, and the parameters were
fitted against the current behaviour, so changing the form without refitting
trades a known bias for an unknown one.

### A1 — the hydrophobic ramp rises instead of decaying

```c
if (distance > contact_cutoff + range) return 0.0;
if (distance < contact_cutoff)         return 1.0;
return (distance - contact_cutoff) / range;   /* rises 0 -> 1 */
```

Its sibling `linear_decay()` — same signature, same doc-comment style — returns
`1.0 - (distance - cutoff)/width`; both commented-out alternatives in the file
also decay; the header comment describes `f(d) = d⁻¹`, a decaying form. As
written the well is non-monotonic with two discontinuities.

**Measured**: restoring `1.0 -` reshuffles the ranking completely, makes
energies systematically more negative (−29.06 vs −28.55 RT), and moves the
top-1 RMSD from **0.66 Å to 0.84 Å**. `kauzmann_param = 0.122` is a
contrastive-divergence fit and the code does not record which functional form
it was fitted against. **Not patched.**

### A4 — log compression applied to `e.map` and `d.map`

`lower_gridenergy()` compresses map values above +2.718 to `ln(E) + 1.718`,
inside `gridmap_initialise()`, so it hits all nine maps — including the
electrostatic and desolvation maps. Because it only fires above +2.718, the
**positive** lobe of the Coulomb field is squashed while the negative lobe is
untouched: the field is no longer antisymmetric in the sign of the charge.

**The fix was written and does not survive measurement.** A `soften` flag,
false for slots 7 and 8, builds clean and passes all 22 tests. Then the 3Q47
redock, four independent blocks of 16 seeds:

| seeds | 1-16 | 101-116 | 201-216 | 301-316 |
|---|---|---|---|---|
| compressed (HEAD) | 0.70 | 0.84 | 0.78 | 0.74 |
| `e.map` exempt | 0.85 | **2.58** | **9.47** | 0.90 |

The search is not what breaks — the near-native pose is still *found* every
time, it just stops *winning* (a 0.94 Å pose loses to a 9.47 Å pose). Only
0.0097% of `e.map` points are affected (1,467 of 15,192,578 across 50 targets),
but that is the same order as the gaps separating adjacent poses in a 0.25/0.75
ranking blend that was fitted **with** the compression in place. `d.map` is
simpler: **0 of 15,192,578 points** exceed the threshold, so exempting it is a
literal no-op. **Not patched.**

Both are parked as one joint recalibration item, with the fitting rig
(`energy_probe_1()`, slot 9 = `kauzmann_param`) and the missing pieces — a
training set disjoint from A/B/C, a discrimination objective, a component
dump — written up in [audit-fixes.md](audit-fixes.md).

---

## 7. What this fork does **not** take from `cyclic` / ADFRsuite

Comparison 5 measured HEAD against the shipped ADFRsuite 1.0 binary and against
upstream's `cyclic` branch, which is that binary's lineage (identical banner,
byte-identical fold trajectory). **HEAD is not a superset of either.**

| not ported | why |
|---|---|
| `cyclic_energy()` rework | it **collapses the macrocycle**: median CA1–CAn 2.52 Å across 432 poses on 18 backbone-cyclic targets, against 3.83 Å crystallographic and 3.99 Å at HEAD. Two alpha carbons at 2.5 Å is not a conformation — C–C vdW contact is ~3.4 Å. Cause: `NCDistance` is never assigned (the line is commented out) and `distance()` returns the **square**, so the harmonic centred on 3.819 Å almost never runs; what executes is a harmonic centred on **zero** |
| `energy2`'s CA–CA harmonic on every adjacent pair | same defect family; likely source of the remaining 4.5× cyclic speed gap. Buying speed by collapsing the ring is not a trade worth making |
| `sbond_energy` single-shortest-pair disulfide search, `0.25 * SSloss` | one design decision, undecided — it changes which disulfides are scored at all |
| `transmutate`'s deleted early-return on FIXED residues | upstream deleted it as a MacOS compile fix; the deletion is a behaviour change |
| `MaxRotamers = 20` **as the default** | shipped by upstream, opt-in here (§4) |

**Measured remainder.** At 10,000 steps, seed sweep 1–16: HEAD agrees with
`cyclic` bit-for-bit on **10/16** seeds, matches energy with a different pose on
2 more, and differs on 4 (seeds 2, 8, 15, 16). Against the shipped binary,
HEAD differs on 16/16 — the shipped binary is a few *internal* commits past the
public branch tip, which no public commit reproduces.

**So: HEAD does not contain everything the ADFRsuite binary does**, and some of
what it lacks it lacks on purpose. What it does have that the shipped binary
does not: no crashes on `master`'s production settings, corrected nested
sampling, the A2/A3/A5/A6 fixes, an uncollapsed macrocycle, and a better
redock on the one target with a crystallographic answer (0.70 Å vs 0.77 Å).

---

## 8. How each change was validated

The standing gates, applied per change:

1. `ctest --test-dir build` — 22 tests, must stay green.
2. `func_adcp_fold_determinism` **12×** for any `energy`/`metropolis`/`main`
   change (MIGRATION.md's standing rule; it is how `68dd3ae`'s 1-in-3 hang was
   caught).
3. `val_3q47_redock` must still pass, and its RMSD is recorded whenever a
   change is expected to move it (0.66 Å at compares/4 → 0.70 Å now;
   `7ea27e1` is the only commit in between that moves any result, per the
   per-commit bisect in [compares/5.md](compares/5.md), though the redock
   itself was not re-bisected).
4. **Inert-by-design changes** must prove it: byte-identical trajectories over
   the 49-target smoke tier (`797eb72`, `c6dcf24`) or over the 16-seed sweep.
5. **Result-moving changes** must state *which* targets move, ideally as a
   prediction made before the run (`7e38b98`: 43 identical / 6 changed, called
   in advance).
6. A full four-way comparison against pristine, `cyclic` and the shipped binary
   — [compares/5.md](compares/5.md), with the reproduction protocol in
   [5-protocol.md](compares/5-protocol.md).

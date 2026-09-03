# Behavioural changes in this fork: a technical rationale

This document is the long form of [science-changes.md](science-changes.md).
Where that document lists *what* changed, this one argues *why each change was
necessary*, in the terms the change belongs to — floating-point analysis,
Markov-chain sampling theory, the statistical mechanics of nested sampling, and
force-field calibration — and gives the concrete evidence and worked examples
behind each claim.

Every quantitative statement below is a measurement made in this repository and
reproducible from [compares/5-protocol.md](compares/5-protocol.md) or the
commit it is attributed to. Where a statement is an inference rather than a
measurement, it is marked as such.

**Contents**

1. [Epistemic standard](#1-epistemic-standard-what-counts-as-evidence-here)
2. [Instruments](#2-instruments)
3. [Class I — numerical fidelity](#3-class-i--numerical-fidelity-of-the-language-port)
4. [Class II — correctness of the sampler](#4-class-ii--correctness-of-the-sampler)
5. [Class III — the energy function](#5-class-iii--the-energy-function)
6. [Class IV — what was measured and deliberately reverted](#6-class-iv--what-was-measured-and-deliberately-reverted)
7. [Class V — what is not adopted from upstream `cyclic`](#7-class-v--what-is-not-adopted-from-upstream-cyclic)
8. [Discussion](#8-discussion)

---

## 1. Epistemic standard: what counts as evidence here

Docking codes are difficult to validate because almost any change produces
*different* numbers, and different is not the same as *wrong* or *better*. The
standard applied throughout this fork is therefore threefold, and each change
below is filed under exactly one of the three:

**(i) Bit-level invariance.** A change claimed to be a correctness fix with no
scientific consequence must produce byte-identical trajectories on a fixed
fixture set. This is a falsifiable, binary criterion and it is used wherever
the claim is "inert" (§4.3, §5.3). It is far stronger than "the tests still
pass": a Monte Carlo trajectory is a chaotic object, so any perturbation, no
matter how small, eventually shows up in the output (§3.1 quantifies how
quickly).

**(ii) Mechanistic proof.** A change claimed to fix a defect must come with a
demonstration that the defect *is real and is reached* — an ASan/MSan/UBSan
report at a named `file:line`, an instrumented counter, or an experiment that
distinguishes the two hypotheses. Reading the code and finding it suspicious is
a reason to investigate, never a reason to commit. §5.2 shows the pattern: the
claim "this map is never read" was tested by *deleting the map file* and
confirming the trajectory did not change.

**(iii) Predict, then measure.** Where a change is expected to move results,
the prediction of *which* systems move is recorded before the run. `7e38b98`
predicted exactly six of forty-nine targets by name, plus three that must stay
identical for a specific reason, and the run returned exactly that partition.
A prediction that survives contact with the data is evidence; a rationalisation
written afterwards is not.

A corollary, applied twice in §6: **a physics argument that a term is wrong is
not sufficient warrant to change it**, if the surrounding parameters were fitted
against the wrong term. That situation calls for a recalibration, which is a
research task, not a patch.

---

## 2. Instruments

| instrument | what it can decide |
|---|---|
| `func_adcp_fold_determinism`, run 12× | reveals nondeterminism and intermittent hangs that a single run hides (§3.2) |
| 16-seed docking sweep at 10,000 steps, `targetE` + ATOM md5 | bit-level equality of trajectories between two binaries |
| crash-rate map (16 seeds × {50k, 150k, 250k, 2.5M} steps) | prevalence of a failure, not merely its existence |
| 49-target benchmark (`tests/validation/`) | per-target deltas against published references |
| `val_3q47_redock` | the one target in the set with a crystallographic answer |
| ASan / UBSan / MSan | mechanistic proof of memory and initialisation defects |
| per-commit bisect over Release builds | attribution of an observed change to one commit |
| four-way comparison (pristine, HEAD, upstream `cyclic`, shipped ADFRsuite binary) | whether a defect is ours, upstream's, or lineage-specific |

The four-way comparison ([compares/5.md](compares/5.md)) deserves emphasis
because it changes the *interpretation* of several findings in §4. Upstream
publishes two branches: `master` (which this fork ports) and `cyclic` (which
`master` never merged). Comparison 5 established behaviourally — identical
`PepTide 2.0` banner, byte-identical fold trajectory, ATOM md5
`abf915bc58b0eaedf4e4d3a1353f723b` — that the binary distributed in ADFRsuite
1.0 belongs to the `cyclic` lineage, not to `master`. Consequently, defects
found in `master` are not necessarily defects in "ADCP as published", and the
distinction is made explicitly wherever it matters.

---

## 3. Class I — numerical fidelity of the language port

### 3.1 Silent precision loss under C++ overload resolution (`55b58b6`)

**The defect.** Five sites in `energy.cpp` compute an inverse norm:

```c
float v[3];
n = 1. / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
```

In C there is exactly one `sqrt`, `double sqrt(double)`; the `float` argument
undergoes the usual arithmetic conversions and the root is evaluated in binary64.
In C++ `<cmath>` introduces the overload set `{float, double, long double}`, the
argument binds to `float sqrt(float)`, and the root is *computed and rounded in
binary32* before the reciprocal is taken in binary64. The source text is
unchanged; the semantics are not. This is the class of change that a
line-by-line review of a C→C++ port is least likely to catch, because there is
nothing on the line to look at.

**Magnitude.** The extra rounding is bounded by half a unit in the last place of
binary32, i.e. a relative error of at most 2⁻²⁴ ≈ 6·10⁻⁸ (machine epsilon
2⁻²³ ≈ 1.19·10⁻⁷). In isolation this is far below any physically meaningful
threshold — it is roughly 10⁻⁷ Å on a bond vector.

**Why an error of 10⁻⁷ is nonetheless fatal.** The quantity feeds a Metropolis
acceptance test. Writing the test as `accept ⟺ ΔE < θ` for a stochastic
threshold θ, a perturbation ε to ΔE changes the decision only when ΔE lies
within ε of θ. If ΔE is distributed with density of order 1 RT⁻¹ near the
threshold, the per-step probability of a flipped decision is O(ε) ≈ 10⁻⁷, so
the expected number of flips in an n-step run is O(10⁻⁷·n): negligible at
n = 10⁴, order unity at n = 10⁷ — and a single flipped acceptance re-seeds the
whole subsequent trajectory, because the chain is deterministic given its state
and its random stream. This is precisely the observed signature: **runs up to
~200,000 steps matched pristine; longer runs diverged.** The defect was
invisible in exactly the regime where testing is cheap.

**How it was isolated** (an example of standard (ii)). Two plausible suspects —
the VLA→`std::vector` conversion and the `View3` parameter flattening — were
each reverted in isolation and each failed to restore agreement. Rebuilding both
sides at `-O0` reproduced the same divergence, excluding optimiser-driven
contraction (FMA, reassociation) as the cause. Tracing intermediates at `-O0`
showed the normal vectors differing on the *first* call, before any coordinate
had been read back. Casting the five arguments to `double` and changing nothing
else reproduced pristine bit-for-bit. The fix as it stands today is visible at
[src/energy.cpp:1501](../src/energy.cpp#L1501) and four sibling sites.

**Why it was necessary.** Without it, this fork would silently disagree with
published ADCP on every production-length run while passing every short test —
the worst possible failure mode for a scientific code, since the disagreement is
invisible in validation and material in use.

### 3.2 A denormal as a non-terminating loop condition (`68dd3ae`)

`func_adcp_fold_determinism` hung on roughly one run in three — same binary,
same seed, same command line. Nondeterminism of that kind cannot come from the
Monte Carlo stream, which is seeded; it must come from reading uninitialised
memory. `currTargetEnergy` acquired a denormal value (~10⁻³⁰⁸, different each
run), and the retry loop at `main.cpp:306` compares every pool entry
`swapEnergy[i] >= currTargetEnergy`. Against a denormal that predicate is
satisfied by every entry, so the loop never exits.

Two methodological points. First, a 1-in-3 intermittency is invisible to a test
suite that runs each test once; the standing rule of running the determinism
fixture **twelve times** after any `energy`/`metropolis`/`main` change exists
because of this bug. Second, the bisect isolated it to the ported energy code by
building `90c63a8` with *only* `energy` reverted to C (0/12 hangs against
4/12) — a differential experiment rather than an inspection.

---

## 4. Class II — correctness of the sampler

The changes in this section do not alter the energy function. They restore the
sampler to the algorithm it is documented to implement.

### 4.1 A one-element array underflow that made production settings unusable (`7e2192d`)

**The defect.** In `main.c` the swap pool was declared as

```c
double swapEnergy[swapLength + 1];   /* 11 elements */
Chain*  swapChains[swapLength];      /* 10 elements — one short */
for (int i = 0; i < swapLength + 1; i++)  /* writes swapChains[10] */
```

Every other piece of evidence in the file says `swapChains` was meant to have
`swapLength + 1` elements: `swapEnergy` is sized `+1`, the comment above the
initialisation loop states that the last element holds the best energy, and
index `swapLength` is dereferenced directly at three further sites. Only the
declaration was short. UBSan reports it as
`main.c:237:14: index 10 out of bounds for type 'Chain *[*]'`.

**Why it mattered.** Writing one pointer past a stack array is undefined
behaviour, and in practice it corrupted an adjacent stack slot. All uses of the
overflowed element are gated behind `swapMutateSteps = 200000`, so the
consequence appears only in long runs — which is to say, in *every* run the
distributed wrapper performs: `runADCP.py`'s own defaults are `-N 50 -n 2500000`,
and those crashed **100% of the time**.

**Why it was not noticed upstream.** The failure is not silent, but it is
*plausible*: the process still wrote a fifty-model PDB with sensible-looking
energies before dying. Only the absence of the final `best target energy` line
and a non-zero exit status distinguish a crashed run from a complete one. A
user reading the output file cannot tell.

**Prevalence, not anecdote.** The crash-rate map turns this from "it crashed
for me" into a measured rate, reproduced identically in five independent
sessions: pristine loses 12/16, 14/16 and 14/16 seeds at 50k, 150k and 250k
steps, and 15/16 at the full 2.5M-step protocol. This fork, upstream `cyclic`,
and the shipped ADFRsuite binary lose none at any budget — which locates the
defect in the `master` lineage specifically.

### 4.2 Proposal states that discarded side-chain dihedrals (`5dd2a23`)

**The defect.** `transmutate()`, `transmove()` and `transopt()` build a proposal
by copying a hand-enumerated subset of fields from the current residue into a
scratch structure `chaint->aat[j]`, and then commit the accepted proposal with a
whole-structure assignment back into `chain->aa[i]`. The enumerated subset omits
`chi1` and `chi2` — the side-chain dihedral angles — and `aat_init()` leaves the
same two fields uninitialised for the same reason. (A commented-out
whole-structure copy sits next to the incomplete list in all three functions,
which suggests the omission was an optimisation that lost two fields.)

The asymmetry is the whole bug: an *incomplete* copy in, a *complete* copy out.
Accepting a move therefore overwrites valid dihedrals with whatever occupied the
scratch slot.

**Why it survived so long.** MSan reports the read in `crankshaft()`
([`metropolis.cpp:727`](../src/metropolis.cpp#L727); pristine `metropolis.c:732`),
which is *not* where the corruption happens — `crankshaft()` always rewrites
those fields before reading them back. The detection site and the defect site
are in different functions, so a reviewer following the report reads innocent
code. Attribution required tracing writers rather than readers.

**Scope.** All three writers are gated on `external_potential_type == 5`, i.e.
they run only when docking against grid maps with `Opt=1` — essentially every
production docking configuration, and none of the folding tests. The defect is
byte-identical in pristine C, twelve commits before the migration: it is
upstream's, not the port's.

**Restraint.** `flex.cpp` and `nested.cpp` contain the same copy pattern but are
backfilled by `copybetween()`/`mpi_rec_chain()` before first use. They were left
untouched rather than bundling a change to the nested-sampling and MPI paths —
which have essentially no test coverage — into a docking fix.

### 4.3 Rejection loops whose acceptance predicate can be unsatisfiable (audit A3, `c6dcf24`)

**The defect.** The stagnation branch of `simulate()` draws a replacement
conformation from the swap pool by rejection:

```c
while (swapEnergy[swapInd] >= currTargetEnergy) swapInd = rand() % 11;
```

The loop is a rejection sampler over an 11-element pool whose acceptance set is
`{i : E_i < E_current}`. Rejection sampling terminates almost surely **iff the
acceptance set is non-empty**, and here it can be empty by construction: this
branch is entered *because the search stagnated*, and the state in which the
search stagnates is typically the one in which the current conformation is
already the best member of the pool. The predicate then holds for every index
and the loop spins forever.

**Evidence** (standard (ii)). Reproduced on the production cyclic path: 1SFI
seeds 2 and 3, 3P8F seed 7, 4KEL seeds 6 and 8 spin at 100% CPU with no output
and no timeout, measured to 1 h 53 min. It is also reachable from the documented
`-A AMPLITUDE,FIX_AMP` flag, which makes `move()` unreachable and therefore
stagnates immediately.

**Why it presents as a heisenbug.** The trigger is the *state* of the search,
not the seed. Pristine and this fork hang on different seeds while running
identical loop code, because their trajectories differ for unrelated reasons
(§5.1); and adding any instrumentation perturbs the trajectory off the offending
state, so the bug "disappears when observed". This is why an operational
reproduction across five targets was required before touching it.

**The fix and its termination argument.** The draw is bounded to ten attempts
per slot ([src/main.cpp:325-332](../src/main.cpp#L325-L332)); if nothing better
is found the chain stays put and annealing handles the stagnation, which is the
behaviour the branch was trying to express. Note that two sibling draws
([:429](../src/main.cpp#L429), [:531](../src/main.cpp#L531)) are deliberately
left unbounded: their predicate is `E_i > E_best + 5 kcal`, and the invariant
that slot `swapLength` holds `E_best` is re-established at every update site
([:405](../src/main.cpp#L405) and [:414](../src/main.cpp#L414)), so the
acceptance set always contains that slot and the number of draws is geometric
with p = 1/11. *(That invariant is verified by inspection of the update sites,
not by instrumentation — it is the weaker of the arguments in this document.)*

**Inertness.** 375 of 375 healthy replica trajectories remain byte-identical,
and the 17 replicas that previously produced nothing now complete. Both upstream
branches still carry the unbounded loop, so this fix is *not* redundant with
`cyclic`.

### 4.4 Nested sampling seeded from the discarded point (`35fb3fb`)

This is the change with the clearest theoretical stake, so it is worth restating
the algorithm before the defect.

**Background.** Nested sampling maintains N "live" points drawn from the prior
π constrained to a likelihood contour, L(θ) > L\*. At each iteration the worst
live point (the one attaining L\*) is removed and recorded, L\* is raised to its
likelihood, and a **new point is drawn from π restricted to the new, strictly
smaller constrained region**. The evidence estimator
Z ≈ Σᵢ Lᵢ (Xᵢ₋₁ − Xᵢ) depends on the prior-mass sequence Xᵢ, whose statistics
are known *only because* each replacement is an independent draw from the
surviving constrained prior: the shrinkage ratio tᵢ = Xᵢ/Xᵢ₋₁ is then
Beta(N,1)-distributed, with E[ln t] = −1/N. Every error bar the method reports
rests on that independence and on the replacement coming from the **surviving**
mass.

**The defect.** In the serial branch the replacement was drawn as

```c
do  copies = 1 + ((int)(rand()/(double)RAND_MAX * N)) % N;   /* 1..N */
while (copies == 1 && N > 1);
```

but `find_worst()` has already swapped the heap root out to
`chainhash[heaplength]` and decremented the count, so with P = 1 the **discarded
worst point sits at index N** while the live set occupies 1..N−1. Drawing
uniformly over 1..N therefore re-seeded the new sample from the point being
discarded with probability 1/N.

**Why that is not a rounding detail.** The discarded point sits exactly *on* the
new constraint boundary, L = L\*, which is the one point excluded from the region
being sampled. Seeding from it (i) violates the support of the target
distribution and (ii) correlates the new live point with the one just removed,
which is precisely the independence the tᵢ statistics assume. The evidence
estimate is biased, and the size of the bias is not something the algorithm's
own error bars can report, since they are derived under the assumption that was
violated.

**Confirmation by instrumentation, not inference.** A run drew `copies == 20`
with `N == 20`, and that entry's log-likelihood equalled `logLstar` exactly —
`logLstar` being, by definition, the discarded point's likelihood. The
`PARALLEL` branch a few hundred lines above already had the correct expression,
`% (N-P)`; the fix is that same expression at P = 1
([src/nested.cpp:987](../src/nested.cpp#L987)). The `while (copies == 1)` guard
was separately shown to be inert (removing it alone left the output
byte-identical) and to be excluding a legal draw — the worst *surviving* point.

**Consequence, and why the changed numbers are the correct ones.** The NS energy
sequence moved and stays moved:

```
before  40.058829  31.225769  10.764044   7.463047  6.393763
after   40.058829 688.128472  81.995422  10.222738  7.086147
```

Individual iterations going "backwards" is normal for nested sampling — the
sequence is a record of a shrinking constrained prior, not a monotone
optimisation trace. Comparison 5 shows pristine, upstream `cyclic` and the
shipped ADFRsuite binary all still emit the *old* sequence, so this fork is
currently alone in producing the corrected sampler.

---

## 5. Class III — the energy function

Changes in this section move docking results by construction. Each therefore
carries a per-target account of what moved, and none is claimed as an
"improvement" without a crystallographic reference to support the claim.

### 5.1 The out-of-box penalty (audit A5, `7ea27e1`)

**The original form.** For an atom whose fractional grid coordinate along an
axis is g with N grid points, the penalty was

    P(g) = (g − N/2)² / 20      (with N/2 in integer arithmetic)

plus an escape returning 10 for values beyond ~10⁹. Three defects follow
directly from the formula:

1. **Discontinuity at the wall.** The penalty is measured from the box *centre*,
   so an atom just outside a 55-point box pays (55/2)²/20 ≈ 38 RT while an atom
   just inside pays whatever the interpolated map says — typically O(1) RT. The
   energy surface has a step of tens of RT at the boundary. Metropolis dynamics
   on a surface with a step that large simply never crosses back.
2. **The gradient does not point home.** Because the penalty grows with distance
   from the centre rather than from the nearest face, its gradient near a face
   is dominated by the centre direction, not by the inward normal. A restraint
   whose gradient does not restore the coordinate is not a restraint.
3. **Non-monotonicity.** An atom flung past the 10⁹ escape paid 10 — *less* than
   an atom that barely left the box. The scoring function was therefore
   maximised, in that region, by throwing the ligand as far away as possible.

**The fix** ([src/energy.cpp:1905-1945](../src/energy.cpp#L1905-L1945)) measures
the distance past the **nearest face**, squares it, and damps and caps the
result (`erg > 10000 ? 10000 : erg/100`). The penalty is then continuous at the
wall, monotone in the excursion, and its gradient is the inward normal. This is
upstream's own fix, taken from the `cyclic` branch, where it had been
independently arrived at.

**One deliberate divergence from upstream.** The grid index computation is moved
*below* the out-of-box return. `getindex()` takes `int` arguments, and
converting an out-of-range `double` to `int` is undefined behaviour in C and
C++ — reachable exactly for the far-flung atoms this branch exists to handle.
Past that return, every axis satisfies 0 < g < N−1, so all eight trilinear
stencil indices are in range by construction rather than by a bounds check.
Upstream keeps a stale guard instead; this fork removes the need for one.

### 5.2 Sulfur typed as carbon for methionine (audit A2, `7e38b98`)

**The defect.** AutoDock scores each ligand atom against a precomputed affinity
map for its atom *type*. `AD_init()` loaded `rigidReceptor.SA.map` into grid
slot 4 — the sulfur map — only when the peptide sequence contained cysteine.
But the rotamer library assigns atom type 4 to two residues, not one
(`canonicalAA.cpp`: CYS `{4,3}`, MET `{4,0,0}`). A peptide containing methionine
and no cysteine therefore scored its Sδ against `rigidReceptor.C.map`, the
**carbon** map, silently and with no diagnostic.

**Proof the map was never read** — an experiment that distinguishes "loaded but
unused" from "never loaded": delete `SA.map` from the working directory. Before
the fix, a methionine peptide runs to completion and produces a *bit-identical*
trajectory (the file was never opened), while a cysteine control dies with
`Missing gridmap_file.map file`. After the fix, the methionine peptide dies too,
and its energy moves (−7.538 vs −7.646 at the same seed and step count). This is
the cleanest available demonstration that the map participates in scoring.

**Magnitude.** Comparing the two maps point-by-point over the 3Q47 pocket,
outside the steric wall: median difference −0.003 kcal/mol, 95th percentile
**+3.04 kcal/mol per atom**; the best pocket point is −1.075 (SA) against −0.838
(C). A methionine sulfur was losing about **28%** of its best available
affinity — a term comparable in size to the energy gaps that separate competing
poses.

**Why the benchmark cannot certify it as an improvement, and why it landed
anyway.** Only 6 of the 49 manifest targets contain methionine without cysteine,
and none sits in a set whose published reference isolates a sulfur contribution.
The commit therefore claims *correctness of typing*, not improved performance:
scoring a sulfur against a carbon map is wrong independently of whether it helps
on this particular benchmark. What the benchmark does establish is the
mechanism, under standard (iii): the prediction — exactly 1CM1, 2F31, 3AVI,
3AVJ, 5UWI and 6CIT move; the three targets containing *both* Met and Cys stay
byte-identical because they already loaded the map; the remaining 40 are
untouched — was recorded before the run and returned exactly: 43 identical, 6
changed.

### 5.3 Domain of the Ramachandran lookup (audit A6, `797eb72`)

`ramaprob`, `alaprob` and `glyprob` are 180×180 tables indexed by
`segphi = ⌊(φ + π)/(2·π/180)⌋` and likewise for ψ. The dihedral is computed with
`atan2`, whose range is the half-open interval (−π, π]. At φ = π exactly,
`segphi` evaluates to 180 and the flattened index reaches 32400 — one element
past the allocation. At ψ = π, the index does not overrun but *wraps into the
next φ row*, silently scoring the conformation against the wrong bin, which is
the more insidious of the two failures.

**Scope, measured before fixing.** Instrumented and run for 200,000 steps on
1SFI: **zero** out-of-range indices. Reaching φ = π exactly requires geometry
more idealised than the Monte Carlo produces. This is therefore a *latent* read,
and the appropriate response is a three-line clamp, not a redesign — the
instrumentation is what justifies the small fix rather than a large one.

**Validation under standard (i).** The change is required to be observationally
inert and is: across the 49-target smoke tier every target it could touch stays
byte-identical. The six targets that do move in the same run move because of
§5.2, and which six was predicted in advance.

### 5.4 Two calibration mismatches fixed with A5 (`7ea27e1`)

**Cysteine bypassed the side-chain machinery.** In `ADenergyNoClash()`, cysteine
was the only residue scored by a bare `gridenergy()` call on the pseudo-gamma
position — no rotamer library, no clash check, no SH hydrogen — while every
other side chain went through `scoreSideChainNoClash()`. There is no physical
argument for the exception; it now takes the same path as the rest.

**The gamma pseudo-atom was placed on the wrong atom.** ADCP represents each
side chain by a pseudo-atom whose radii — hydrophobic contact, van der Waals,
hydrogen-bond — are calibrated in `aadict.cpp` against the side-chain
**centroid**, near the charged tip: 4.700 Å from CB for lysine, 4.900 Å for
arginine. The code wrote `tc[SCRot][0]`, the first rotamer atom (CG), roughly
1.5 Å short of that. The parameters were therefore being applied at a position
they were not fitted for — a calibration mismatch rather than an outright bug,
and the sort that produces plausible-looking wrong answers. `bestSideChainCenter`
was already being computed and discarded.

**Attribution.** Comparison 5's per-commit bisect shows that these three changes
together (§5.1 and §5.4) account for the entire movement of this fork away from
pristine's frozen 16-seed table: `c6dcf24` leaves seeds 1 and 10 at −11.0274 and
−11.9092; `7ea27e1` moves them to −9.21059 and −11.968; the three subsequent
commits move nothing. The redocking result moved with them, from 0.66 Å to
0.70 Å top-1 backbone RMSD — a change well inside the ensemble's own spread and
far inside the 2.5 Å acceptance threshold.

---

## 6. Class IV — what was measured and deliberately reverted

Two class-A findings were implemented, measured, and **not** committed. They are
the clearest illustration of the corollary in §1: a correct criticism of a term
does not license changing it when the surrounding parameters absorbed its error.

### 6.1 The hydrophobic ramp rises where it should decay (A1)

At [src/energy.cpp:891](../src/energy.cpp#L891):

```c
if (distance > contact_cutoff + range) return 0.0;
if (distance < contact_cutoff)         return 1.0;
return (distance - contact_cutoff) / range;   /* rises 0 → 1 */
```

The evidence that this is a sign error is strong and entirely internal to the
file: the sibling `linear_decay()`, with the same signature and doc-comment
style, returns `1.0 − (distance − cutoff)/width`; both commented-out
alternatives, `hydrophobic_low_recip` and `hydrophobic_low_spline`, decay; and
the header comment above the function describes f(d) = d⁻¹, a decaying form.
Three of four candidate forms in the file decay. As written, the potential is
non-monotonic with two discontinuities: a pair in contact weighs 1, just past
contact ≈ 0, then climbs back to 1 at the far edge of the 2.8 Å range before
dropping to 0 — i.e. the model rewards partners for separating.

**Measurement.** Restoring `1.0 −` and re-running the full validation reshuffles
the pose ranking completely (top-5 becomes 15, 14, 8, 3, 12 instead of 15, 12,
13, 4, 10), makes energies systematically more negative (−29.06 vs −28.55 RT),
and moves the top-1 backbone RMSD on the one crystallographically-answerable
target **from 0.66 Å to 0.84 Å**.

**Interpretation.** `kauzmann_param = 0.122` is a contrastive-divergence fit,
and the code does not record which functional form it was fitted against. Since
the parameter is a pure linear prefactor on the hydrophobic term, a fit
performed against the rising form will have absorbed part of that form's error
into its magnitude. Correcting the shape without refitting the prefactor
exchanges a *known, characterised* bias for an *uncharacterised* one — and on
the only system with an experimental answer, it measurably loses. **Not
patched**; reopened as a recalibration.

### 6.2 Logarithmic compression applied to the electrostatic map (A4)

`lower_gridenergy()` ([src/energy.cpp:352](../src/energy.cpp#L352)) maps
E ↦ ln E + 1.718 for E > 2.718, and is applied inside `gridmap_initialise()`,
hence to all nine maps — including `e.map` (electrostatic potential) and `d.map`
(desolvation).

**The physical objection.** Because the transform fires only above +2.718, it
compresses the **positive** lobe of the Coulomb field while leaving the negative
lobe untouched. The resulting field is no longer antisymmetric under charge
inversion: a positive ligand atom near a positive receptor patch is penalised
substantially less than a negative atom near a negative patch — a systematic,
sign-dependent distortion of the electrostatic term. That reasoning is correct
and is not withdrawn.

**The measurement that stopped it.** A `soften` flag, false for the electrostatic
and desolvation slots, builds cleanly and passes all 22 tests. Four independent
blocks of 16 seeds each on the 3Q47 redock:

| seed block | 1–16 | 101–116 | 201–216 | 301–316 |
|---|---|---|---|---|
| compressed (this fork) | 0.70 | 0.84 | 0.78 | 0.74 |
| `e.map` exempt | 0.85 | **2.58** | **9.47** | 0.90 |

The failure is diagnostic: ranking every pose in each ensemble shows the
near-native pose is still **found** in every block — it stops **winning**. Under
the uncompressed map a 0.94 Å pose loses to a 9.47 Å pose. The search is
unaffected; the discrimination is destroyed.

**Why a 0.01% perturbation does this.** Across all 50 prepared targets, 1,467 of
15,192,578 `e.map` points exceed the threshold — 0.0097%, in 32 of 50 targets,
largest difference 1.043 kcal/mol·e⁻¹, about 0.36 kcal/mol once multiplied by
ADCP's backbone charges. That is the same order as the gaps separating adjacent
poses in a ranking formed as 0.25·E_total + 0.75·E_external. The compression was
in place when those weights and `kauzmann_param` were fitted, so it is not an
independent bias but part of the fitted model.

**And a genuine no-op.** For `d.map`, **0 of 15,192,578 points** across all 50
targets exceed the threshold. The compression provably never fires on the
desolvation map, so exempting it would change nothing at all — worth recording,
because it removes one third of the finding from the decision entirely.

**Status.** A1 and A4 are parked as a single joint recalibration: two continuous
knobs (`kauzmann_param`, `hydrophobic_cutoff_range`) and two discrete ones
(hydrophobic form, `e.map` compression). The fitting machinery is already in the
tree — `energy_probe_1()` is CRANKITE's contrastive-divergence rig, with slot 9
already bound to `kauzmann_param`. What is missing is not code but a **training
set disjoint from the 49-target benchmark**, and a *discriminative* objective:
contrastive divergence maximises the likelihood of native conformations, whereas
what failed here is ranking. [audit-fixes.md](audit-fixes.md) records the full
design, including the trap — the published ADCP numbers were produced *with*
A1's rising form, so "recovers the paper's numbers" and "correct" are not the
same test.

---

## 7. Class V — what is not adopted from upstream `cyclic`

Since comparison 5 established that the shipped ADFRsuite binary descends from
`cyclic`, the question "should this fork simply adopt `cyclic`?" is a real one.
The answer is no, and the reason is measurable.

**Upstream's `cyclic_energy()` collapses the macrocycle.** The intended restraint
is harmonic about the crystallographic Cα1–Cαn distance of 3.819 Å. The code on
that branch reads, in essence:

```c
double NCDistance = 0.0;
//NCDistance = distance(a->n, b->c);        // assignment commented out
if (CaDistance < 5) {
        ans += (sqrt(CaDistance) - 3.819)*(sqrt(CaDistance) - 3.819);
        ans += (sqrt(NCDistance) - 1.345)*(sqrt(NCDistance) - 1.345);
} else ans += CaDistance;
```

Two defects compound. First, `NCDistance` is never assigned, so the second term
is the constant (0 − 1.345)² = 1.809, independent of geometry. Second — and
decisively — `distance()` returns the **square** of the distance
(`vector.cpp`), while `CaDistance < 5` compares it against a *linear* threshold;
the test therefore means d < 2.24 Å, so the harmonic branch centred on 3.819 Å
almost never executes. What runs is `else ans += CaDistance`, i.e. d² — a
harmonic centred on **zero**, with no equilibrium distance at all. It pulls the
termini together until van der Waals repulsion stops them.

**The prediction this makes, and its test.** Such a term should compress every
backbone-cyclic macrocycle to the contact limit. Measured across all 18
backbone-cyclised targets in the benchmark, 432 poses each:

| source | median Cα1–Cαn | n |
|---|---|---|
| crystallographic ligand | **3.83 Å** | 18 |
| shipped ADFRsuite 1.0 binary | **2.52 Å** | 432 |
| this fork at HEAD | 3.99 Å | 432 |

The shipped binary closes the ring ~1.3 Å too tight on **every one of the 18
targets**, with no overlap with the crystallographic range. Two α-carbons at
2.5 Å is not a conformation: the C–C van der Waals contact distance is ≈3.4 Å.
This fork, which retains `master`'s harmonic centred on 3.819 Å, lands within
0.16 Å of the experimental median.

**Consequences for the port.** The remaining large hunks — `cyclic_energy()` and
`energy2()`'s new Cα–Cα harmonic on every sequence-adjacent pair — are the
likely source of the residual 4.5× speed advantage on cyclic targets. Buying
that speed by collapsing the macrocycle is not a trade worth making, so the port
stops there. Also unported: the single-shortest-pair disulfide search together
with the `0.25 · SSloss` reweighting (one design decision, undecided), and the
deletion of `transmutate()`'s early return on FIXED residues (upstream removed
it as a macOS compile fix; the removal is a behaviour change). Upstream's
rotamer cap *is* available here, as the opt-in `MaxRotamers` flag, defaulting to
the exhaustive scan this code has always performed — measured at 2.25× wall
clock (4825 s → 2140 s over the benchmark) at the cost of halving the
correlation with published per-target values (0.27 → 0.17), which is a trade the
user should make explicitly rather than inherit silently.

**A caveat this raises about the reference values themselves.** The cyclic
entries in the published tables — 38 per-target references transcribed into
`tests/validation/manifest.tsv` — were produced by a binary that compresses the
macrocycle by 1.3 Å. That is a concrete, testable explanation for the weak
agreement measured here (r = 0.19 between this fork and the shipped binary;
r = −0.07 between the shipped binary and the authors' own published values), and
it is more specific than "insufficient sampling". *It remains a hypothesis*:
confirming it requires checking whether the deposited poses carry the same
compression.

---

## 8. Discussion

Three patterns run through the changes above and are worth stating explicitly.

**Small numerical errors are not small in a Monte Carlo code.** §3.1 is the
sharpest case — a 10⁻⁷ relative error in one intermediate produces a completely
different trajectory after ~10⁵ steps, because the acceptance test is a
threshold and the chain is deterministic given its state. The practical
consequence is that bit-level trajectory comparison, not tolerance-based
comparison, is the only useful regression test for this class of code; the
16-seed ATOM-md5 sweep exists for that reason.

**Defects that only appear at production scale are systematically
under-detected.** The array underflow (§4.1) manifests only past 200,000 steps;
the precision loss (§3.1) only past ~200,000; the non-terminating draws (§4.3)
only when a long run stagnates. All three are invisible to short tests, and two
of the three produce *plausible output* rather than an obvious failure. Any
validation protocol for this code has to include at least one full-length,
production-settings run — the crash-rate map and the 2.5M-step redock exist to
cover that regime.

**"The term is wrong" and "the model is better without it" are different
claims.** §6 documents two cases where the first is well supported and the second
is false as measured. In a force field whose parameters were fitted against the
existing behaviour, an isolated correction moves the model *off* its calibration
point; the correct unit of work is the correction plus a refit, and the correct
thing to do in the absence of a training set disjoint from the benchmark is to
record the finding with its measurements and stop. Recording it — rather than
either patching it or quietly dropping it — is what keeps the finding actionable
when a fitting set does become available.

Finally, the four-way comparison changed the framing of this whole effort. Before
it, defects found in pristine `master` were naturally read as defects in "ADCP".
After it, two of them (§4.1's crash and §5.1's out-of-box penalty) are known to
be absent from the lineage that was actually shipped; one (§4.3) is present in
both upstream branches and remains unfixed upstream today; and one (§4.4's
nested-sampling draw) is present in *every* other binary measured — pristine,
`cyclic` and the shipped ADFRsuite build all still emit the uncorrected
sequence. Knowing which of those a given
finding is changes what should be done with it — ported, reported, or simply
noted — and that is why the comparison protocol is maintained as a standing,
re-runnable artefact rather than a one-off exercise.

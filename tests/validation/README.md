# ADCP validation dataset

A benchmark with **published reference values per target**, so a change to
`energy.cpp` / `metropolis.cpp` / `main.cpp` can be answered with a number
instead of an opinion.

Before this existed the repo validated on exactly one target (3Q47, redocking,
16 runs). That is a regression test, not a benchmark: it measures no success
rate, covers no long peptides, and never once exercised the cyclic-peptide path.

## The sets

| set | what | n | metric | reference |
|---|---|---|---|---|
| A | linear peptides, 16–20 aa | 11 | fnc | Bioinformatics 2019, Table 4 |
| B | backbone-cyclised peptides | 18 | fnc | JCTC 2019, Table 1 |
| C | disulfide-cyclised peptides | 20 | fnc | JCTC 2019, Table 2 |
| D | LEADS-PEP, ≥5 aa subset | 42 | bbRMSD ≤ 2.5 Å | JCIM 2016 SI *(not yet acquired)* |
| E | Glide SP-PEP, holo + apo | 19+10 | iRMSD | JCIM 2013 SI *(not yet acquired)* |
| F | peptiDB | 102 | — | no published ADCP number; coverage only |
| G | apo cross-docking of B and C | 28 | fnc | JCTC 2019, Table S3 *(not yet acquired)* |

Sets **A, B and C are transcribed verbatim** from the published tables, with a
value for every target at every rank the paper reports. `check_manifest.py`
recomputes every aggregate those tables print and asserts it — 40 checks. Run it
after any edit to `manifest.tsv`; it is the only thing that can catch a
transcription typo, because the per-target values have no second source.

It documents one genuine inconsistency in the source: Bioinformatics 2019
Table 4's footer reports 90.9% for top-3, while its own body gives 9/11 = 81.8%
(5N4B is 0.24 and 4RS9 is 0.49). The body is kept; the footer is not.

## Metric definitions

Read out of the authors' `clusterADCP.py`, not inferred:

- **ranking energy** `0.25·totalE + 0.75·extE`, × 0.59219 → kcal/mol
- **contact clustering** peptide CB (CA for Gly) vs receptor CB within 8 Å,
  Jaccard ≥ 0.8, leader clustering seeded by the best energy
- **fnc** non-hydrogen atom pairs within 5 Å of the reference,
  `|intersection| / |reference pairs|`
- **side chains are rebuilt** from the `Rotamers:` indices before contacts are
  counted; the docked output is backbone + pseudo-gamma only, so skipping the
  rebuild undercounts every contact
- **top N** means the best value among the first N *clusters* — not the Nth,
  and not the top N runs

## Protocol tiers

| tier | what | cost, 16 threads |
|---|---|---|
| `smoke` | 8 replicas × 200k steps | minutes |
| `reduced` | half the replicas, a quarter of the steps (⅛ CPU) | ~2.7 days for all sets |
| `full` | the published protocol, verbatim | ~22 days for all sets |

Published protocol: linear 80 × 3M steps/aa; long 80 × 7M steps/aa;
cyclic 300 × 1M steps/aa; box = crystallographic peptide + 4 Å padding.

**Only `full` is comparable with a published number.** `report.py` prints that
warning on every report that is not full-tier.

## Running it

```sh
./prepare_all.sh  targets ABC              # fetch + build .trg for sets A,B,C
./run_set.sh      ../../build/src/adcp targets results smoke ABC
for d in results/*/; do ./score_ref.sh "$d" nc 0.8; done
python3 report.py results --tsv report.tsv
```

`prepare_all.sh` and `run_set.sh` are both idempotent and resumable, so a
multi-day run survives an interruption.

## Target preparation, and why it does not reproduce the authors' targets

`prepare_target.sh` fetches the RCSB entry, picks the peptide chain, and drives
`prepare_receptor` + `agfr` from ADFRsuite. Three things it handles that are
not obvious, each found by a target that failed without it:

- **altLoc must be blanked, not just filtered.** Keeping the `A` in column 17
  makes `prepare_receptor` build 5-character hydrogen names (`HE2A` → `HE2A1`)
  that overflow the atom-name field and shift the residue name one column left.
  agfr then dies with "invalid or missing coordinate(s)" hundreds of lines from
  the cause.
- **Heteroatoms AutoDockTools cannot type** (3WNF's cadmium) do not make
  `prepare_receptor` fail — it writes a chargeless PDBQT line that agfr rejects
  later. The script asks its log which residues it could not type, drops those,
  and retries. They have no AutoGrid map anyway, so they could not have
  contributed.
- **`translationPoints.npy` dtype is not fixed.** The published 3Q47 target is
  float64, which is all `tests/fetch_target.sh` knows how to read; the agfr on
  this machine emits float32.
- **Every staged map is checked against its own `NELEMENTS` header.** agfr can
  die inside `addGradientToMaps` *after* writing a `.trg`, leaving maps that are
  truncated or padded. Accepting "the .trg exists and is non-empty" passed 6 of
  49 targets whose maps disagreed with their own header by thousands of values.
  A corrupt map does not crash the docking, it silently changes the potential —
  the worst failure mode a benchmark can have. On failure the stage and the
  `.trg` are deleted so a partial target cannot be mistaken for a good one.
  (All 6 succeeded on retry, so those failures were transient contention, not
  bad receptors. The check stays regardless.)

### Known limitation: metals are dropped

The untypeable-heteroatom retry removes what `prepare_receptor` cannot charge,
and that currently includes **Mg (5UWI, 6CIT), Zn (4OU3), Ca, Li and Cd**.
`provenance.txt` records the dropped residues per target.

This is a real modelling decision, not just plumbing. ADCP only ever reads nine
map types (C, A, SA, N, NA, OA, HD, d, e) — there is no Mg or Zn map — so a
metal could never be felt directly. But `e.map` and `d.map` are computed from
the whole receptor, and a Zn²⁺ contributes strongly to the electrostatic map.
Dropping it removes that contribution. Giving these ions a formal charge instead
of deleting them is the obvious improvement; it is not done yet.

Calibration against the one target we can compare: our 3Q47 build gives box
`54 × 54 × 52` with 193 AutoSite points against the published `54 × 56 × 56`
with 106, and redocks to **0.78 Å** where the published target gives **0.66 Å**.
Our `-b user` box is byte-identical to what `agfr -b ligand -P 4` produces, so
the difference is in the authors' preparation, not ours. `provenance.txt` in
every target records the box, the point count, the dropped heteroatoms and the
chain that was chosen.

The published length is not always the chain's residue count — JCTC 2019
Table 1 calls the group-I HIV integrase binders 6-mers while 3AV9 chain X is the
8-mer `SAKIDNLD`, apparently counting the macrocycle. A mismatch is recorded in
`peptide.chain` and reported, never fatal.

## Two things to fix before spending the full-tier budget

- **`scripts/runADCP.py` does no clustering.** It ranks the N runs by energy and
  prints the top 5; the upstream wrapper calls `clusterADCP` with `-nc`/`-rmsd`/
  `-ref`. Until that is closed, "top 10" here is top 10 *runs* and is not the
  paper's "top 10". `score_ref.sh` works around it for the benchmark, but the
  production wrapper still ships without it.
- **The Met→C.map bug** (audit finding A2): `hasCYS` in `main.cpp` only tests
  `'C'`, so a peptide with methionine and no cysteine scores its sulfur against
  the carbon map. A large fraction of any real benchmark contains Met. Running
  the full tier with that live measures the bug, not the engine.

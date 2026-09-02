# Porting upstream's `cyclic` branch into this codebase

## Context

This repository is a C++17 port of `ccsb-scripps/ADCP` **`master`** (`1c1a330`).
Upstream has a second public branch, **`cyclic`** (97 commits, tip `1d34035`,
2019-06-11), which `master` never merged. Evidence that `cyclic` — not `master` —
is the lineage of the binary the authors actually ship in ADFRsuite 1.0:

- the binary's `DW_AT_comp_dir` is `/data/mgltools/MakeInstallers2/src/Crankite_AD`;
  the `cyclic` branch's merge commits push to `orca:/mnt/raid/services/git/Crankite_AD`
- binary is dated July 2019, `cyclic` tip is 2019-06-11 (`master`'s last code
  change is 2019-02-28)
- measured behaviour matches (below); `master` does not

Measured on this machine, same targets, same seeds, 200k steps:

| target | cyclization | official binary | `cyclic` built here | our HEAD (from `master`) |
|---|---|---|---|---|
| 2ivz | none | 9.7 s | 10.8 s | 10.5 s |
| 3av9 | backbone | 6.1 s | 9.1 s | **25.5 s** |
| 4m1d | disulfide | 14.4 s | 11.8 s | **54.8 s** |
| 1sfi | backbone+ss | 15.5 s | 16.4 s | **64.5 s** |

Linear docking is unaffected. **Cyclic docking is 3–4× slower in our tree**, and
that gap is the `cyclic` branch. It is not a regression introduced by the C++
port.

DWARF `decl_line` drift between the shipped binary and each branch confirms the
same three files carry the difference — `energy.c`, `main.c`, `metropolis.c` —
and that `cyclic` is much closer to the binary than `master` (max drift +9 vs
+65 lines in `energy.c`). No public commit matches exactly, so the shipped
binary is a few internal commits past the public branch tip.

## What actually differs

`git diff master origin/cyclic` — 6 files, 184 insertions, 178 deletions.
Sorted by what the change *is*, not by file:

### Group N — noise, port verbatim, no validation needed

| file | change |
|---|---|
| `canonicalAA.c` / `.h` | the 19-line LGPL header is **absent** in `cyclic` (it branched before `master`'s "added license stuff"). **Keep our headers.** Do not port. |
| `main.c` @880 | whitespace only inside `AD_init` |
| `metropolis.c` @71 | comment only |
| `energy.c` @347 | comment only, inside `lower_gridenergy` |

### Group P — performance, no intended science change

| where | change |
|---|---|
| `energy.c` `scoreSideChain` / `scoreSideChainNoClash` (6 hunks) | rotamer scoring rework |
| `energy.c` `gridenergy` @1817, @1901 | index/bounds handling around the trilinear stencil |
| `energy.c` `ADenergyNoClash` @2153 | driver-level change |
| `energy.c` `indMoved` @1927 | −3 lines |

These are the likely source of the 3–4×. They must still be validated as
science-neutral, because "intended" is not "verified".

### Group S — real science changes, each needs its own decision

| where | change | note |
|---|---|---|
| `energy.c` `sbond_energy` @1217 (+31, the largest hunk) | the greedy nearest-neighbour disulfide search is replaced by picking **only the single shortest Cys–Cys pair**, with an early `return 0.0` when none is found | changes which disulfides are scored |
| `metropolis.c` @182 | `loss = loss + SSloss + externalloss` → `loss = (loss + 0.25 * SSloss) + externalloss` | **scales the disulfide term in the Metropolis test by 0.25**; commit is literally "0.25 for -S-S-" |
| `energy.c` `gridenergy` @1828 (+15) | out-of-box penalty rewritten: measured from the **nearest face** instead of the box centre, plus `if (outofBox) return erg > 10000 ? 10000 : erg/100;` | **this is audit finding A5, fixed** |
| `energy.c` `cyclic_energy` @2550 | cyclic bond term reworked | |
| `energy.c` `energy2` @2505, @2515 | pair-term dispatch | |
| `energy.c` `lowlevel_sbond` @1156, @1179 | disulfide geometry | |
| `energy.c` `hbond` @758 | one line | |
| `metropolis.c` `transmutate` @258 | `void` → `int`, and the early `return` on FIXED residues is **deleted** | MacOS compile fix; the deletion is a behaviour change |
| `metropolis.c` `transmove` @373, @542 | −4 and −5 lines | |
| `metropolis.c` `crankshaft` @895 | one line | |
| `metropolis.h` @8 | `transmutate` signature | follows the `void`→`int` change |
| `main.c` @32, @210, @378, @386 | `simulate()` internals | |

## What the port does NOT fix

Re-audited every class-A finding against `cyclic`. Only one is addressed:

| finding | `master` (ours) | `cyclic` |
|---|---|---|
| A1 inverted hydrophobic ramp | present | **identical** — no hunk in that region |
| A2 methionine scored on the carbon map | present | **identical** — the `AD_init` hunk is whitespace |
| A3 unbounded "strictly better" swap loops | present | **identical code**; the hang is latent, not fixed |
| A4 log compression applied to `e.map`/`d.map` | present | **identical** — the `lower_gridenergy` hunk is a comment |
| A5 out-of-box penalty from the box centre | present | **FIXED** |
| A6 `ramabias` index without clamp | present | **identical** |

A3 deserves emphasis: the five `while (swapEnergy[swapInd] ...)` loops are
byte-identical in both branches. `cyclic` did not hang in 48 replicas
(16 seeds × 3 targets) only because its different energy code takes different
trajectories. **The fix in this repo is not redundant with upstream.**

## Approach

A patch cannot be applied: our tree is `.cpp`, uses `std::vector`, `View3<>`,
`Chain::aa` as a container, and has ~15 memory/UB fixes upstream lacks. Every
hunk is ported **semantically**, by hand, against our current code.

Order, chosen so each step is separately validatable:

1. **Group N** — skip entirely except confirming our license headers stay.
2. **Group S / A5 only** — the out-of-box rewrite. Small, self-contained, and
   independently justified by the audit. Land it first as its own commit.
3. **Group P** — the performance hunks, one function at a time. Gate: the
   fold-determinism test and a bit-level docking probe must stay green, and the
   cyclic timing must actually drop. If a hunk changes results, it moves to
   Group S.
4. **Group S / disulfide** — `sbond_energy`, `lowlevel_sbond`, and the
   `0.25 * SSloss` weight, together, since they are one design decision.
5. **Group S / remainder** — `cyclic_energy`, `energy2`, `transmutate`,
   `transmove`, `crankshaft`, `simulate`.

Steps 4 and 5 change docking results by construction. They need a **decision**,
not just a port: upstream's shipped binary behaves this way, so matching it is
defensible, but it invalidates every number this repo has produced so far.

## Validation

The benchmark built in `tests/validation/` is the instrument. Per step:

- `ctest` (smoke + functional), and `func_adcp_fold_determinism` 12×, per
  MIGRATION.md's standing rule for `energy`/`metropolis`/`main` changes
- `val_3q47_redock` must still reach 0.66 Å
- **Group P gate**: run `tests/validation/run_set.sh` at `smoke` over sets A/B/C
  before and after, and require `report.py` to show *zero* changed fnc values.
  Anything that moves is not a performance change.
- **Group S gate**: the same run, but the expectation is inverted — results
  *will* move, and the diff must be recorded per target in a
  `docs/compares/`-style document
- cyclic timing measured on 3av9 / 4m1d / 1sfi, targeting the official binary's
  6.1 / 14.4 / 15.5 s
- differential oracle: `~/ADFRsuite-1.0/bin/adcp_Linux-x86_64` on the same
  targets and seeds; after Group P the numbers should converge toward it

## Risks

| risk | mitigation |
|---|---|
| `cyclic` is 4 commits behind the shipped binary, so parity is not reachable from public source | accept; target behavioural parity (timing, no hang), not bit parity |
| Group S changes make every result in `report_fix.tsv` obsolete | land Groups N/P first, freeze a report, then treat Group S as a deliberate re-baseline |
| The disulfide rework may interact with our A3 fix in `simulate()` | port A3's bound first (already done), re-verify inertness after each disulfide hunk |
| `transmutate` losing its FIXED-residue guard | our tree has fixed-residue support paths upstream may not exercise; keep the guard and port only the `void`→`int` change unless a test demands otherwise |
| No upstream test coverage for any of this | the 49-target benchmark is the coverage; it did not exist before this week |

## Out of scope

- Merging `cyclic`'s missing LGPL headers (we keep ours)
- Chasing the 4 internal commits between the public tip and the shipped binary
- Fixing A1/A2/A4/A6, which this port does not touch and which remain open

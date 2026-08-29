# AutoDock CrankPep

AutoDock CrankPep or ADCP is an AutoDock docking engine specialized for docking peptides.
It combines technology from the protein folding field with an efficient representation of
a rigid receptor as affinity grids to fold the peptide in the context of the energy landscape
created by the receptor. A Monte-Carlo search is used to fold the peptide while concurrently
optimizing the interaction between the peptide and the receptor molecule, yielding docked
peptides. The program can dock peptides provided as 3D structures in PDB files or from a
sequence string. It has been shown to successfully re-dock peptides with up to 20 amino acids in length.

ADCP is developed based on CRANKITE.
Podtelezhnikov, A.A. and Wild, D.L. (2008) *Source Code Biol. Med.*, 3, 12.

ADCP is available under the GNU LGPL v2.0 OpenSource license.
Please visit [adcp.scripps.edu](https://adcp.scripps.edu) for more details.

## Usage

There are two ways of using ADCP.

The first and recommended way is to use the prepared Python wrapper `runADCP.py`.
`runADCP` works directly with a zipped target file (`.trg`) prepared by AGFR.
Please see [adcp.scripps.edu/tutorial](https://adcp.scripps.edu/tutorial) for more details.

The other way is to use the compiled binary directly.
Note this approach requires unzipped map files within the execution folder.

Example usage:

```
adcp -r 10000x3000 -t 2 apgvgvapgvgv -p Bias=NULL,external=5,con8,2,1.0,external2=4,con8,2,1.0,Opt=1,0.5,0.5,-0.5 -o output.pdb
```

(The binary is built as `adcp`; older docs call it `adcp_Linux-x86_64`.)
The parameter set above documents the syntax but is not what the Python wrapper
runs in practice -- see `README_AD` for the production command line.

### `external=5,con8,2,1.0`

This calls the autodock grid maps. It will look for `rigidReceptor.*.map`.
- `con8` is the file indicating the residues to score, usually just all residue numbers.
- `2` is a thermo factor that is an addition to the original temperature.
- `1.0` calls the reconstruction of side chains during the scoring.

### `external2=4,con8,2,1.0`

This calls the cyclic procedure to create an artificial peptide bond between the first and last residue.
- `con8` is meaningless right now, but can be used to indicate the location of the cyclic bond in the future.
- `2` is the thermo factor, and it will overwrite the previous one in `external`.
- `1.0` is meaningless right now, but can be used to flag the cyclic bond type, etc.

### `Opt=1,0.5,0.5,-0.5`

This calls the optimizing/docking procedure. `Opt=0` is to use regular MC.
The following 3 numbers are the weight of the target energy in the optimizing procedure.
In the example above:

```
targetE = 0.5 * totalE + 0.5 * externalE + (-0.5) * firstlastE
```

Note that `totalE` includes `internalE`, `externalE` and `firstlastE`. So the above weights
scale down the internal energy by a factor of 2 and remove the energy between the first and
last residues. The default weights are set to be `1, 0, 0`.

## Building

Requires CMake >= 3.16 and a C99 compiler.

```
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

Binaries land in `build/src/adcp` and `build/tools/`.

Build options:

| Option | Default | Description |
|---|---|---|
| `ADCP_MPI` | OFF | MPI parallel tempering (builds `adcp_mpi`) |
| `ADCP_OPENMP` | OFF | OpenMP support (used by cdlearn) |
| `ADCP_TOOLS` | ON | build the auxiliary tools |
| `ADCP_DOCKING_TESTS` | OFF | enable the docking tests (see below) |
| `ADCP_LEGACY_COMMON` | OFF | restore `-fcommon`; only for legacy compilers |

```
cmake -B build -DCMAKE_BUILD_TYPE=Release -DADCP_MPI=ON -DADCP_OPENMP=ON
```

## Testing

```
ctest --test-dir build --output-on-failure
```

That runs the 16 tests that need nothing but the source tree. Add
`-DADCP_DOCKING_TESTS=ON` at configure time for the 3 that also exercise
docking, for 19 total.

Tests are grouped by label:

| Label | Tests | Time | Network |
|---|---|---|---|
| smoke | 11 | <1 s | no — tool binaries vs a committed PDB |
| functional | 5 | ~2 min | no — folding, nested sampling, checkpointing, and the analytic-evidence check |
| docking | 2 | ~20 s | first run — 3Q47 grid path |
| validation | 2 | ~2 min | first run — full redocking vs the crystal pose |

```
ctest --test-dir build -L smoke
ctest --test-dir build -L "smoke|functional"
```

The docking and validation rows both count `dock_fetch_target`: CTest pulls
that setup test in automatically as a fixture. There are 3 distinct docking
tests, not 4.

Useful flags:

```
ctest --test-dir build -N                    # list tests without running
ctest --test-dir build -R dock               # run tests matching a regex
ctest --test-dir build -V -R val_3q47_redock # verbose: show the RMSD table
ctest --test-dir build -j4                   # run in parallel

cmake --build build --target check           # rebuild, then run everything
```

### Docking tests

ADCP's docking path needs AutoGrid maps, which normally come from AGFR in
the ADFRsuite. The ADCP documentation publishes an already-prepared target
for PDB 3Q47, so the suite is not required:

```
cmake -B build -DCMAKE_BUILD_TYPE=Release -DADCP_DOCKING_TESTS=ON
cmake --build build -j$(nproc)
ctest --test-dir build --output-on-failure -L "docking|validation"
```

On first run this downloads `3Q47.trg` (~10 MB) into `build/tests/target/` and
caches it there. Nothing large is added to the repository. If the download
is unavailable the docking tests report SKIP rather than failing, so an
offline or firewalled machine still gets a green run of the other 12.

`tests/fetch_target.sh` stages the maps exactly the way `scripts/runADCP.py`
does, so the tests exercise the real production layout.

What they check:

- **`dock_3q47_smoke`** — the grid path runs to completion, the maps and all
  106 translation points load, the external energy is not zero (the folding
  path reports exactly 0.000000), and a fixed seed reproduces byte for byte.

  It runs 250k steps on purpose. A stack overflow in the swap pool used to
  crash every production docking run, but only past 200k steps -- a shorter
  test passes against the broken binary.

- **`val_3q47_redock`** — the scientific check. Runs 16 independent searches
  at 2.5M steps, ranks them by target energy the way `runADCP.py` does, and
  measures the backbone RMSD of the top-ranked pose against the
  crystallographic ligand. No superposition is applied: ADCP docks into a
  fixed receptor frame, so both are already in the same frame.

  Currently redocks the native peptide (`npisdvd`) to 0.66 Å, against a
  2.5 Å threshold.

Reading the validation output:

```
ctest --test-dir build -V -R val_3q47_redock
```

It prints the top five runs by energy, the top-ranked pose's RMSD, and the
best RMSD found anywhere in the ensemble. Those last two are worth comparing:
if the ensemble contains a good pose but the top-ranked one is poor, the
search is working and the scoring function is what needs attention.

#!/usr/bin/env python3
"""CA1-CAn distance for backbone-cyclised targets: what the ring actually closes to.

cyclic_energy() restrains the first and last alpha carbons to 3.819 A, the
trans-peptide CA-CA distance. Whether the search honours that is not something
the energy function can tell you -- it has to be read off the docked poses.

Reported against the crystallographic ligand of the same target, so a
systematic compression shows up as a shift of the whole distribution rather
than as a per-target discrepancy.

Usage: ring_geometry.py <results_dir> [pdb ...]
       ring_geometry.py <results_dir_a> --vs <results_dir_b>
"""
import sys, os, glob, math, statistics


def last_model_atoms(path):
    """Atoms of the final MODEL in a multi-model pdb, or all ATOMs if unmodelled.

    The final MODEL is the pose ADCP settled on; earlier ones are the
    improvement trajectory. Matches how run_redock_validation.sh reads a run.
    """
    cur, keep, saw_model = [], [], False
    for l in open(path):
        if l.startswith('MODEL'):
            saw_model, cur = True, []
        elif l.startswith(('ATOM', 'HETA')):
            cur.append(l)
        elif l.startswith('ENDMDL'):
            keep = cur
    return (keep or cur) if saw_model else cur


def ca_span(path):
    """Distance between the CA of the lowest- and highest-numbered residue.

    Residues are ordered numerically, not by appearance: pdbqt writes atoms in
    torsion-tree order, so first-seen is not first-in-sequence.
    """
    ca = {}
    for l in last_model_atoms(path):
        if l[12:16].strip() != 'CA':
            continue
        ca[int(l[22:26])] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    if len(ca) < 2:
        return None
    lo, hi = ca[min(ca)], ca[max(ca)]
    return math.dist(lo, hi)


def target_spans(tdir):
    """(crystal, [pose, ...]) for one target directory."""
    native = os.path.join(tdir, 'native_pep.pdb')
    crystal = ca_span(native) if os.path.exists(native) else None
    poses = []
    for run in sorted(glob.glob(os.path.join(tdir, 'run_*.pdb'))):
        rc = os.path.join(tdir, 'run_%s.rc' % os.path.basename(run)[4:-4])
        if os.path.exists(rc) and open(rc).read().strip() != '0':
            continue          # crashed or timed out: not a pose
        d = ca_span(run)
        if d is not None:
            poses.append(d)
    return crystal, poses


def collect(results, only=None):
    rows = []
    for tdir in sorted(glob.glob(os.path.join(results, '*/'))):
        pdb = os.path.basename(tdir.rstrip('/'))
        if only and pdb not in only:
            continue
        if not os.path.exists(os.path.join(tdir, 'done')):
            continue
        crystal, poses = target_spans(tdir)
        if poses:
            rows.append((pdb, crystal, poses))
    return rows


def med(xs):
    return statistics.median(xs) if xs else float('nan')


def main():
    args = [a for a in sys.argv[1:]]
    other = None
    if '--vs' in args:
        i = args.index('--vs')
        other = args[i + 1]
        args = args[:i] + args[i + 2:]
    results, only = args[0], set(args[1:]) or None

    rows = collect(results, only)
    if not rows:
        sys.exit('no finished targets with poses under %s' % results)
    brows = {p: (c, q) for p, c, q in collect(other, only)} if other else {}

    hdr = '%-8s %8s %8s %5s' % ('target', 'crystal', 'poses', 'n')
    if other:
        hdr += ' %8s %8s' % ('other', 'delta')
    print(hdr)
    print('-' * len(hdr))

    all_crystal, all_poses, all_other = [], [], []
    for pdb, crystal, poses in rows:
        line = '%-8s %8s %8.2f %5d' % (
            pdb,
            '%.2f' % crystal if crystal else '   n/a',
            med(poses), len(poses))
        if other:
            oc, op = brows.get(pdb, (None, []))
            line += ' %8s %8s' % (
                '%.2f' % med(op) if op else '   n/a',
                '%+.2f' % (med(poses) - med(op)) if op else '     n/a')
            all_other += op
        print(line)
        if crystal:
            all_crystal.append(crystal)
        all_poses += poses

    print('-' * len(hdr))
    print('median over all poses      %8.2f A  (n=%d)' % (med(all_poses), len(all_poses)))
    if all_other:
        print('median, other arm          %8.2f A  (n=%d)' % (med(all_other), len(all_other)))
    print('median, crystal ligands    %8.2f A  (n=%d)' % (med(all_crystal), len(all_crystal)))
    print('restraint equilibrium         3.82 A')


if __name__ == '__main__':
    main()

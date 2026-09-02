#!/usr/bin/env python3
"""Turn scored results into a per-target comparison against the published numbers.

Reads manifest.tsv for the reference values and each <results>/<pdb>/score_ref_nc.txt
for ours, then prints one row per target plus the aggregate the papers report.

"Top N" means what the papers mean by it: the best value among the first N
clusters, not the Nth cluster. Verified against JCTC 2019 Table 2 semantics.

Usage: report.py <results_dir> [--tsv out.tsv]
"""
import sys, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
# Ranks the papers report. Bioinformatics 2019 Table 4 uses 1/3/5/20/all;
# JCTC 2019 Tables 1-2 use 1/5/10/100. Compute all of them, compare only where
# the manifest actually carries a published value.
RANKS = (1, 3, 5, 10, 20, 100)


def load_manifest():
    rows = {}
    with open(os.path.join(HERE, 'manifest.tsv')) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if f[0] == 'set':
                continue
            ref = {}
            if f[6] != '-':
                ref = {k: float(v) for k, v in (p.split('=') for p in f[6].split(';'))}
            rows[f[1]] = dict(set=f[0], pdb=f[1], length=int(f[3]), cyc=f[4],
                              metric=f[5], ref=ref, source=f[7])
    return rows


def parse_clusters(path):
    """fnc of each cluster, in rank order, from clusterADCP's -nc table."""
    out = []
    for line in open(path):
        m = re.match(r'\s*(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+)\s', line)
        if m:
            out.append((int(m.group(1)), float(m.group(2)), float(m.group(3)),
                        int(m.group(4))))
    out.sort(key=lambda r: r[0])
    return out


def best_by_rank(clusters):
    """Best fnc within the top N clusters, for each N the papers report."""
    got = {}
    for n in RANKS:
        head = [c[2] for c in clusters[:n]]
        got['top%d' % n] = max(head) if head else None
    got['all'] = max((c[2] for c in clusters), default=None)
    return got


def main():
    results = sys.argv[1] if len(sys.argv) > 1 else 'results'
    tsv_out = None
    if '--tsv' in sys.argv:
        tsv_out = sys.argv[sys.argv.index('--tsv') + 1]

    man = load_manifest()
    rows, missing = [], []
    for pdb, meta in man.items():
        score = os.path.join(results, pdb, 'score_ref_nc.txt')
        if not os.path.exists(score):
            missing.append(pdb)
            continue
        clusters = parse_clusters(score)
        if not clusters:
            missing.append(pdb)
            continue
        info = {}
        ipath = os.path.join(results, pdb, 'run_info.txt')
        if os.path.exists(ipath):
            for l in open(ipath):
                k, _, v = l.partition(' ')
                info[k.strip()] = v.strip()
        rows.append(dict(meta, got=best_by_rank(clusters), nclust=len(clusters),
                         tier=info.get('tier', '?'),
                         sequence=info.get('sequence', ''),
                         replicas=info.get('replicas', '?'),
                         steps=info.get('steps', '?')))

    if not rows:
        print('no scored targets found under %s' % results)
        print('run run_set.sh then score_ref.sh first')
        return 1

    tiers = sorted({r['tier'] for r in rows})
    print('ADCP validation report — %d targets scored, tier(s): %s'
          % (len(rows), ', '.join(tiers)))
    if tiers != ['full']:
        print('NOTE: only the "full" tier runs the published protocol. Numbers from')
        print('      any other tier are not comparable with the published columns.')
    print()

    keys = [k for k in ['top%d' % n for n in RANKS] + ['all']
            if any(k in r['ref'] for r in rows)]
    hdr = '%-5s %-6s %-20s %3s ' % ('set', 'pdb', 'sequence', 'aa')
    for k in keys:
        hdr += ' %13s' % ('%s o/pub' % k)
    print(hdr)
    print('-' * len(hdr))

    agg = {k: {'ours': [], 'pub': []} for k in keys}
    for r in sorted(rows, key=lambda x: (x['set'], x['pdb'])):
        line = '%-5s %-6s %-20s %3d ' % (r['set'], r['pdb'],
                                         r.get('sequence', '')[:20], r['length'])
        for k in keys:
            ours, pub = r['got'].get(k), r['ref'].get(k)
            if ours is not None and pub is not None:
                agg[k]['ours'].append(ours)
                agg[k]['pub'].append(pub)
            line += ' %6s/%-6s' % ('%.3f' % ours if ours is not None else '-',
                                   '%.2f' % pub if pub is not None else 'NA')
        print(line)

    print()
    print('%-34s' % 'aggregate (matched targets only)', end='')
    for k in keys:
        print(' %13s' % k, end='')
    print()
    for label, key in (('avg fnc, ours', 'ours'), ('avg fnc, published', 'pub')):
        print('%-34s' % label, end='')
        for k in keys:
            v = agg[k][key]
            print(' %13s' % ('%.3f' % (sum(v) / len(v)) if v else '-'), end='')
        print()
    for label, key in (('fnc >= 0.5, ours', 'ours'), ('fnc >= 0.5, published', 'pub')):
        print('%-34s' % label, end='')
        for k in keys:
            v = agg[k][key]
            print(' %13s' % ('%d/%d' % (sum(1 for x in v if x >= 0.5), len(v))
                             if v else '-'), end='')
        print()

    if missing:
        print()
        print('not scored yet (%d): %s' % (len(missing), ' '.join(sorted(missing))))

    if tsv_out:
        with open(tsv_out, 'w') as fh:
            fh.write('set\tpdb\tlength\ttier\treplicas\tsteps\tnclusters\t'
                     + '\t'.join('ours_%s\tpub_%s' % (k, k) for k in keys) + '\n')
            for r in sorted(rows, key=lambda x: (x['set'], x['pdb'])):
                cells = [r['set'], r['pdb'], str(r['length']), r['tier'],
                         r['replicas'], r['steps'], str(r['nclust'])]
                for k in keys:
                    cells.append('%.4f' % r['got'][k] if r['got'].get(k) is not None else '')
                    cells.append('%.2f' % r['ref'][k] if r['ref'].get(k) is not None else '')
                fh.write('\t'.join(cells) + '\n')
        print('\nwrote %s' % tsv_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())

"""Does the 50-site chrX floor matter? Subsampling test on CEPH 1463.

The sex check calls a sample only when the whole Wilson 95% interval of its
chrX heterozygous fraction lies on one side of the cohort boundary; samples
with fewer than ``min_x_sites`` usable chrX variant sites are skipped outright.
The floor (50) is a convention, not a derived value.  This script asks what it
changes: each sample's non-reference chrX sites are subsampled to n sites, and
sex is inferred with the tool's own ``sex_evidence`` under the default floor
(50) and with no floor (1).  Every call is compared with the published pedigree.

Data (not vendored):  git clone --depth 1 https://github.com/brentp/peddy /tmp/peddy
Run:   python validation/sex_floor_subsampling.py [--data /tmp/peddy/data] [--reps 200]
Output: validation/results/sex_floor_subsampling.csv and a printed summary.
"""
from __future__ import annotations

import argparse, csv, random, sys
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from admissible.checks.identity import IdentityConfig, SampleSets, build_sets, sex_evidence  # noqa: E402
from admissible.contigs import normalize_contig  # noqa: E402
from admissible.ped import Sex, read_ped  # noqa: E402
from admissible.vcfio import load_cohort  # noqa: E402

SIZES = [10, 20, 30, 50, 100, 200, 400]
FLOORS = {"floor 50 (default)": 50, "no floor": 1}
BUILDS = {"build undeclared (as in the paper)": None, "GRCh37, PAR and XTR excluded": "GRCh37"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="/tmp/peddy/data")
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20261111)
    args = ap.parse_args()
    data = Path(args.data)
    matrix = load_cohort([data / "ceph1463.peddy.vcf.gz"])
    ped = read_ped(data / "ceph1463.ped")
    truth = {s: ped.individuals[s].sex for s in matrix.samples if s in ped.individuals}
    sets = build_sets(matrix)
    x_ids = {i for i, (chrom, *_r) in enumerate(matrix.index.keys) if normalize_contig(chrom) == "X"}

    rows, rng = [], random.Random(args.seed)
    for bname, build in BUILDS.items():
        # sites the tool would use under this build: run once at full depth to learn them
        full, _, _ = sex_evidence(matrix, sets, IdentityConfig(), build)
        for n in SIZES + ["all"]:
            for fname, floor in FLOORS.items():
                cfg = replace(IdentityConfig(), min_x_sites=floor)
                tally = {"correct": 0, "wrong": 0, "undetermined": 0}
                for _ in range(1 if n == "all" else args.reps):
                    sub = {}
                    for s, ss in sets.items():
                        xs = sorted(ss.nonref & x_ids)
                        keep = set(xs) if n == "all" or len(xs) <= n else set(rng.sample(xs, n))
                        drop = (ss.nonref & x_ids) - keep
                        sub[s] = SampleSets(ss.assessed, ss.het - drop, ss.homalt - drop, ss.homref)
                    ev, _, _ = sex_evidence(matrix, sub, cfg, build)
                    for s, e in ev.items():
                        if s not in truth or truth[s] is Sex.UNKNOWN:
                            continue
                        k = "undetermined" if e.inferred is Sex.UNKNOWN else ("correct" if e.inferred is truth[s] else "wrong")
                        tally[k] += 1
                total = sum(tally.values())
                rows.append({"build": bname, "sites_per_sample": n, "floor": fname, "calls": total, **tally,
                             "wrong_pct": round(100 * tally["wrong"] / total, 2),
                             "undetermined_pct": round(100 * tally["undetermined"] / total, 2)})
        fr = sorted(round(e.x_het_frac, 3) for s, e in full.items() if truth.get(s) is Sex.FEMALE and e.n_x_sites)
        print(f"{bname}: usable chrX sites per sample {min(e.n_x_sites for e in full.values())}-"
              f"{max(e.n_x_sites for e in full.values())}; female het fractions {fr}")

    out = ROOT / "validation" / "results" / "sex_floor_subsampling.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(f"\n{'build':38s} {'n':>5s} {'floor':20s} {'wrong%':>7s} {'undet%':>7s}")
    for r in rows:
        print(f"{r['build']:38s} {str(r['sites_per_sample']):>5s} {r['floor']:20s} {r['wrong_pct']:7.2f} {r['undetermined_pct']:7.2f}")
    print(f"\nwritten: {out}")


if __name__ == "__main__":
    main()

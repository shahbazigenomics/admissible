"""Calibrate the identity thresholds of admissible on many families (issues #5 and #6).

Two questions, one script:

* #6  Where should the parent-offspring IBS0 limit sit?  For real parent-offspring
      pairs, real full sibs, simulated full sibs, second-degree and unrelated pairs,
      IBS0 per het call is measured exactly as admissible measures it, under the
      conditions a real delivery can have: whole chromosome or exome target only,
      all calls or PASS only, and 0 / 0.1 / 0.5 / 1% added genotype error.
* #5  What separates a consanguineous female from a male on chrX?  For every sample,
      chrX heterozygosity (PAR and XTR excluded), autosomal heterozygosity and the
      inbreeding coefficient F are measured, and candidate sex rules are scored.

The metrics are re-implemented with numpy so that thousands of samples fit in
memory; ``selfcheck`` proves the re-implementation gives the same numbers as
admissible's own code on CEPH 1463 before anything else is trusted.

Subcommands
-----------
selfcheck  --vcf ceph1463.peddy.vcf.gz --ped ceph1463.ped --build GRCh37
pairs      --vcf <genotypes> --ped <ped> [--phased-vcf <phased>] [--bed exome.bed] --out DIR
sex        --x-vcf <chrX> --auto-vcf <chr22> [--auto-vcf ...] --ped <ped> [--bed exome.bed] --out DIR

Requires: numpy, cyvcf2 (pip install cyvcf2), and admissible importable (run from the
repository root, or pip install -e .).
"""
from __future__ import annotations

import argparse, csv, gzip, math, random, sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

HOMREF, HET, MISS, HOMALT = 0, 1, 2, 3          # cyvcf2 gt_types
ERRORS = (0.0, 0.001, 0.005, 0.01)
CHUNK = 20000


# ----------------------------------------------------------------- inputs
def read_ped(path):
    """Accepts a 6-column PED or the 1000 Genomes table with a header
    (FamilyID SampleID FatherID MotherID Sex Population Superpopulation)."""
    rows = {}
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f if l.strip()]
    head = lines[0].split()
    if "SampleID" in head:
        idx = {h: i for i, h in enumerate(head)}
        for l in lines[1:]:
            t = l.split()
            rows[t[idx["SampleID"]]] = dict(
                father=t[idx["FatherID"]], mother=t[idx["MotherID"]], sex=t[idx["Sex"]],
                pop=t[idx["Superpopulation"]] if "Superpopulation" in idx else "all")
    else:
        for l in lines:
            if l.startswith("#"):
                continue
            t = l.split()
            rows[t[1]] = dict(father=t[2], mother=t[3], sex=t[4], pop="all")
    for r in rows.values():
        for k in ("father", "mother"):
            if r[k] in ("0", "", "NA"):
                r[k] = None
        r["sex"] = {"1": "male", "2": "female"}.get(r["sex"], "unknown")
    return rows


def read_bed(path):
    """Merged, sorted intervals per contig (0-based half-open, as in BED)."""
    if not path:
        return None
    iv = defaultdict(list)
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        for l in f:
            if l.startswith(("#", "track", "browser")):
                continue
            c, s, e = l.split()[:3]
            iv[c.removeprefix("chr")].append((int(s), int(e)))
    out = {}
    for c, xs in iv.items():
        xs.sort(); m = [list(xs[0])]
        for s, e in xs[1:]:
            if s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        a = np.array(m, dtype=np.int64)
        out[c] = (a[:, 0], a[:, 1])
    return out


def in_bed(bed, chrom, pos):
    """pos: 1-based positions (numpy). True where inside the target."""
    if bed is None:
        return np.ones(len(pos), bool)
    iv = bed.get(chrom.removeprefix("chr"))
    if iv is None:
        return np.zeros(len(pos), bool)
    starts, ends = iv
    k = np.searchsorted(starts, pos - 1, side="right") - 1
    ok = k >= 0
    res = np.zeros(len(pos), bool)
    res[ok] = (pos[ok] - 1) < ends[k[ok]]
    return res


def encode(g):
    """admissible's encode_gt on a cyvcf2 allele array (n, 2): any missing allele -> missing,
    two different alleles -> het (1/2 included), haploid 0 -> hom-ref, haploid alt -> hom-alt."""
    a0, a1 = g[:, 0], g[:, 1]
    hap = a1 == -2                      # cyvcf2 pads haploid calls with -2
    miss = (a0 < 0) | ((a1 < 0) & ~hap)
    out = np.where(hap | (a0 == a1), np.where(a0 == 0, HOMREF, HOMALT), HET)
    return np.where(miss, MISS, out).astype(np.int8)


def chunks(vcf_path, samples, phased=False, autosomes_only=False, chrom=None, snvs_only=False):
    """Every variant record, as admissible reads it (indels and multiallelic sites included,
    reference blocks skipped), in chunks: (chrom, pos[], pass[], codes[sites,samples], haps or None)."""
    from cyvcf2 import VCF
    v = VCF(vcf_path, samples=samples, lazy=True, threads=2)
    order = [v.samples.index(s) for s in samples]
    pos, pas, codes, haps, cur = [], [], [], [], None
    for rec in v:
        c = rec.CHROM.removeprefix("chr")
        if chrom and c != chrom:
            continue
        if autosomes_only and not c.isdigit():
            continue
        alts = [x for x in rec.ALT if x != "<NON_REF>"]
        if not alts:
            continue
        if snvs_only and (not rec.is_snp or len(alts) != 1):
            continue
        if cur is not None and c != cur and pos:
            yield cur, np.array(pos), np.array(pas), np.array(codes, np.int8), (np.array(haps, np.int8) if phased else None)
            pos, pas, codes, haps = [], [], [], []
        cur = c
        g = rec.genotype.array()[order, :2]
        pos.append(rec.POS); pas.append(rec.FILTER is None)
        codes.append(encode(g))
        if phased:
            haps.append(np.where(g < 0, -1, g))
        if len(pos) == CHUNK:
            yield cur, np.array(pos), np.array(pas), np.array(codes, np.int8), (np.array(haps, np.int8) if phased else None)
            pos, pas, codes, haps = [], [], [], []
    if pos:
        yield cur, np.array(pos), np.array(pas), np.array(codes, np.int8), (np.array(haps, np.int8) if phased else None)


def add_error(codes, rate, rng):
    """Each called genotype is replaced, with probability `rate`, by one of the two other genotypes."""
    if rate == 0:
        return codes
    out = codes.copy()
    called = out != MISS
    flip = called & (rng.random(out.shape) < rate)
    alt = np.array([[HET, HOMALT], [HOMREF, HOMALT], [0, 0], [HOMREF, HET]], np.int8)
    pick = rng.integers(0, 2, size=out.shape)
    out[flip] = alt[out[flip], pick[flip]]
    return out


# ----------------------------------------------------------------- pair metrics
def pair_counts(codes, I, J):
    """admissible's dense-pair arithmetic, for many pairs at once (autosomal sites)."""
    a, b = codes[:, I], codes[:, J]
    sa, sb = a != MISS, b != MISS
    both = sa & sb
    n_inf = both.sum(0)
    n_het = ((a == HET) & sb).sum(0) + ((b == HET) & sa).sum(0)
    n_hethet = ((a == HET) & (b == HET)).sum(0)
    n_ibs0 = (((a == HOMALT) & (b == HOMREF)) | ((a == HOMREF) & (b == HOMALT))).sum(0)
    return np.stack([n_inf, n_het, n_hethet, n_ibs0])


def define_pairs(ped, samples, n_unrelated, rng):
    s = set(samples); P = []
    for c, r in ped.items():
        for p in (r["father"], r["mother"]):
            if c in s and p in s:
                P.append(("parent-offspring", p, c))
    kids = defaultdict(list)
    for c, r in ped.items():
        if c in s and r["father"] and r["mother"]:
            kids[(r["father"], r["mother"])].append(c)
    for ks in kids.values():
        for i in range(len(ks)):
            for j in range(i + 1, len(ks)):
                P.append(("full sibling", ks[i], ks[j]))
    for c, r in ped.items():
        for p in (r["father"], r["mother"]):
            if p and p in ped and c in s:
                for g in (ped[p]["father"], ped[p]["mother"]):
                    if g in s:
                        P.append(("second-degree", g, c))
    founders = defaultdict(list)
    related = {x for t in P for x in t[1:]}
    for x in samples:
        if x in ped and not ped[x]["father"] and not ped[x]["mother"]:
            founders[ped[x]["pop"]].append(x)
    pops = [p for p, v in founders.items() if len(v) >= 2]
    seen = set()
    possible = sum(len(founders[p]) * (len(founders[p]) - 1) // 2 for p in pops)
    n_unrelated = min(n_unrelated, possible)
    tries = 0
    while pops and len(seen) < n_unrelated and tries < 50 * max(n_unrelated, 1):
        tries += 1
        p = rng.choice(pops); a, b = rng.sample(founders[p], 2)
        if (a, b) in seen or (b, a) in seen:
            continue
        if a in related and b in related and any({a, b} == {t[1], t[2]} for t in P):
            continue
        seen.add((a, b)); P.append(("unrelated", a, b))
    return P


def trios(ped, samples):
    s = set(samples)
    return [(r["father"], r["mother"], c) for c, r in ped.items()
            if c in s and r["father"] in s and r["mother"] in s]


def conditions(bed):
    regs = ["whole chromosome"] + (["exome target"] if bed is not None else [])
    return [(reg, flt, e) for reg in regs for flt in ("all calls", "PASS only") for e in ERRORS]


def cmd_pairs(a):
    rng = random.Random(a.seed); nrng = np.random.default_rng(a.seed)
    ped = read_ped(a.ped); bed = read_bed(a.bed)
    from cyvcf2 import VCF
    vs = VCF(a.vcf).samples
    pairs = define_pairs(ped, vs, a.unrelated, rng)
    used = sorted({x for t in pairs for x in t[1:]})
    ix = {x: i for i, x in enumerate(used)}
    I = np.array([ix[t[1]] for t in pairs]); J = np.array([ix[t[2]] for t in pairs])
    conds = conditions(bed)
    acc = {c: np.zeros((4, len(pairs)), np.int64) for c in conds}
    for chrom, pos, pas, codes, _ in chunks(a.vcf, used, autosomes_only=True):
        tgt = in_bed(bed, chrom, pos)
        for reg, flt, e in conds:
            m = (tgt if reg == "exome target" else np.ones_like(tgt)) & (pas if flt == "PASS only" else True)
            if m.any():
                acc[(reg, flt, e)] += pair_counts(add_error(codes[m], e, nrng), I, J)
    rows = [dict(source="observed", type=t, a=x, b=y, region=reg, filter=flt, error=e,
                 n_informative=int(v[0, k]), n_het_total=int(v[1, k]), n_ibs0=int(v[3, k]))
            for (reg, flt, e), v in acc.items() for k, (t, x, y) in enumerate(pairs)]

    if a.phased_vcf:  # simulated full sibs and parent-child pairs from phased trio parents
        T = trios(ped, VCF(a.phased_vcf).samples)
        par = sorted({x for f, m, _ in T for x in (f, m)}); pix = {x: i for i, x in enumerate(par)}
        lo_hi = {}
        sims = []  # per trio: two children, each with per-parent crossover positions and start phase
        for f, m, _ in T:
            sims.append([{p: (sorted(nrng.uniform(0, 1, nrng.poisson(a.crossovers))), int(nrng.integers(2)))
                          for p in (f, m)} for _ in range(2)])
        sacc = {c: np.zeros((4, 3 * len(T)), np.int64) for c in conds}
        for chrom, pos, pas, codes, haps in chunks(a.phased_vcf, par, phased=True, autosomes_only=True):
            if chrom not in lo_hi:
                lo_hi[chrom] = (a.chrom_start, a.chrom_end) if a.chrom_end else (0, 51_000_000)
            frac = (pos - lo_hi[chrom][0]) / (lo_hi[chrom][1] - lo_hi[chrom][0])
            kids = []
            for (f, m, _), pair in zip(T, sims):
                for child in pair:
                    al = []
                    for p in (f, m):
                        xo, start = child[p]
                        h = (start + np.searchsorted(np.array(xo), frac)) % 2
                        al.append(haps[np.arange(len(pos)), pix[p], h])
                    kids.append(encode(np.stack(al, 1)))
            K = np.stack(kids, 1)                       # sites x 2*len(T)
            Pc = codes[:, [pix[f] for f, m, _ in T]]    # observed father calls
            allc = np.concatenate([K, Pc], 1)
            n = len(T)
            I2 = np.r_[np.arange(0, 2 * n, 2), np.arange(0, 2 * n, 2)]
            J2 = np.r_[np.arange(1, 2 * n, 2), 2 * n + np.arange(n)]
            tgt = in_bed(bed, chrom, pos)
            for reg, flt, e in conds:
                msk = (tgt if reg == "exome target" else np.ones_like(tgt)) & (pas if flt == "PASS only" else True)
                if msk.any():
                    sacc[(reg, flt, e)][:, :2 * n] += pair_counts(add_error(allc[msk], e, nrng), I2, J2)
        for (reg, flt, e), v in sacc.items():
            n = len(T)
            for k in range(2 * n):
                f, m, c = T[k % n]
                rows.append(dict(source="simulated", type="full sibling" if k < n else "parent-offspring",
                                 a=f"{f}x{m}_child1", b=f"{f}x{m}_child2" if k < n else f,
                                 region=reg, filter=flt, error=e, n_informative=int(v[0, k]),
                                 n_het_total=int(v[1, k]), n_ibs0=int(v[3, k])))
    write_pairs(rows, Path(a.out), a.ibs0_max, a.min_count)


def write_pairs(rows, out, thr, min_count):
    out.mkdir(parents=True, exist_ok=True)
    for r in rows:
        r["ibs0_per_het"] = r["n_ibs0"] / r["n_het_total"] if r["n_het_total"] else float("nan")
        r["blocked_as_not_parent"] = int(r["n_ibs0"] >= min_count and r["ibs0_per_het"] > thr)
    with gzip.open(out / "pairs.csv.gz", "wt", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    grp = defaultdict(list)
    for r in rows:
        grp[(r["source"], r["type"], r["region"], r["filter"], r["error"])].append(r)
    qs = (0, 0.01, 0.05, 0.5, 0.95, 0.99, 0.999, 1)
    summ = []
    for k, rs in sorted(grp.items()):
        x = np.array([r["ibs0_per_het"] for r in rs if not math.isnan(r["ibs0_per_het"])])
        if not len(x):
            continue
        d = dict(zip(("source", "type", "region", "filter", "error"), k), n_pairs=len(x),
                 median_n_het_total=int(np.median([r["n_het_total"] for r in rs])),
                 **{f"q{q}": round(float(np.quantile(x, q)), 5) for q in qs},
                 above_current_limit=int(sum(r["blocked_as_not_parent"] for r in rs)))
        summ.append(d)
    with open(out / "pairs_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summ[0])); w.writeheader(); w.writerows(summ)
    for d in summ:
        if d["type"] in ("parent-offspring", "full sibling"):
            print(f"{d['source']:9s} {d['type']:16s} {d['region']:16s} {d['filter']:9s} err={d['error']:<5} "
                  f"n={d['n_pairs']:4d}  median={d['q0.5']:.4f}  q99.9={d['q0.999']:.4f}  max={d['q1']:.4f}  "
                  f"min={d['q0']:.4f}  over {thr}: {d['above_current_limit']}")


# ----------------------------------------------------------------- sex metrics
def wilson(k, n, z=1.96):
    if n <= 0:
        return float("nan"), float("nan")
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d
    h = (z / d) * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return max(0.0, c - h), min(1.0, c + h)


def x_usable(pos, build):
    from admissible.contigs import PAR, XTR
    keep = np.ones(len(pos), bool)
    for s, e in PAR.get(build, {}).get("X", []) + XTR.get(build, {}).get("X", []):
        keep &= ~((pos >= s) & (pos <= e))
    return keep


def cmd_sex(a):
    ped = read_ped(a.ped); bed = read_bed(a.bed)
    from cyvcf2 import VCF
    S = [s for s in VCF(a.x_vcf).samples if s in ped and ped[s]["sex"] != "unknown"]
    pops = sorted({ped[s]["pop"] for s in S}); pid = np.array([pops.index(ped[s]["pop"]) for s in S])
    founders = np.array([not ped[s]["father"] and not ped[s]["mother"] for s in S])
    obs_het = np.zeros(len(S)); exp_het = np.zeros(len(S)); a_het = np.zeros(len(S)); a_nr = np.zeros(len(S))
    for path in a.auto_vcf:   # autosomes: inbreeding F (method of moments, within superpopulation) and het fraction
        for chrom, pos, pas, codes, _ in chunks(path, S, autosomes_only=True):
            m = in_bed(bed, chrom, pos) if a.bed_autosomes else np.ones(len(pos), bool)
            c = codes[m]; called = c != MISS
            a_het += (c == HET).sum(0); a_nr += ((c == HET) | (c == HOMALT)).sum(0)
            for k in range(len(pops)):
                cols = (pid == k) & founders
                if cols.sum() < 20:
                    continue
                g = np.where(called[:, cols], np.where(c[:, cols] == HOMALT, 2, (c[:, cols] == HET).astype(int)), 0)
                n = called[:, cols].sum(1)
                p = np.divide(g.sum(1), 2 * n, out=np.zeros(len(n)), where=n > 0)
                common = (p >= 0.05) & (p <= 0.95)
                e2 = 2 * p * (1 - p)
                sel = pid == k
                cc = c[common][:, sel]; cl = called[common][:, sel]
                obs_het[sel] += (cc == HET).sum(0)
                exp_het[sel] += (e2[common][:, None] * cl).sum(0)
    F = 1 - np.divide(obs_het, exp_het, out=np.full(len(S), np.nan), where=exp_het > 0)
    x_nr = np.zeros(len(S)); x_het = np.zeros(len(S))
    for chrom, pos, pas, codes, _ in chunks(a.x_vcf, S, chrom="X"):
        m = x_usable(pos, a.build) & in_bed(bed, chrom, pos)
        c = codes[m]
        x_het += (c == HET).sum(0); x_nr += ((c == HET) | (c == HOMALT)).sum(0)
    rows = []
    for i, s in enumerate(S):
        lo, hi = wilson(int(x_het[i]), int(x_nr[i]))
        rows.append(dict(sample=s, sex=ped[s]["sex"], pop=ped[s]["pop"], F=round(float(F[i]), 4),
                         x_nonref=int(x_nr[i]), x_het=int(x_het[i]),
                         x_het_frac=round(x_het[i] / x_nr[i], 4) if x_nr[i] else float("nan"),
                         ci_low=round(lo, 4), ci_high=round(hi, 4),
                         auto_het_frac=round(a_het[i] / a_nr[i], 4) if a_nr[i] else float("nan")))
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    with open(out / "sex_samples.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    score_sex_rules(rows, out, a)


def score_sex_rules(rows, out, a):
    """Score the current rule and candidate fixes against the pedigree sex."""
    from admissible.stats import largest_gap
    ok = [r for r in rows if r["x_nonref"] >= 50]
    gap = largest_gap([r["x_het_frac"] for r in ok], 0.40, 0.58)
    boundary_cohort = (gap[0] + gap[1]) / 2 if gap and gap[1] - gap[0] >= 0.08 else 0.45
    auto_med = {p: np.nanmedian([r["auto_het_frac"] for r in rows if r["pop"] == p]) for p in {r["pop"] for r in rows}}

    def current(r, b):
        if r["x_nonref"] < 50: return "unknown"
        return "female" if r["ci_low"] > b else "male" if r["ci_high"] < b else "unknown"

    def three_way(r, b, ceiling):
        call = current(r, b)
        return "unknown" if call == "male" and r["ci_high"] >= ceiling else call

    def demote(r, b, ceiling, auto_ratio):
        call = three_way(r, b, ceiling)
        low_auto = r["auto_het_frac"] < auto_ratio * auto_med[r["pop"]]
        return "unknown (possible homozygous X)" if call == "male" and low_auto else call

    rules = [("current, boundary 0.45", lambda r: current(r, 0.45)),
             (f"current, cohort boundary {boundary_cohort:.3f}", lambda r: current(r, boundary_cohort))]
    rules += [(f"three-way, male ceiling {c}", (lambda c: lambda r: three_way(r, 0.45, c))(c)) for c in a.ceilings]
    rules += [(f"three-way {c} + autosomal demotion {q}", (lambda c, q: lambda r: demote(r, 0.45, c, q))(c, q))
              for c in a.ceilings for q in a.auto_ratios]
    res = []
    for name, rule in rules:
        for grp, sel in (("females, F >= 0.03", lambda r: r["sex"] == "female" and r["F"] >= 0.03),
                         ("females, F < 0.03", lambda r: r["sex"] == "female" and r["F"] < 0.03),
                         ("males", lambda r: r["sex"] == "male")):
            rs = [r for r in rows if sel(r)]
            calls = [rule(r) for r in rs]
            wrong = sum(c in ("male", "female") and c != r["sex"] for c, r in zip(calls, rs))
            und = sum(c.startswith("unknown") for c in calls)
            res.append(dict(rule=name, group=grp, n=len(rs), wrong=wrong, undetermined=und,
                            correct=len(rs) - wrong - und))
    with open(out / "sex_rules.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(res[0])); w.writeheader(); w.writerows(res)
    for d in res:
        print(f"{d['rule']:45s} {d['group']:20s} n={d['n']:4d} wrong={d['wrong']:3d} undetermined={d['undetermined']:3d}")


# ----------------------------------------------------------------- self-check
def cmd_selfcheck(a):
    """Re-implementation must equal admissible on CEPH 1463, pair by pair and sample by sample."""
    from admissible.checks.identity import IdentityConfig, build_sets, pair_evidence, sex_evidence
    from admissible.vcfio import load_cohort
    m = load_cohort([a.vcf]); sets = build_sets(m)
    ref = {(e.a, e.b): (e.n_het_total, e.n_ibs0) for e in pair_evidence(m, sets)}
    S = m.samples; ix = {s: i for i, s in enumerate(S)}
    I = np.array([ix[x] for x, _ in ref]); J = np.array([ix[y] for _, y in ref])
    acc = np.zeros((4, len(ref)), np.int64)
    for chrom, pos, pas, codes, _ in chunks(a.vcf, S, autosomes_only=True):
        acc += pair_counts(codes, I, J)
    bad = [(k, ref[k], (int(acc[1, n]), int(acc[3, n]))) for n, k in enumerate(ref) if ref[k] != (acc[1, n], acc[3, n])]
    print(f"pairs: {len(ref) - len(bad)}/{len(ref)} identical (n_het_total, n_ibs0)")
    for b in bad[:5]:
        print("  differs:", b)
    ev, _, _ = sex_evidence(m, sets, IdentityConfig(), a.build)
    xn = np.zeros(len(S)); xh = np.zeros(len(S))
    for chrom, pos, pas, codes, _ in chunks(a.vcf, S, chrom="X"):
        k = x_usable(pos, a.build) if a.build else np.ones(len(pos), bool)
        xh += (codes[k] == HET).sum(0); xn += ((codes[k] == HET) | (codes[k] == HOMALT)).sum(0)
    same = sum(ev[s].n_x_sites == xn[i] and ev[s].n_x_het == xh[i] for i, s in enumerate(S))
    print(f"chrX: {same}/{len(S)} samples identical (non-ref sites, het sites)")
    if bad or same != len(S):
        sys.exit("self-check FAILED: the re-implementation does not match admissible")
    print("self-check passed")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("selfcheck"); s.add_argument("--vcf", required=True); s.add_argument("--ped", required=True)
    s.add_argument("--build", default=None)
    q = sub.add_parser("pairs")
    for k in ("--vcf", "--ped", "--out"):
        q.add_argument(k, required=True)
    q.add_argument("--phased-vcf"); q.add_argument("--bed")
    q.add_argument("--unrelated", type=int, default=1000); q.add_argument("--seed", type=int, default=20261114)
    q.add_argument("--crossovers", type=float, default=0.7, help="mean crossovers per meiosis on the chromosome (chr22 ~0.7)")
    q.add_argument("--chrom-start", type=int, default=0); q.add_argument("--chrom-end", type=int, default=0)
    q.add_argument("--ibs0-max", type=float, default=0.005); q.add_argument("--min-count", type=int, default=20)
    x = sub.add_parser("sex")
    for k in ("--x-vcf", "--ped", "--out"):
        x.add_argument(k, required=True)
    x.add_argument("--auto-vcf", action="append", required=True); x.add_argument("--bed")
    x.add_argument("--bed-autosomes", action="store_true"); x.add_argument("--build", default="GRCh38")
    x.add_argument("--ceilings", type=float, nargs="+", default=[0.25, 0.30, 0.35])
    x.add_argument("--auto-ratios", type=float, nargs="+", default=[0.8, 0.9])
    a = p.parse_args()
    {"selfcheck": cmd_selfcheck, "pairs": cmd_pairs, "sex": cmd_sex}[a.cmd](a)


if __name__ == "__main__":
    main()

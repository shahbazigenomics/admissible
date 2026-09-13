"""Command line interface.

Every check runs standalone (``admissible identity ...``) and the full audit
runs as ``admissible audit ...``.  Exit codes are meant for pipelines:

    0  no blocking problem found
    1  at least one check FAILed or a blocking finding was raised
    2  the tool could not read its inputs well enough to answer

A malformed input is a finding, not a crash: it exits 1 or 2 with a report,
never with a traceback.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .checks import (
    IdentityConfig,
    check_callability,
    check_genotype,
    check_identity,
    check_models,
    check_provenance,
)
from .checks.callability import load_coverage
from .checks.models import ModelConfig
from .model import Report, Status
from .ped import read_ped
from .report import render_json, render_text
from .vcfio import load_cohort


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="admissible",
        description="Audit whether exome data can support an interpretation.",
    )
    p.add_argument("--version", action="version", version=f"admissible {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    def common(
        sp: argparse.ArgumentParser, need_ped: bool = True, need_vcf: bool = True
    ) -> None:
        # Check 4 answers a question a VCF cannot answer, so it must not demand
        # one: coverage BEDs and a target are the whole of its input.
        sp.add_argument(
            "vcf", nargs="+" if need_vcf else "*", type=Path,
            help="VCF/VCF.gz files, or ANNOVAR *_multianno.txt tables",
        )
        if need_ped:
            sp.add_argument("--ped", type=Path, required=True, help="PED file")
        sp.add_argument("--family", default=None, help="family id to label the report")
        sp.add_argument("--json", type=Path, default=None, help="write the JSON report here")
        sp.add_argument("--build", default=None, choices=["GRCh37", "GRCh38"])
        sp.add_argument(
            "--coverage", action="append", default=[], metavar="LABEL=BED",
            help="callable regions for one sample; repeat per sample. Accepts a "
                 "mosdepth --quantize output directly.",
        )
        sp.add_argument(
            "--target", type=Path, default=None,
            help="capture kit target BED - the denominator for 'fraction searched'",
        )
        sp.add_argument(
            "--af-field", default=None, metavar="KEY",
            help="name of an allele-frequency field already present in your own "
                 "annotation (INFO key or ANNOVAR column), used to rarity-filter "
                 "the model sweep. Nothing is fetched or bundled.",
        )
        sp.add_argument(
            "--max-af", type=float, default=0.01,
            help="rarity threshold applied to --af-field (default 0.01)",
        )
        sp.add_argument("-v", "--verbose", action="store_true", help="print every finding")

    common(sub.add_parser("audit", help="run every check and print the one-page verdict"))
    common(sub.add_parser("identity", help="check 1 only: sample identity"))
    common(sub.add_parser("provenance", help="check 3 only: content provenance"), need_ped=False)
    common(sub.add_parser("genotype", help="check 2 only: genotype confidence"), need_ped=False)
    common(
        sub.add_parser("callability", help="check 4 only: callable fraction"),
        need_vcf=False,
    )
    common(sub.add_parser("models", help="check 5 only: inheritance model sweep"))
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    missing = [str(p) for p in args.vcf if not Path(p).exists()]
    if missing:
        print(f"admissible: no such file: {', '.join(missing)}", file=sys.stderr)
        return 2

    matrix = load_cohort(list(args.vcf), af_key=args.af_field) if args.vcf else None
    scans = (
        list({id(s): s for s in matrix.scans.values()}.values()) if matrix else []
    )
    ped = read_ped(args.ped) if getattr(args, "ped", None) else None

    family = args.family
    if family is None and ped is not None and ped.families:
        fams = sorted(ped.families)
        family = fams[0] if len(fams) == 1 else f"{len(fams)} families: {','.join(fams)}"
    report = Report(family_id=family or "unspecified", tool_version=__version__)
    report.inputs = {
        "vcfs": [str(p) for p in args.vcf],
        "ped": str(args.ped) if getattr(args, "ped", None) else None,
        "n_samples": len(matrix.samples) if matrix else 0,
        "n_sites": matrix.n_sites if matrix else 0,
    }

    cfg = IdentityConfig()
    if args.command in ("audit", "identity"):
        report.checks.append(check_identity(matrix, ped, cfg, build=args.build))
    if args.command in ("audit", "provenance"):
        report.checks.append(check_provenance(scans, target_bed=args.target))
    if args.command in ("audit", "genotype"):
        report.checks.append(check_genotype(list(args.vcf)))
    if args.command in ("audit", "callability", "models"):
        spec: dict[str, str] = {}
        bad = [c for c in getattr(args, "coverage", []) if "=" not in c]
        if bad:
            print(f"admissible: --coverage needs LABEL=PATH, got {bad}", file=sys.stderr)
            return 2
        for item in getattr(args, "coverage", []):
            label, path = item.split("=", 1)
            spec[label] = path
        report.checks.append(
            check_callability(
                load_coverage(spec) if spec else None,
                ped,
                target_bed=args.target,
                # --family is a free-text report label; only use it to select
                # affected samples when it actually names a family in the PED.
                family_id=(
                    args.family
                    if ped and args.family in ped.families
                    else None
                ),
            )
        )
    if args.command in ("audit", "models"):
        cov_res = report.get("callability")
        report.checks.append(
            check_models(
                matrix,
                ped,
                callability=(cov_res.metrics.get("models") if cov_res else None),
                cfg=ModelConfig(
                    max_af=args.max_af if args.af_field else None
                ),
                family_id=(
                    args.family if ped and args.family in ped.families else None
                ),
            )
        )

    sys.stdout.write(render_text(report, verbose=args.verbose))
    if args.json:
        try:
            args.json.write_text(render_json(report))
        except OSError as exc:
            print(f"admissible: could not write {args.json}: {exc}", file=sys.stderr)
            return 2

    if any(c.status is Status.FAIL or c.blocking for c in report.checks):
        return 1
    if all(c.status in (Status.UNKNOWN, Status.SKIPPED) for c in report.checks):
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

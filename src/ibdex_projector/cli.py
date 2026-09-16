from __future__ import annotations

import argparse

from .project import project_counts_to_json


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="ibdex")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("project", help="project raw counts into the frozen IBDEX reference")
    p.add_argument("--counts", required=True, help="raw count matrix CSV/TSV")
    p.add_argument("--metadata", required=True, help="sample metadata CSV")
    p.add_argument("--out", required=True, help="output JSON path")
    p.add_argument("--sample-col", default="sample_id")
    p.add_argument("--tissue-col", default="tissue")
    p.add_argument("--gene-col", default=None, help="gene identifier column for genes x samples input")
    p.add_argument("--device", default="cpu")
    args = parser.parse_args(argv)

    if args.command == "project":
        out = project_counts_to_json(
            args.counts,
            args.metadata,
            sample_col=args.sample_col,
            tissue_col=args.tissue_col,
            gene_col=args.gene_col,
            device=args.device,
            out=args.out,
        )
        print(out)
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())

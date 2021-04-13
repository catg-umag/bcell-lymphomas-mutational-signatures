#!/usr/bin/env python3
import argparse
import csv
import re

from twobitreader import TwoBitFile


# AID patterns (3 nucleotides, reference C/T only), using "<R><A> <CCC>"" format
AID_PATTERNS = {
    "RCY": re.compile(r"^CT [AG].[CT]$"),
    "WA": re.compile(r"^T[CGA] [ACGT].[AT]$"),
    "RCG": re.compile(r"^CT [AG].G$"),
}


def main():
    args = parse_arguments()

    ig_locus = load_ig_locus(args.ig_locus) if args.ig_locus else None
    reference = TwoBitFile(args.reference)

    columns = [
        "sample",
        "group",
        "chrom",
        "pos",
        "ref",
        "alt",
        "context",
        "mutation_type",
        "aid_pattern",
        "ig",
    ]

    with open(args.input) as inf, open(args.output, "w") as outf:
        reader = csv.DictReader(inf)
        writer = csv.DictWriter(outf, fieldnames=columns)

        writer.writeheader()
        for row in reader:
            writer.writerow(
                annotate_variant(variant=row, reference=reference, ig_locus=ig_locus)
            )


def annotate_variant(
    variant: dict,
    reference: TwoBitFile,
    ig_locus: dict,
) -> dict:
    # format chromosome using UCSC convention
    fmt_chrom = lambda x: x if x.startswith("chr") else f"chr{x}"
    chrom = fmt_chrom(variant["chrom"])
    pos = int(variant["pos"])

    # get contextggg
    context = reference[chrom][pos - 2 : pos + 1].upper()

    # get signature "canonic" mutation (ref C/T only)
    if variant["ref"] in "AG":
        mutation = (
            f"{reverse_complement(variant['alt'] + variant['ref'])} "
            + reverse_complement(context)
        )
    else:
        mutation = f"{variant['ref']}{variant['alt']} {context}"

    # get AID pattern
    for name, pattern in AID_PATTERNS.items():
        if re.match(pattern, mutation):
            aid_pattern = name
            break
    else:
        aid_pattern = None

    # check if variant is in IG locus
    in_ig = (
        chrom in ig_locus and ig_locus[chrom][0] <= pos <= ig_locus[chrom][1]
    )

    # reformat context and mutation
    context = re.sub(r"(.).(.)", "\\1.\\2", context)
    mutation = re.sub(r"(.)(.) (.).(.)", "\\3[\\1>\\2]\\4", mutation)

    annotated_variant = {
        **variant,
        "chrom": chrom,
        "context": context,
        "mutation_type": mutation,
        "aid_pattern": aid_pattern,
        "ig": in_ig,
    }

    return annotated_variant


def load_ig_locus(filename: str) -> dict:
    locus = {}
    with open(filename) as f:
        for line in f.read().splitlines():
            if not line.startswith("#"):
                chrom, start, end, _ = line.split("\t")
                locus[chrom] = [int(start), int(end)]

    return locus


def reverse_complement(seq: str):
    """
    Basic reverse complement
    """
    complements = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complements[x] for x in seq[::-1])


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Annotates a list of variants context, C/T> mutation, AID motifs, and IG locus"
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="list containing variants to annotate",
    )
    parser.add_argument("--output", "-o", required=True, help="output CSV file")
    parser.add_argument(
        "--reference",
        "-r",
        required=True,
        help="VCF reference (to get variant contexts), in 2bit format",
    )
    parser.add_argument(
        "--ig-locus",
        "-g",
        required=True,
        help="IG locus (in bed format)",
        dest="ig_locus",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import csv
from typing import List

from cyvcf2 import VCF


def main():
    args = parse_arguments()

    # load all variants
    variants = []
    with open(args.inputs) as f:
        reader = csv.DictReader(f)
        for row in reader:
            variants += get_records(
                vcf_path=row["file"],
                sample_name=row["name"],
                group=row["group"],
                no_alt=args.no_alt,
                autosomes_only=args.autosomes_only,
            )

    # write csv
    columns = ["sample", "group", "chrom", "pos", "ref", "alt"]
    with open(args.output, "w") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(variants)


def get_records(
    vcf_path: str,
    sample_name: str,
    group: str,
    no_alt: bool = False,
    autosomes_only: bool = True,
) -> List[dict]:
    """
    Load records from VCF and create a list of dictionaries with their relevant info
    """
    hs_autosomes = [str(x) for x in range(1, 23)]
    hs_standard = hs_autosomes + ["X", "Y"]

    records = []
    for variant in VCF(vcf_path):
        chrom_ncbi = variant.CHROM.replace("chr", "")
        if (no_alt and chrom_ncbi not in hs_standard) or (
            autosomes_only and chrom_ncbi not in hs_autosomes
        ):
            continue

        record = {
            "sample": sample_name,
            "group": group,
            "chrom": variant.CHROM,
            "pos": variant.POS,
            "ref": variant.REF,
            "alt": variant.ALT[0],
        }

        records.append(record)

    return records


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Reads a set of VCF files and creates list with the variants"
    )
    parser.add_argument(
        "--inputs",
        "-i",
        required=True,
        help="list containing VCF files to be processed (CSV file with columns name,group,file)",
    )
    parser.add_argument(
        "--no-alt",
        "-n",
        action="store_true",
        help="ignore variant from alternate sequences",
        dest="no_alt",
    )
    parser.add_argument(
        "--autosomes-only",
        "-a",
        action="store_true",
        help="keep only variatns from autosomes (1-22), imples --no-alt",
        dest="autosomes_only",
    )
    parser.add_argument("--output", "-o", required=True, help="output CSV file")

    return parser.parse_args()


if __name__ == "__main__":
    main()

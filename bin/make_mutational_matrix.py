#!/usr/bin/env python3
import pandas as pd
import argparse


def main():
    args = parse_arguments()

    df = pd.read_csv(args.input)
    table = (
        df.groupby(["sample", "mutation_type"])
        .size()
        .to_frame("n")
        .pivot_table(index="mutation_type", columns="sample", values="n", fill_value=0)
        .reset_index()
        .assign(
            X=lambda x: x.mutation_type.str.replace(
                r"(.)\[(.)>(.)\](.)", "\\2\\3 \\1.\\4", regex=True
            )
        )
        .sort_values("X")
        .drop(columns=["X"])
    )
    table.to_csv(args.output, index=False, sep="\t" if args.use_tabs else ",")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Makes mutational matrix from list of variants"
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="list of variants (CSV)",
    )
    parser.add_argument("--output", "-o", required=True, help="output CSV file")
    parser.add_argument(
        "--tab-delimiter",
        "-t",
        action="store_true",
        help="use tab as separator instead of comma",
        dest="use_tabs",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()

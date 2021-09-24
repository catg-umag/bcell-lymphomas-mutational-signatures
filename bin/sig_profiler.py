#!/usr/bin/env python3
import argparse
import re
from os import listdir, makedirs

import pandas as pd
from SigProfilerExtractor import sigpro as sig


def main():
    args = parse_arguments()

    # run sigprofiler
    sig.sigProfilerExtractor(
        input_type="matrix",
        output="sigprofiler_out",
        input_data=args.input,
        minimum_signatures=args.min_signatures,
        maximum_signatures=args.max_signatures,
        cpu=args.ncpus,
        gpu=args.use_gpu,
    )

    signature_dir = (
        f"All_Solutions/SBS96_{args.forced_signatures}_Signatures"
        if args.forced_signatures
        else "Suggested_Solution/SBS96_De-Novo_Solution"
    )
    base_dir = f"sigprofiler_out/SBS96/{signature_dir}/Signatures"

    makedirs(args.output_dir, exist_ok=True)

    # copy signatures file
    signatures_file = next(
        x for x in listdir(base_dir) if re.match(r"SBS96_.*_Signatures.txt", x)
    )
    signatures_df = pd.read_csv(f"{base_dir}/{signatures_file}", sep="\t")
    signatures_df.columns = ["mutation"] + [
        f"S{i}" for i in range(1, len(signatures_df.columns))
    ]
    signatures_df.to_csv(f"{args.output_dir}/signatures.csv", index=False)

    # copy solution stats
    stats_df = pd.read_csv("sigprofiler_out/SBS96/All_solutions_stat.csv")
    stats_df = stats_df[
        ["Signatures", "Stability (Avg Silhouette)", "Mean Cosine Distance"]
    ]
    stats_df.columns = ["signatures", "stability", "mean_cosine_distance"]
    stats_df.to_csv(f"{args.output_dir}/statistics.csv", index=False)

    # copy contributions
    contrib_base_dir = f"sigprofiler_out/SBS96/{signature_dir}/Activities"
    contrib_file = next(
        x for x in listdir(contrib_base_dir) if re.match(r"SBS96_.*_Activities.txt", x)
    )
    contrib_df = pd.read_csv(f"{contrib_base_dir}/{contrib_file}", sep="\t")
    contrib_df.columns = ["sample"] + [
        f"S{i}" for i in range(1, len(contrib_df.columns))
    ]
    contrib_df.to_csv(f"{args.output_dir}/contributions.csv", index=False)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Executes signature profiler from matrix and extract relevant info"
    )
    parser.add_argument("--input", "-i", required=True, help="Input matrix (.csv)")
    parser.add_argument(
        "--output-dir",
        "-o",
        default="output",
        help="Output directory",
        dest="output_dir",
    )
    parser.add_argument(
        "--min-signatures",
        type=int,
        default=1,
        help="Minimum number of signatures to extract",
        dest="min_signatures",
    )
    parser.add_argument(
        "--max-signatures",
        type=int,
        required=True,
        help="Maximum number of signatures to extract",
        dest="max_signatures",
    )
    parser.add_argument(
        "--force-nsignatures",
        "-f",
        type=int,
        help="Force the selection of a given number of signatures",
        dest="forced_signatures",
    )
    parser.add_argument(
        "--n-cpus",
        "-c",
        type=int,
        default=1,
        help="Number of CPUs to use",
        dest="ncpus",
    )
    parser.add_argument(
        "--use-gpu",
        "-g",
        action="store_true",
        help="Use GPU for signature extraction",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

import argparse
import sys

"""
EXAMPLE:

in_string: (fixed_seq_1)(cell_barcode_1)(umi)(cell_barcode_2)(fixed_seq_2)(sequence)

in_string_config:

fixed_seqs: GATTACA, ACATTAG
barcode_lengths: 4, 3
umi_lengths: 4

"""


def get_args():
    """Parse the commandine arguments"""
    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument("-i", "--in_string", required=True)
    arg_parser.add_argument(
        "-c",
        "--cb_lengths",
        help="The comma delimited "
        "list of cell barcodes length. The position in "
        "the list corresponds to its order in in_string"
    )
    arg_parser.add_argument(
        "-u",
        "--umi_lengths",
        help="The comma delimited "
        "list of umi lengths. The position in the list "
        "corresponds to its order in in_string"
    )
    arg_parser.add_argument(
        "-f",
        "--fixed_seqs",
        help="The comma delimited "
        "list of fixed sequences. The position in the "
        "corresponds to its order in in_string"
    )
    arg_parser.add_argument(
        "-o", "--out_file", help="The out file for the regex and umi_tools pattern to be output."
    )

    args = arg_parser.parse_args()

    return args


def convert_regex(in_string, cb_lengths, umi_lengths, fixed_seqs, out_file):
    """Converts the human_readable string to proper regex"""

    feature_annos = convert_seq_features(cb_lengths, umi_lengths, fixed_seqs)

    regex = convert_to_regex(in_string, feature_annos)
    umi_tools_pattern = get_umi_tools_pattern(feature_annos["cell_barcode"], feature_annos["umi"])

    with open(out_file, "w") as f:
        f.write("REGEX" + "\t" + regex + "\n")
        f.write("UMI_TOOLS" + "\t" + umi_tools_pattern + "\n")
        f.write("BC_LENGTH" + "\t" + str(umi_tools_pattern.count("C")) + "\n")
        f.write("UMI_LENGTH" + "\t" + str(umi_tools_pattern.count("N")) + "\n")

    return


def get_umi_tools_pattern(cell_barcode_info, umi_info):
    umi_tools_pattern = ""

    for length in cell_barcode_info:
        umi_tools_pattern += "C" * int(length)

    for length in umi_info:
        umi_tools_pattern += "N" * int(length)

    return umi_tools_pattern


def convert_seq_features(cb_lengths, umi_lengths, fixed_seqs):
    seq_features = {}

    seq_features["cell_barcode"] = cb_lengths.strip().split(",") if cb_lengths else ""
    seq_features["umi"] = umi_lengths.strip().split(",") if umi_lengths else ""
    seq_features["fixed_seq"] = fixed_seqs.strip().split(",") if fixed_seqs else ""

    return seq_features


def convert_to_regex(in_string, feature_annos):
    """Converts the in_string to regex"""
    regex = ""

    for feature in in_string.split(","):
        if feature:
            feature = feature.strip()
            feature_name = "_".join(feature.split("_")[:-1])
            feature_idx = int(feature.split("_")[-1])

            if feature_name in feature_annos:
                regex += add_feature(feature_name, feature_idx, feature_annos[feature_name])
            else:
                print(f"Unknown feature: {in_string}")
                sys.exit()

    return regex


def add_feature(feature_name, feature_idx, feature_info):
    regex = "(?P<{}>{}{{{}}})"

    string_char = ""
    string_extra = ""

    if feature_name == "fixed_seq":
        string_char = "(" + feature_info[feature_idx - 1].strip() + ")"
        string_extra = "e<=3"
        feature_name = "discard"
    else:
        string_char = "."
        string_extra = feature_info[feature_idx - 1].strip()

    return regex.format("_".join([feature_name, str(feature_idx)]), string_char, string_extra)


def main():
    """Main Subroutine"""

    args = get_args()

    convert_regex(args.in_string, args.cb_lengths, args.umi_lengths, args.fixed_seqs, args.out_file)


if __name__ == "__main__":
    main()

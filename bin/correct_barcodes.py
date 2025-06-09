#!/usr/bin/env python3

""" This will take in a bam and will check barcodes that have been added by
    umi_tools and will check them against a whitelist barcode and correct
    those barcodes that are not on the list
"""

import argparse
import itertools
import math
import pysam as ps
import pygtrie
from Bio import SeqIO

QUAL_OFFSET = 33


def parse_args():
    """Parse the commandline arguments"""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--infile",
        default=None,
        type=str,
        required=True,
        help="The input csv containing read-id and estimated barcode")

    parser.add_argument(
        "-o",
        "--outfile",
        default=None,
        type=str,
        required=True,
        help="The output csv containing the corrected barcodes.")

    parser.add_argument(
        "-w",
        "--whitelist",
        default=None,
        type=str,
        required=True,
        help="The whitelist file containing cell barcodes"
    )

    parser.add_argument(
        "-b",
        "--barcode_count",
        default=None,
        type=str,
        required=True,
        help="The barcode count file, contains the union between the umi_tools generated "
        "whitelist and the proper barcode whitelist."
    )

    parser.add_argument(
        "--max_edit_dist",
        default=2,
        type=int,
        required=False,
        help="The maximum edit distance that a barcode can be from a barcode on the whitelist"
    )

    parser.add_argument(
        "--min_post_prob",
        default=0.975,
        type=int,
        required=False,
        help="The minimum posterior probability a barcode on the whitelist must have to replace "
        "the barcode detected by the pipeline"
    )

    #parser.add_argument(print_header, default=True)
    parser.set_defaults(print_header=True)

    parser.add_argument(
        "--skip_header",
        action='store_false',
        dest='print_header',
        help="Indicates whether to skip printing the header"
    )

    args = parser.parse_args()

    return args


def get_filtered_outfile(outfile):
    """Get the name of the filtered bam file

    Args:
        outfile (str): The path to the file that information will be output to

    Returns:
        filtered_outfile (str): The path to the file that filtered information
            will be output to
    """

    filtered_outfile = "filtered.tmp"

    if outfile:
        filtered_outfile = outfile + ".filtered"

    return filtered_outfile


def correct_barcode(infile, outfile, whitelist, barcode_count_file, min_post_prob, max_edit_dist, print_header):
    """Iterate through each read in the bam, determine the most likely barcode
        for each read and output it.

    Args:
        infile (str): The path to the input bam
        outfile (str): The path to the output bam
        whitelist (str): The path to the barcode whitelist file
        barcode_count_file (str): The path to the file containing each detected
            barcode and the counts for that barcode
        min_post_prob (float): The minimum probability a barcode must be to be
            considered a 'correct' barcode
        max_edit_dist (int): The maximum edit distance the detected barcode can
            be from a barcode on the whitelist in order to be corrected
        print_header (bool): Indicates whether to print the header to the output
            file or not

    Output: None

    """

    with open(infile, 'r') as infile_h, open(outfile, 'w') as outfile_h:

        # Turn the whitelist into a trie, as it allows for very quick access to check for a barcode
        whitelist_trie = read_whitelist(whitelist)

        # Turn the barcode abundances into percentages
        bc_probabilities = calculate_bc_ratios(barcode_count_file)

        header = ""

        for inline in infile_h.readlines():
            inline = inline.strip('\n')
            if not header:
                header = inline
                if print_header:
                    outfile_h.write(header + "\tcorrected_bc\n")
                continue

            read_id, bc, bc_qual, _, _ = inline.split('\t')

            if not bc:
                continue
            else:
                corrected_bc = get_read_bc(
                    bc, bc_qual, whitelist_trie, bc_probabilities, max_edit_dist, min_post_prob
                )
                if corrected_bc:
                    outfile_h.write(inline + '\t' + corrected_bc + '\n')

def read_whitelist(in_whitelist):
    """Read the barcode whitelist into a trie. While this does have a hit on
        memory, this allows for faster lookup during the barcode correction
        process.

    Args:
        in_whitelist (str): The file containing the whitelisted barcodes

    Returns:
        whitelist_trie (trie): A trie containing the barcodes in the whitelist

    """

    whitelist_trie = pygtrie.CharTrie()

    # TODO: Put in try-catch
    with open(in_whitelist, "r", encoding="UTF-8") as f_handle:
        # TODO: Can this be a one-liner? Is it worth reducing down?
        for line in f_handle.readlines():
            # The data stored on the leaf is irrelevant
            whitelist_trie[line.strip("\n")] = True

    return whitelist_trie


def calculate_bc_ratios(barcode_count_file):
    """Read the file that contains the barcodes and their amounts, and convert
        the counts into percentages

    Args:
        barcode_count_file (str): The path to the file containing barcode and
            their amounts

    Returns:
        bc_ratios (dict): The dictionary containing the ratio of barcodes as a
            percentage
    """

    # TODO: Reconsider the variable names
    bc_ratios = {}

    raw_bc_counts = {}
    bc_totals = 0

    with open(barcode_count_file, "r", encoding="UTF-8") as f_handle:
        for line in f_handle.readlines():
            cell_bc, bc_count = line.split(",")

            raw_bc_counts[cell_bc] = int(bc_count)
            bc_totals += int(bc_count)

    for cell_bc, bc_count in raw_bc_counts.items():
        # bc_ratios[cell_bc] = bc_count / bc_totals
        bc_ratios[cell_bc] = bc_count + 1

    return bc_ratios


def get_read_bc(barcode, barcode_qual, whitelist_trie, bc_probabilities, max_edit_dist, min_post_prob):
    """This will determine the 'correct' barcode for a given read

    Args:
        barcode (str): The estimated barcode for a read
        barcode_qual (str): The quality of barcode
        whitelist_trie (trie): The barcode whitelist
        bc_probabilities (dict): THe dictionary containing the barcode as the
            key with the value being how much of the total amount of barcodes
            the key composes
        min_post_prob (float): The minimum probability a barcode must be to be
            considered a 'correct' barcode
        max_edit_dist (int): The maximum edit distance the detected barcode can
            be from a barcode on the whitelist in order to be corrected

    Returns:
        barcode (str): The corrected barcode for the read

    """

    corrected_barcode = ""

    if not whitelist_trie.has_key(barcode):
        similar_bcs = get_similar_bcs(barcode, barcode_qual, whitelist_trie, bc_probabilities, max_edit_dist)
        if similar_bcs:
            corrected_barcode = find_correct_barcode(similar_bcs, min_post_prob)
    else:
        corrected_barcode = barcode

    return corrected_barcode


def get_similar_bcs(query_bc, query_bc_qual, bc_trie, bc_probabilities, max_edit_dist):
    """Return all bcs in bc_list that are within edit_distance

    Args:
        query_bc (str): The barcode that we are searching for
        query_bc_qual (str): The quality of query_bc
        bc_trie (trie): The barcode whitelist
        bc_probabilities (dict): THe dictionary containing the barcode as the
            key with the value being how much of the total amount of barcodes
            the key composes
        max_edit_dist (int): The maximum edit distance the detected barcode can
            be from a barcode on the whitelist in order to be corrected

    Return:
        similar_bcs (dict): The barcodes that are the same or are within the
            acceptable edit distance. Each are kept under separate keys
    """

    similar_bcs = []

    for cell_bc in list(get_mutated_bcs(query_bc, max_edit_dist)):
        if bc_trie.has_key(cell_bc):
            bc_probability = get_bc_probability(query_bc, query_bc_qual, cell_bc, bc_probabilities)
            print("INFO: BARCODE_CHECK - ({} - {}) Probability {}".format(query_bc, cell_bc, bc_probability))

            if bc_probability > 0:
                similar_bcs.append((cell_bc, bc_probability))

    return similar_bcs


def get_mutated_bcs(query_bc, max_edit_dist):
    """Generate all strings that have up to max_edit_dist mutations in them,
        allowing us to determine if the barcode needs to be error corrected

    Args:
        query_bc (str): The barcode to be mutated
        max_edit_dist (int): The maximum amount of mutations to be applied to
            the barcode

    Returns: None, yields the mutated barcodes
    """

    for edit_dist in range(1, max_edit_dist + 1):
        for locs in itertools.combinations(range(len(query_bc)), edit_dist):
            query_bc_list = [[base] for base in query_bc]

            for loc in locs:
                orig_base = query_bc_list[loc][0]
                query_bc_list[loc] = [b for b in "ACGT" if b != orig_base]

            for poss in set(itertools.product(*query_bc_list)):
                yield "".join(poss)


def get_bc_probability(query_bc, query_bc_qual, potential_bc, bc_probabilities):
    """Calculates the probability that the changed bases between query_bc and
        potential_bc are a result of sequencing error

    Args:
        query_bc (str): The barcode that we are searching for
        query_bc_qual (str): The quality of query_bc
        potential_bc (str): The whitelisted barcode that could be query_bc's
            origin.
        bc_probabilities (dict): THe dictionary containing the barcode as the
            key with the value being how much of the total amount of barcodes
            the key composes

    Return:
        likelihood (float): The probability that query_bc originates from
            potential_bc.

    """
    tot_edit_probability = 1
    for idx in get_mismatch_locs(query_bc, potential_bc):
        edit_probability = get_edit_probability(query_bc_qual[idx])
        tot_edit_probability *= edit_probability

    potential_bc_prob = 0
    if potential_bc in bc_probabilities:
        potential_bc_prob = bc_probabilities[potential_bc]

    print("TOTAL_EDIT_PROBABILITY: {}".format(tot_edit_probability))
    print("POTENTIAL_BC_PROB: {}".format(potential_bc_prob))
    return tot_edit_probability * potential_bc_prob


def get_mismatch_locs(query_bc, potential_bc):
    """Get the location of bases that are mismatched between the original
        and potential barcodes

    Args:
        query_bc (str): The barcode that we are searching for
        potential_bc (str): The whitelisted barcode that could be query_bc's
            origin.

    Returns:
        mismatch_locs (list): The list of locations where query_bc and
            potential_bc differ. The length of this list is <= edit_distance
            that was passed to this program
    """
    mismatch_locs = []
    for index, bc_chars in enumerate(zip(query_bc, potential_bc)):
        if bc_chars[0] != bc_chars[1]:
            mismatch_locs.append(index)
    return mismatch_locs


def get_edit_probability(q_score):
    """Calculates the probability that the position is edited based on q score

    Args:
        q_score (str): The q_score of the base being evaluated

    Returns:
        edit_probability (float): The probability of how likely the base
            reported is the "correct" base
    """

    # 33 in formula below is linked  standard (e.g.: Illumina) qual offset
    # also applicable to nanopore, but reminder to continue to check
    return math.pow(10, -1 * (ord(q_score) - QUAL_OFFSET) / 10)


def find_correct_barcode(potential_bcs, min_prob):
    """Determines which barcode is the most likely origin barcode

    Args:
        potential_bcs (dict): The list of whitelisted barcodes and their
            qualities that could be the origin barcode
        min_prob (float): The minimum posterior probability that a barcode
            should be before considering it to be an option for correction

    Returns:
        max_likelihood_bc (str): The barcode deemed to be the best option
            to correct the source barcode to

    """
    # TODO: Add Docstring
    max_likelihood = 0
    max_likelihood_bc = ""

    total_probability = sum([p for s, p in potential_bcs])

    for potential_bc_info in potential_bcs:
        potential_bc, potential_bc_qual = potential_bc_info
        likelihood = potential_bc_qual / total_probability
        print("INFO: BARCODE_CHECK - ({}) Likelihood {}".format(potential_bc, likelihood))

        if likelihood > max_likelihood and likelihood > min_prob:
            max_likelihood = likelihood
            max_likelihood_bc = potential_bc

    return max_likelihood_bc


def main():
    """Main Subroutine"""
    args = parse_args()
    correct_barcode(
        args.infile,
        args.outfile,
        args.whitelist,
        args.barcode_count,
        args.min_post_prob,
        args.max_edit_dist,
        args.print_header
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

""" This will take in a bam and will check barcodes that have been added by
    umi_tools and will check them against a whitelist barcode and correct
    those barcodes that are not on the list
"""

import argparse
import gzip
import math
import editdistance
import pysam
import pygtrie
import sys
from Bio import SeqIO

# TODO: Add logging statements in various places

MAX_EDIT_DIST = 1
MIN_PROB = 0.975
QUAL_OFFSET = 33

def parse_args():
    """ Parse the commandline arguments """

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--infile", default=None, type=str,
                        required=True, help="The input bam file")

    parser.add_argument("-o", "--outfile", default=None, type=str,
                        required=True, help="The output bam file")

    parser.add_argument("-w", "--whitelist", default=None, type=str,
                        required=True, help="The whitelist file containing " \
                                            "cell barcodes")

    parser.add_argument("-b", "--barcode_count", default=None, type=str,
                        required=True, help="The barcode count file, " \
                                            "contains the union between the " \
                                            "umi_tools generated whitelist " \
                                            "and the proper barcode whitelist.")

    parser.add_argument("--max_edit_dist", default=MAX_EDIT_DIST, type=int,
                        required=False, help="The maximum edit distance that " \
                                             "a barcode can be from a " \
                                             "barcode on the whitelist")

    parser.add_argument("--min_post_prob", default=MIN_PROB, type=int,
                        required=False, help="The minimum posterior " \
                                             "probability a barcode on the " \
                                             "whitelist must have to replace " \
                                             "the barcode detected by the " \
                                             "pipeline")

    parser.add_argument("--write_filtered_reads", required=False, type=str,
                        default="", help="Output the filtered reads into a " \
                                         "the file name specified by this " \
                                         "parameter.")

    args = parser.parse_args()

    return args

def get_filtered_outfile(outfile):
    """ Gets the name of the filtered bam file """

    filtered_outfile = "filtered.tmp"

    if outfile:
        filtered_outfile = outfile

    return filtered_outfile

def tag_barcodes(infile, outfile, whitelist, barcode_count_file, write_filtered_reads):
    """ Will correct and tag reads """
    # TODO: Refactor this function
    # Open up all the necessary files

    with pysam.AlignmentFile(infile, "rb") as bam_in, \
        pysam.AlignmentFile(outfile, "wb", template=bam_in) as bam_unfilt_out, \
        pysam.AlignmentFile(get_filtered_outfile(write_filtered_reads), "wb", template=bam_in) as bam_filt_out:
        # Get the list of acceptable barcodes
        whitelist_trie = read_lines(whitelist)

        # Get the barcode count file
        bc_probabilities = read_bc_probabilities(barcode_count_file)

        for read in bam_in:
            # We assume that the Cell Barcode information is tagged with the
            #   read. Could add in an option to support UMI-tools but would
            #   still need to get the quality by parsing the fastq
            cell_bc = read.get_tag('CR')
            cell_bc_qual = read.get_tag('CY')

            # This gets all barcodes with a Levenstein distance <= MAX_EDIT_DIST
            similar_bcs, total_likelihood = get_similar_bcs(cell_bc, cell_bc_qual, whitelist_trie, bc_probabilities)

            # Because downstream umitools uses the sequence name, we need to
            #   adjust it to contain the new barcode
            seq_name = [read.query_name.split('_')[0]]

            needs_filtering = False

            if similar_bcs['same']:
                # If the barcode matches one in the list, we go ahead and tag
                # the read and continue in the loop
                read.tags += [('CB', cell_bc)]
                read.tags += [('NT', 'EXACT_MATCH')]
                
                seq_name.append(cell_bc)

            elif similar_bcs['similar']:
                # If there are inexact matches, let's try and determine which
                #    one makes the most sense to have been the original barcode
                corrected_bc = correct_barcode(similar_bcs['similar'], total_likelihood)

                # If we found an inexact match, let's update the tags
                # Otherwise we are going to filter the read out
                if corrected_bc:
                    read.tags += [('CB', corrected_bc)]
                    read.tags += [('NT', 'PARTIAL_MATCH')]

                    seq_name.append(corrected_bc)
                else:
                    # We were unable to find a partial match
                    needs_filtering = True
                    read.tags += [('NT', 'UNABLE_TO_CALC_PARTIAL_MATCH')]
            else:
                # There were no barcodes that are within the hamming distance
                needs_filtering = True
                read.tags += [('NT', 'NO_BCS_IN_EDIT_DIST')]

            if not needs_filtering:
                # Lets fix the read name to contain the correct barcode
                seq_name.append(read.query_name.split('_')[-1])
                read.query_name = '_'.join(seq_name)

                bam_unfilt_out.write(read)

            elif write_filtered_reads:
                bam_filt_out.write(read)

def read_lines(file_to_read):
    """ Will read a file and return the lines as a list """

    #file_lines = []

    whitelist_trie = pygtrie.CharTrie()

    with open(file_to_read, 'r', encoding = 'UTF-8') as f_handle:
        for line in f_handle.readlines():
            #file_lines.append(line)
            whitelist_trie[line.strip('\n')] = True

    return whitelist_trie

def read_bc_probabilities(barcode_count_file):
    """ Parse the data from the barcode_count_file into a dictionary """

    bc_probabilities = {}

    raw_bc_counts = {}
    total_counts = 0

    with open(barcode_count_file, 'r', encoding = 'UTF-8') as f_handle:
        for line in f_handle.readlines():
            cell_bc, _, bc_count, _ = line.split('\t')

            raw_bc_counts[cell_bc] = int(bc_count)
            total_counts += int(bc_count)

    for cell_bc, bc_count in raw_bc_counts.items():
        bc_probabilities[cell_bc] = bc_count / total_counts

    return bc_probabilities

def read_bc_umi_quals(bc_umi_qual_file):
    """ Will parse the fastq containing the bc_umi sequences and get their
        qualities
    """

    bc_umi_quals = {}

    with gzip.open(bc_umi_qual_file, "rt") as unzip_file:
        for record in SeqIO.parse(unzip_file, "fastq"):

            seq_name = record.id.split(' ')[0]
            phred_quals = record.letter_annotations['phred_quality']

            bc_umi_quals[seq_name] = [record.seq, phred_quals]

    return bc_umi_quals

def parse_bc_umi(read_name):
    """ Will parse the barcode and umi from the read name """

    cell_bc = read_name.split('_')[-2]
    umi = read_name.split('_')[-1]

    return cell_bc, umi

def get_similar_bcs(query_bc, query_bc_qual, bc_trie, bc_probabilities):
    """ Return all bcs in bc_list that are either the same or similar

    Args:
        query_bc (str): The barcode that we are searching for
        bc_list (str): The list of barcodes

    Return:
        similar_bcs (dict): The barcodes that are the same or are within the
            acceptable edit distance. Each are kept under separate keys
    """

    # TODO: Update the docstring

    similar_bcs = {}
    similar_bcs['same'] = []
    similar_bcs['similar'] = []

    total_likelihood = 0

    if bc_trie.has_key(query_bc):
        similar_bcs['same'] = [query_bc]
    else:
        all_bcs = get_all_bcs(query_bc)

        for cell_bc in all_bcs:
            if bc_trie.has_key(cell_bc):
                bc_probability = get_bc_probability(query_bc, query_bc_qual, cell_bc, bc_probabilities)
                total_likelihood += bc_probability

                similar_bcs['similar'].append([cell_bc, bc_probability])

    #for cell_bc in bc_list:
    #    cell_bc = cell_bc.strip()
    #    edit_dist = editdistance.eval(query_bc, cell_bc)

    #    if edit_dist == 0:
    #        similar_bcs = {}
    #        similar_bcs['same'] = [cell_bc]
    #        break

    #    elif edit_dist <= MAX_EDIT_DIST:
    #        bc_probability = get_bc_probability(query_bc, query_bc_qual, cell_bc, bc_probabilities)
    #        total_likelihood += bc_probability

    #        similar_bcs['similar'].append([cell_bc, bc_probability])

    return similar_bcs, total_likelihood

def get_all_bcs(query_bc):
    dna_bases = ['A', 'C', 'G', 'T']

    all_bcs = []

    for i in range(0, len(query_bc)):
        changed_bc = list(query_bc)

        for base in dna_bases:
            if base != query_bc[i]:
                changed_bc[i] = base
                all_bcs.append(''.join(changed_bc))

    return all_bcs


def get_bc_probability(query_bc, query_bc_qual, potential_bc, bc_probabilities):
    likelihood = 0

    # TODO: Add Docstring
    for idx in get_mismatch_locs(query_bc, potential_bc):
        edit_probability = get_edit_probability(query_bc_qual[idx])

        potential_bc_prob = 0
        if potential_bc in bc_probabilities:
            potential_bc_prob = bc_probabilities[potential_bc]

        likelihood = (potential_bc_prob + bc_probabilities[query_bc]) * edit_probability

    return likelihood

def correct_barcode(potential_bcs, total_likelihood):
    # TODO: Add Docstring
    max_likelihood = 0
    max_likelihood_bc = ""

    for potential_bc_info in potential_bcs:
        potential_bc, potential_bc_qual = potential_bc_info
        likelihood = potential_bc_qual / total_likelihood

        if likelihood > max_likelihood and likelihood > MIN_PROB:
            max_likelihood = likelihood
            max_likelihood_bc = potential_bc

    return potential_bc

def get_mismatch_locs(original_bc, potential_bc):
    """ Get the location of bases that are mismatched between the original
        and potential barcodes
    """
    mismatch_locs = []
    for index, bc_chars in enumerate(zip(original_bc,
                                         potential_bc)):
        if bc_chars[0] != bc_chars[1]:
            mismatch_locs.append(index)
    return mismatch_locs

# 33 in formula below is linked  standard (e.g.: Illumina) qual offset
# also applicable to nanopore, but reminder to continue to check
def get_edit_probability(q_score):
    """ Calculates the probability that the position is edited based on q score
    """
    return math.pow(10, -1 * (ord(q_score) - QUAL_OFFSET) / 10)

def main():
    """ Main Subroutine """
    args = parse_args()
    tag_barcodes(args.infile, args.outfile, args.whitelist,
                 args.barcode_count, args.write_filtered_reads)

if __name__ == '__main__':
    main()

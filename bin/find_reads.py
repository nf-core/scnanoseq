"""
This is modified script from the `support-5prim-kit` branch of the BLAZE repo, link below:
https://github.com/shimlab/BLAZE/tree/support-5prim-kit

The only modifications that have been made are formatting changes in order to meet the
nf-core linting requirements. The script has been copied into the pipeline temporarily
in order to support 5-prime data, and will be removed once this feature has been
incorporated into an official BLAZE release.
"""

# find specific reads with given read id and output to a new fastq file


import argparse
import Bio.SeqIO
import multiprocessing as mp
import textwrap
from pathlib import Path


import helper

def parse_arg():
    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
        '''
        Find specific reads with given read id and output to a new fastq file
        '''))

    # Required positional argument
    parser.add_argument('input_fastq_dir', type=str,
                        help='Fastq directory, Note that this should be a folder.')

    # required name argment
    requiredNamed = parser.add_argument_group('These Arguments are required')
    requiredNamed.add_argument(
        '--output_file', type=str, required = True,
        help='Filename for the output fastq.')
    requiredNamed.add_argument('--id_file', type=str, required = True,
                        help='A file containing all the read ids to look for.')

    parser.add_argument('--threads', type=int,
                        help='Number of threads. Default: # of CPU - 1')


    args = parser.parse_args()

    return args

def find_reads(in_fastq, id_list):
    fastq = Bio.SeqIO.parse(in_fastq, "fastq")
    read_list = [r for r in fastq if r.id in id_list]
    return read_list
    # check id
def main(args):

    # get ids (from args.id_file)
    ids = []
    with open (args.id_file, 'r') as f:
        for line in f:
            ids.append(line.strip())


    fastq_fns = list(Path(args.input_fastq_dir).rglob('*.fastq'))
    rst_futures = helper.multiprocessing_submit(find_reads,
                                                fastq_fns,
                                                n_process=args.threads,
                                                id_list = ids)

    rst_ls = []
    for f in rst_futures:
        rst_ls+=f.result()

    print(len(rst_ls))
    Bio.SeqIO.write(rst_ls, args.output_file, 'fastq')

if __name__ == '__main__':
    args = parse_arg()
    main(args)


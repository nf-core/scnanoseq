#!/usr/bin/env python3

import argparse
import gzip
import shutil
import os
from pathlib import Path

def cat_files(input_files: list[str], output_file: str) -> None:
    """Concatenate gzipped files."""
    with open(output_file, 'wb') as f_out:
        for f_path in input_files:
            try:
                with gzip.open(f_path, 'rb') as f_in:
                    shutil.copyfileobj(f_in, f_out)
            except gzip.BadGzipFile:
                print(f"Warning: {f_path} is not a valid gzip file. Attempting to read as plain text.")
                try:
                    with open(f_path, 'rb') as f_in_plain: # Read as binary for consistency
                        shutil.copyfileobj(f_in_plain, f_out)
                except Exception as e:
                     print(f"Error processing file {f_path}: {e}")
            except Exception as e:
                print(f"Error processing file {f_path}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Concatenate FASTQ files.")
    parser.add_argument("--prefix", type=str, required=True, help="Output file prefix.")
    parser.add_argument("--single_end", action="store_true", help="Input files are single-end.")
    parser.add_argument("--reads", nargs='+', required=True, help="List of input FASTQ files.")

    args = parser.parse_args()

    output_dir = Path(".")

    read_paths = [Path(f) for f in args.reads]

    if args.single_end:
        output_file = output_dir / f"{args.prefix}.merged.fastq.gz"
        if len(read_paths) == 1:
            print(f"Symlinking {read_paths[0]} to {output_file}...")
            os.symlink(read_paths[0], output_file)
            print("Symlink complete.")
        elif len(read_paths) > 1:
            print(f"Concatenating {len(read_paths)} single-end files to {output_file}...")
            cat_files([str(p) for p in read_paths], str(output_file))
            print("Concatenation complete.")
        else:
            print("Warning: No input files provided for single-end processing.")
    else: # Paired-end
        output_file_1 = output_dir / f"{args.prefix}_1.merged.fastq.gz"
        output_file_2 = output_dir / f"{args.prefix}_2.merged.fastq.gz"

        if len(read_paths) == 2:
            print(f"Symlinking {read_paths[0]} to {output_file_1}...")
            os.symlink(read_paths[0], output_file_1)
            print("R1 symlink complete.")
            print(f"Symlinking {read_paths[1]} to {output_file_2}...")
            os.symlink(read_paths[1], output_file_2)
            print("R2 symlink complete.")
        elif len(read_paths) > 2:
            if len(read_paths) % 2 != 0:
                print("Error: Paired-end reads require an even number of files.")
                return # Or raise error

            read1_paths = [str(read_paths[i]) for i in range(0, len(read_paths), 2)]
            read2_paths = [str(read_paths[i]) for i in range(1, len(read_paths), 2)]

            print(f"Concatenating {len(read1_paths)} R1 files to {output_file_1}...")
            cat_files(read1_paths, str(output_file_1))
            print("R1 concatenation complete.")

            print(f"Concatenating {len(read2_paths)} R2 files to {output_file_2}...")
            cat_files(read2_paths, str(output_file_2))
            print("R2 concatenation complete.")
        else: # len(read_paths) < 2
            print("Warning: Less than 2 input files provided for paired-end processing.")


if __name__ == "__main__":
    main() 
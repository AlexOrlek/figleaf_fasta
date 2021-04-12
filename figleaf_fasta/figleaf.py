#! /usr/bin/env python3
"""Main module"""
# Usage:
# import module in python script... from figleaf-fasta import figleaf
# run script from the command line... run figleaf.py -h to see help on command-line options


def figleaf(fasta_input, fasta_output, ranges_list=None, task='hard_mask', hard_mask_letter='N', ranges_path=None, inverse_mask=False):
    """Mask or extract FASTA sequences according to ranges specified by ranges_list or ranges_path (0-indexed start, end positions with non-inclusive end position

    Args:
      fasta_input (str): Input FASTA filepath
      fasta_output (str): Output FASTA filepath
      ranges_list (list, None): List with nested lists or tuples containing start, end positions (Default value = None)
      task (str): Task performed. Must be hard_mask/soft_mask/exclude/extract (Default value = 'hard_mask')
      hard_mask_letter (str): Must be N/X (Default value = 'N')
      ranges_path (str, None): Filepath to two-column tab-separated file with start, end positions (Default value = None)
      inverse_mask (bool): Specifies whether to inverse ranges for hard/soft masking (Default value = False)

    Returns:
      None

    """
    # import packages and define functions
    from Bio import SeqIO
    from Bio.Seq import MutableSeq
    import os, sys, re, gzip, operator
    from functools import reduce
    from figleaf_fasta.utils import nested_list_elements_to_int, wrap_seq_string, mask_to_complement_indices, indices_to_slices

    fastapathregex = re.compile(r'.+\.(?:fa|fasta|fna)(?:\.gz)?$', re.IGNORECASE)

    # check masking arguments (task and hard_mask_letter)
    assert task in ['hard_mask', 'soft_mask', 'exclude', 'extract'], 'Error: task must be either "hard_mask", "soft_mask", "exclude", or "extract"'
    assert hard_mask_letter in ['N', 'X'], 'Error: hard_mask_letter must be either "N" or "X"'
    # check fasta_input
    assert os.path.isfile(fasta_input), 'Error: fasta_input must be a valid fasta file filepath'
    assert fastapathregex.match(fasta_input), 'Error: fasta_input must be a filepath to input fasta file, ending with .fasta/.fa/.fna(.gz)'
    # check fasta_output and make output directory
    output_path = os.path.dirname(fasta_output)
    assert fastapathregex.match(fasta_input), 'Error: fasta_output must be a filepath to output fasta file, ending with .fasta/.fa/.fna(.gz)'
    try:
        os.makedirs(output_path, exist_ok=True)
    except:
        sys.exit('Error: invalid fasta_output filepath')
    # get mask from ranges_path/ranges_list
    assert sum(i is not None for i in [ranges_path, ranges_list]) == 1, 'Error: mask must be provided either as ranges_path OR as ranges_list'
    if ranges_path is not None:
        # parse mask ranges file
        assert os.path.isfile(ranges_path), 'Error: ranges_path must be a valid filepath'
        mask = []
        with open(ranges_path) as f:
            for line in f:
                data = line.strip().split('\t')
                assert len(data) == 2, 'Error: ranges_path must be a two-column tsv file'
                start, stop = data
                mask.append((start, stop))
    else:
        mask = ranges_list
    # check mask is a list of len==2 tuples and get maximum position in mask
    assert isinstance(mask, list) and all(isinstance(item, (list, tuple)) and len(item) == 2 for item in mask), 'Error: mask must be a list of tuples/sublists of length 2'
    try:
        mask = nested_list_elements_to_int(mask)
    except:
        sys.exit('Error: elements of mask must be integers (or coerceable to integers)')
    assert all(item[1] > item[0] for item in mask), 'Error: end positions in mask ranges must be > than start positions (ranges must be 0-indexed, with non-inclusive end position)'
    reduced_mask = reduce(operator.concat, mask)
    mask_max = max(reduced_mask)
    # parse fasta and apply mask
    # read input
    if fasta_input.endswith('.gz'):
        input_handle = gzip.open(fasta_input, 'rt')
    else:
        input_handle = open(fasta_input, 'r')
    # parse input and write masked fasta to file
    if fasta_output.endswith('.gz'):
        output_handle = gzip.open(fasta_output, 'wt')
    else:
        output_handle = open(fasta_output, 'w')
    for recordid, recordseq in SeqIO.FastaIO.SimpleFastaParser(input_handle):
        recordseqlen = len(recordseq)
        assert recordseqlen >= mask_max, 'Error: mask range positions must not exceed the length of the fasta sequence; \
            max mask position is %s; sequence length for record %s is %s'%(mask_max, recordid, recordseqlen)
        # convert to MutableSeq
        Mutable_seq = MutableSeq(recordseq)

        # apply mask
        if task in ['hard_mask', 'soft_mask']:
            # inverse mask
            if inverse_mask is True:
                complement_indices = mask_to_complement_indices(mask, recordseqlen)
                mask = indices_to_slices(complement_indices)
            # apply hard/soft mask
            for mask_range in mask:
                start, end = mask_range
                if task == 'hard_mask':
                    Mutable_seq[start:end] = hard_mask_letter * (end - start)
                else:
                    Mutable_seq[start:end] = str(Mutable_seq[start:end]).lower()
            masked_seq = str(Mutable_seq)
        
        else:
            if task == 'exclude':
                # exclude mask ranges (select complement of mask ranges)
                complement_indices = mask_to_complement_indices(mask, recordseqlen)
                mask = indices_to_slices(complement_indices)
            extracted_seq = []
            for mask_range in mask:
                start, end = mask_range
                extracted_seq.append(str(Mutable_seq[start:end]))
            masked_seq = ''.join(extracted_seq)

        # wrap sequence to 60 characters per line, and write to file
        masked_seq = wrap_seq_string(masked_seq)
        output_handle.write(">%s\n%s\n" % (recordid, masked_seq))
    input_handle.close()
    output_handle.close()


def main():
    """Run this function if script is called from the command-line"""
    import argparse
    import os
    parser = argparse.ArgumentParser(description="figleaf-fasta: apply hard/soft mask to FASTA file or exclude/extract sub-sequences")
    # Input options
    input_group = parser.add_argument_group('Input')
    input_group.add_argument('-fi', '--fasta_input', help='Filepath to input fasta file to be masked (required)', required=True)
    input_group.add_argument('-r', '--ranges_path', help='Two-column tsv file with rows containing 0-indexed end-exclusive ranges to be masked/excluded/extracted (required)', required=True)
    # Output options
    output_group = parser.add_argument_group('Output')
    output_group.add_argument('-fo', '--fasta_output', help='Filepath for masked output fasta file (required)', required=True)
    # Task options
    task_group = parser.add_argument_group('Task')
    task_group.add_argument('--task', help='"hard_mask","soft_mask","exclude","extract" (default: hard_mask)', default='hard_mask')
    # Mask options
    mask_group = parser.add_argument_group('Mask')
    mask_group.add_argument('--hard_mask_letter', help='Letter to represent hard_mask regions (N or X) (default: N)', default='N')
    mask_group.add_argument('--inverse_mask', action='store_true', help='If flag is provided, all except mask ranges will be masked')

    args = parser.parse_args()

    # get relative file paths
    cwdir = os.getcwd()
    fasta_input_path = os.path.relpath(args.fasta_input, cwdir)
    ranges_input_path = os.path.relpath(args.ranges_path, cwdir)
    fasta_output_path = os.path.relpath(args.fasta_output, cwdir)

    # run figleaf
    figleaf(fasta_input=fasta_input_path, fasta_output=fasta_output_path,
            ranges_list=None, task=args.task, hard_mask_letter=args.hard_mask_letter,
            ranges_path=ranges_input_path, inverse_mask=args.inverse_mask)


if __name__ == "__main__":
    main()

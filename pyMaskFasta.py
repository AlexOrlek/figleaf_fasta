#! /usr/bin/env python3

#Usage: 
# import module in python script... from pyMaskFasta import maskfasta
# run script from the command line... run pyMaskFasta.py -h to see help on command-line options

def maskfasta(fasta_input,fasta_output,mask_list=None,mask_type='hard',hard_mask_letter='N',mask_path=None):
    """reads fasta_input filepath; applies mask provided either as a list of sublists/tuples containing start, end positions (mask_list) or as a filepath to a two-column tsv file containing rows of start \t end positions; writes masked fasta_output
        - mask positions must comprise 0-indexed start and end positions, with non-inclusive end position i.e. 1, 4 represents positions 2 - 4
        - fasta_input and fasta_output must have extensions .fasta/.fa/.fna(.gz),
        ... if there is a .gz extension, gzip will be used to read/write fasta_input/fasta_output
        - mask_type must be hard/soft/exclude; hard_mask_letter must be N/X"""
    ###import packages and define functions
    from Bio import SeqIO
    from Bio.Seq import MutableSeq
    import os, sys, re, gzip, operator
    from functools import reduce
    import numpy as np

    fastapathregex=re.compile(r'.+\.(?:fa|fasta|fna)(?:\.gz)?$',re.IGNORECASE)

    def nested_list_elements_to_int(x):
        """recursive function to convert elements of mask (list of nested tuples/sublists) to integers (run with try except)"""
        if isinstance(x, (list,tuple)):
            return list(map(nested_list_elements_to_int, x))
        else:
            return int(x)

    def wrap_seq_string(seq_string, n=60):
        """wrap nucleotide string so that each line of fasta has length n characters (default 60 characters)"""
        chunks = [seq_string[i:i+n] for i in range(0, len(seq_string), n)]
        return('\n'.join(chunks))

    ###check masking arguments (mask_type and hard_mask_letter)
    assert mask_type in ['hard','soft','exclude'], 'Error: mask_type must be either "hard", "soft" or "exclude"'
    assert hard_mask_letter in ['N','X'], 'Error: hard_mask_letter must be either "N" or "X"'
    ###check fasta_input
    assert os.path.isfile(fasta_input), 'Error: fasta_input must be a valid fasta file filepath'
    assert fastapathregex.match(fasta_input), 'Error: fasta_input must be a filepath to input fasta file, ending with .fasta/.fa/.fna(.gz)'
    ###check fasta_output and make output directory
    output_path=os.path.dirname(fasta_output)
    assert fastapathregex.match(fasta_input), 'Error: fasta_output must be a filepath to output fasta file, ending with .fasta/.fa/.fna(.gz)'
    try:
        os.makedirs(output_path,exist_ok=True)
    except:
        sys.exit('Error: invalid fasta_output filepath')
    ###get mask from mask_path/mask_list
    assert sum(i is not None for i in [mask_path,mask_list])==1, 'Error: mask must be provided either as mask_path OR as mask_list'
    if mask_path is not None:
        #parse mask file
        assert os.path.isfile(mask_path), 'Error: mask_path must be a valid filepath'
        mask=[]
        with open(mask_path) as f:
            for line in f:
                data=line.strip().split('\t')
                assert len(data)==2, 'Error: mask_path must be a two-column tsv file'
                start,stop=data
                mask.append((start,stop))
    else:
        mask=mask_list
    ###check mask is a list of len==2 tuples and get maximum position in mask
    assert isinstance(mask,list) and all(isinstance(item, (list,tuple)) and len(item)==2 for item in mask), 'Error: mask must be a list of tuples/sublists of length 2'
    try:
        mask=nested_list_elements_to_int(mask)
    except:
        sys.exit('Error: elements of mask must be integers (or coerceable to integers)')
    assert all(item[1]>item[0] for item in mask), 'Error: end positions in mask ranges must be > than start positions (ranges must be 0-indexed, with non-inclusive end position)'
    reduced_mask=reduce(operator.concat, mask)
    mask_max=max(reduced_mask)    
    ###parse fasta and apply mask
    #read input
    if fasta_input.endswith('.gz'):
        input_handle=gzip.open(fasta_input,'rt')
    else:
        input_handle=open(fasta_input,'r')
    #parse input and write masked fasta to file
    if fasta_output.endswith('.gz'):
        output_handle=gzip.open(fasta_output, 'wt')
    else:
        output_handle=open(fasta_output,'w')
    for recordid,recordseq in SeqIO.FastaIO.SimpleFastaParser(input_handle):
        assert len(recordseq)>=mask_max, 'Error: mask range positions must not exceed the length of the fasta sequence; \
            max mask position is %s; sequence length for record %s is %s'%(mask_max,recordid,len(recordseq))
        #convert to MutableSeq
        Mutable_seq = MutableSeq(recordseq)

        if mask_type in ['hard','soft']:
            #apply hard/soft mask
            for mask_range in mask:
                start,end=mask_range
                if mask_type=='hard':
                    Mutable_seq[start:end]=hard_mask_letter*(end-start)
                else:
                    Mutable_seq[start:end]=str(Mutable_seq[start:end]).lower()
        else:
            #exclude mask ranges
            mask_set=set()
            complete_set=set()
            complete_set.update(range(0,len(recordseq)))
            for mask_range in mask:
                start,end=mask_range
                mask_set.update(range(start,end))
            include_indices=list(complete_set.difference(mask_set))
            Mutable_seq=np.array(Mutable_seq)
            Mutable_seq=''.join(Mutable_seq[include_indices])

        #wrap sequence to 60 characters per line, and write to file
        masked_seq=wrap_seq_string(str(Mutable_seq))
        output_handle.write(">%s\n%s\n" % (recordid, masked_seq))
    output_handle.close()



if __name__ == "__main__":
    #this code is run if script is executed from command-line
    import argparse, os
    parser = argparse.ArgumentParser(description="pyMaskFasta: apply mask to FASTA file")
    #Input options                                                               
    input_group = parser.add_argument_group('Input')
    input_group.add_argument('-fi','--fasta_input', help='Filepath to input fasta file to be masked (required)', required=True)
    input_group.add_argument('-m','--mask_path', help='Two-column tsv file with rows containing 0-indexed end-exclusive ranges to be masked (required)', required=True)
    #Output options                                               
    output_group = parser.add_argument_group('Output')
    output_group.add_argument('-fo','--fasta_output', help='Filepath for masked output fasta file (required)', required=True)
    #Mask options
    mask_group = parser.add_argument_group('Mask')
    mask_group.add_argument('--mask_type', help='"hard" masking, "soft" masking, "exclude" (default: hard)', default='hard')
    mask_group.add_argument('--hard_mask_letter', help='Letter to represent hard masking (N or X) (default: N)', default='N')

    args = parser.parse_args()

    ###get relative file paths
    cwdir=os.getcwd()
    fasta_input_path=os.path.relpath(args.fasta_input, cwdir)
    mask_input_path=os.path.relpath(args.mask, cwdir)
    fasta_output_path=os.path.relpath(args.fasta_output, cwdir)

    ###run maskfasta
    maskfasta(fasta_input=fasta_input_path,fasta_output=fasta_output_path,mask_list=None,mask_type=args.mask_type,hard_mask_letter=args.hard_mask_letter,mask_path=mask_input_path)

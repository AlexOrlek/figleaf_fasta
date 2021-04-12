#!/usr/bin/env bash

fasta_input='example.fasta'
ranges_path='mask_ranges.tsv'

# hard-mask fasta (default)
fasta_output='output_sh/hard_masked.fasta'
figleaf -fi ${fasta_input} -r ${ranges_path} -fo ${fasta_output}

# output gzipped fasta
fasta_output='output_sh/hard_masked.fasta.gz' #output will be gzipped because of .gz extension
figleaf -fi ${fasta_input} -r ${ranges_path} -fo ${fasta_output}

# inverse mask
fasta_output='output_sh/hard_masked_inverse.fasta'
figleaf -fi ${fasta_input} -r ${ranges_path} -fo ${fasta_output} --inverse_mask

# soft-mask fasta
fasta_output='output_sh/soft_masked.fasta'
figleaf -fi ${fasta_input} -r ${ranges_path} -fo ${fasta_output} --task 'soft_mask'

# exclude mask ranges from fasta
fasta_output='output_sh/exclude_ranges.fasta'
figleaf -fi ${fasta_input} -r ${ranges_path} -fo ${fasta_output} --task 'exclude'

# extract mask ranges from fasta
fasta_output='output_sh/extract_ranges.fasta'
figleaf -fi ${fasta_input} -r ${ranges_path} -fo ${fasta_output} --task 'extract'

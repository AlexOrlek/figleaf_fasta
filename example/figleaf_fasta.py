from figleaf_fasta.figleaf import figleaf

mask = []
with open('mask_ranges.tsv') as f:
    for line in f:
        data = line.strip().split('\t')
        mask.append((data[0], data[1]))  # figleaf will convert strings to integers

fasta_input = 'example.fasta'

# hard-mask fasta (default)
fasta_output = 'output_py/hard_masked.fasta'
figleaf(fasta_input, fasta_output, mask)

# output gzipped fasta
fasta_output = 'output_py/hard_masked.fasta.gz'  # output will be gzipped because of .gz extension
figleaf(fasta_input, fasta_output, mask)

# mask ranges input from tsv file
fasta_output = 'output_py/hard_masked_from_tsv.fasta'
figleaf(fasta_input, fasta_output, ranges_path='mask_ranges.tsv')

# inverse mask
fasta_output = 'output_py/hard_masked_inverse.fasta'
figleaf(fasta_input, fasta_output, mask, inverse_mask=True)

# soft-mask fasta
fasta_output = 'output_py/soft_masked.fasta'
figleaf(fasta_input, fasta_output, mask, task='soft_mask')

# exclude mask ranges from fasta
fasta_output = 'output_py/exclude_ranges.fasta'
figleaf(fasta_input, fasta_output, mask, task='exclude')

# extract mask ranges from fasta
fasta_output = 'output_py/extract_ranges.fasta'
figleaf(fasta_input, fasta_output, mask, task='extract')

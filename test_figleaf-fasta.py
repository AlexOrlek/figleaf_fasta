from figleaf import maskfasta

mask=[]
with open('testing/mask_ranges.tsv') as f:
    for line in f:
        data=line.strip().split('\t')
        mask.append((data[0],data[1]))  #maskfasta will convert strings to integers

fasta_input='testing/test.fasta'

#hard-mask fasta (default)
fasta_output='testing/test.hard_masked.fasta'
maskfasta(fasta_input,fasta_output,mask)

#output gzipped fasta
fasta_output='testing/test.hard_masked.fasta.gz' #output will be gzipped because of .gz extension
maskfasta(fasta_input,fasta_output,mask)

#mask ranges input from tsv file
fasta_output='testing/test.hard_masked_from_tsv.fasta'
maskfasta(fasta_input,fasta_output,ranges_path='testing/mask_ranges.tsv')

#inverse mask
fasta_output='testing/test.hard_masked_inverse.fasta'
maskfasta(fasta_input,fasta_output,mask,inverse_mask=True)

#soft-mask fasta
fasta_output='testing/test.soft_masked.fasta'
maskfasta(fasta_input,fasta_output,mask,task='soft_mask')

#exclude mask ranges from fasta
fasta_output='testing/test.exclude_ranges.fasta'
maskfasta(fasta_input,fasta_output,mask,task='exclude')

#extract mask ranges from fasta
fasta_output='testing/test.extract_ranges.fasta'
maskfasta(fasta_input,fasta_output,mask,task='extract')

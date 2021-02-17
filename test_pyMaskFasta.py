from pyMaskFasta import maskfasta

mask=[]
with open('testing/mask_ranges.tsv') as f:
    for line in f:
        data=line.strip().split('\t')
        mask.append((data[0],data[1]))  #maskfasta will convert strings to integers

fasta_input='testing/test.fasta'

#hard-mask fasta
fasta_output='testing/test.hard_masked.fasta'
maskfasta(fasta_input,fasta_output,mask)

#output gzipped fasta
fasta_output='testing/test.hard_masked.fasta.gz' #output will be gzipped because of .gz extension
maskfasta(fasta_input,fasta_output,mask)

#mask input from tsv file
fasta_output='testing/test.hard_masked_from_tsv.fasta'
maskfasta(fasta_input,fasta_output,mask_path='testing/mask_ranges.tsv')

#soft-mask fasta
fasta_output='testing/test.soft_masked.fasta'
maskfasta(fasta_input,fasta_output,mask,mask_type='soft')

#exclude mask ranges from fasta
fasta_output='testing/test.exclude_masked.fasta'
maskfasta(fasta_input,fasta_output,mask,mask_type='exclude')


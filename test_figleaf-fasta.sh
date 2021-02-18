

fasta_input='testing/test.fasta'
ranges_path='testing/mask_ranges.tsv'

#hard-mask fasta (default)
fasta_output='testing/test.hard_masked.fasta'
python figleaf.py -fi ${fasta_input} -r ${ranges_path} -fo ${fasta_output}

#output gzipped fasta
fasta_output='testing/test.hard_masked.fasta.gz' #output will be gzipped because of .gz extension
python figleaf.py -fi ${fasta_input} -r ${ranges_path} -fo ${fasta_output}

#inverse mask
fasta_output='testing/test.hard_masked_inverse.fasta'
python figleaf.py -fi ${fasta_input} -r ${ranges_path} -fo ${fasta_output} --inverse_mask

#soft-mask fasta
fasta_output='testing/test.soft_masked.fasta'
python figleaf.py -fi ${fasta_input} -r ${ranges_path} -fo ${fasta_output} --task 'soft_mask'

#exclude mask ranges from fasta
fasta_output='testing/test.exclude_ranges.fasta'
python figleaf.py -fi ${fasta_input} -r ${ranges_path} -fo ${fasta_output} --task 'exclude'

#extract mask ranges from fasta
fasta_output='testing/test.extract_ranges.fasta'
python figleaf.py -fi ${fasta_input} -r ${ranges_path} -fo ${fasta_output} --task 'extract'

[![DOI](https://zenodo.org/badge/339805615.svg)](https://zenodo.org/badge/latestdoi/339805615)

figleaf_fasta applies hard/soft masking to a FASTA file or excludes/extracts sub-sequences from a FASTA file.<br>
* hard_mask: replace sequence with Ns or Xs
* soft_mask: convert sequence to lowercase
* exclude: exclude sub-sequences and concatenate non-excluded remainder
* extract: extract and concatenate sub-sequences
<br>

Other tools for handling FASTA files (e.g. `bedtools maskfasta`) often require a sequence name to be specified (in addition to range information), and the sequence name must match headers in the FASTA file.<br>

figleaf_fasta is a simple lightweight tool that takes as input a (multi-)FASTA and range start, end positions; masking/exclusion/extraction will be applied to sequence(s) within the (multi-)FASTA, regardless of FASTA header names.<br>

# Installation

## From pypi
```bash
pip3 install figleaf_fasta
```
## From GitHub repository
```bash
git clone https://github.com/AlexOrlek/figleaf_fasta.git
cd figleaf_fasta
pip3 install .
```

# Options and usage

figleaf_fasta can be run from a Linux command-line as follows:<br>
 `figleaf [`*`arguments...`*`]`

figleaf_fasta can be used within a Python script as follows:<br>
`from figleaf_fasta.figleaf import figleaf`<br>
`figleaf([`*`arguments...`*`])`<br>
<br>
Running `figleaf -h` on the command-line produces a summary of the command-line options:

```
usage: figleaf [-h] -fi FASTA_INPUT -r RANGES_PATH -fo FASTA_OUTPUT [--task TASK] [--hard_mask_letter HARD_MASK_LETTER] [--inverse_mask]

figleaf_fasta: apply hard/soft mask to FASTA file or exclude/extract sub-sequences

optional arguments:
  -h, --help            show this help message and exit

Input:
  -fi FASTA_INPUT, --fasta_input FASTA_INPUT
                        Filepath to input fasta file to be masked (required)
  -r RANGES_PATH, --ranges_path RANGES_PATH
                        Two-column tsv file with rows containing 0-indexed end-exclusive ranges to be masked/excluded/extracted (required)

Output:
  -fo FASTA_OUTPUT, --fasta_output FASTA_OUTPUT
                        Filepath for masked output fasta file (required)

Task:
  --task TASK           "hard_mask","soft_mask","exclude","extract" (default: hard_mask)

Mask:
  --hard_mask_letter HARD_MASK_LETTER
                        Letter to represent hard_mask regions (N or X) (default: N)
  --inverse_mask        If flag is provided, all except mask ranges will be masked
```

The same arguments are required when using the figleaf function within a Python script, except that start, end positions can be provided either as a filepath (`ranges_path`), OR as a Python list (`ranges_list`).


# Example

To generate example output in the example/ directory, run:<br>
`python figleaf_fasta.py` or `bash figleaf_fasta.sh`


# License

[MIT License](https://en.wikipedia.org/wiki/MIT_License)

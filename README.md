[![DOI](https://zenodo.org/badge/339718171.svg)](https://zenodo.org/badge/latestdoi/339718171)

figleaf-fasta applies hard/soft masking to a FASTA file or excludes/extracts sub-sequences from a FASTA file.<br>
* hard_mask: replace sequence with Ns or Xs
* soft_mask: convert sequence to lowercase
* exclude: exclude sub-sequences and concatenate non-excluded remainder
* extract: extract and concatenate sub-sequences
<br>

Other tools for handling FASTA files (e.g. `bedtools maskfasta`) often require a sequence name to be specified (in addition to range information), and the sequence name must match headers in the FASTA file.<br>

figleaf-fasta is a simple lightweight tool that takes as input a (multi-)FASTA and range start, end positions; masking/exclusion/extraction will be applied to sequence(s) within the (multi-)FASTA, regardless of FASTA header names.<br>

# Requirements

* [Python 3](https://www.python.org/)
* [Biopython](https://biopython.org/)


# Installation

```bash
git clone https://github.com/AlexOrlek/figleaf-fasta.git
cd ATCG
```

# Options and usage

First add the path of the figleaf-fasta directory to the [$PATH](https://www.computerhope.com/issues/ch001647.htm) and [$PYTHONPATH](https://bic-berkeley.github.io/psych-214-fall-2016/using_pythonpath.html) variables so that the software can be run from any location.<br>

figleaf-fasta can be run from a Linux command-line as follows:<br>
 `figleaf.py [`*`arguments...`*`]`

figleaf-fasta can be used within a Python script as follows:<br>
`from figleaf-fasta import maskfasta`<br>
`maskfasta([`*`arguments...`*`])`<br>
<br>
Running `figleaf.py -h` on the command-line produces a summary of the command-line options:

```
usage: figleaf.py [-h] -fi FASTA_INPUT -r RANGES_PATH -fo FASTA_OUTPUT [--task TASK] [--hard_mask_letter HARD_MASK_LETTER] [--inverse_mask]

figleaf-fasta: apply hard/soft mask to FASTA file or exclude/extract sub-sequences

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

The same arguments are required when using the figleaf-fasta.maskfasta() function within a Python script, except that start, end positions can be provided either as a filepath (`ranges_path`), OR as a Python list (`ranges_list`).


# Example

To generate example output in the testing/ directory, run:<br>
`python test_figleaf-fasta.py` or `bash test_figleaf-fasta.sh`


# License

[MIT License](https://en.wikipedia.org/wiki/MIT_License)

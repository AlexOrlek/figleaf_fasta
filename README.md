[![DOI](https://zenodo.org/badge/339718171.svg)](https://zenodo.org/badge/latestdoi/339718171)

pyMaskFasta masks specified positions in a FASTA file. The mask can be a hard-mask (replace with Ns), soft-mask (convert to lowercase), or the positions can be excluded.<br>

Other tools for masking FASTAs (e.g. `bedtools maskfasta`) generally require a sequence name to be specified as well as mask positions, and the sequence name must match headers in the FASTA file.<br>

pyMaskFasta is a lightweight simple tool that takes as input a (multi-)FASTA and mask positions; unlike bedtools maskfasta, the mask will be applied to sequence(s) within the (multi-)FASTA, regardless of FASTA header names.<br>

# Requirements

* [Python 3](https://www.python.org/)
* [Biopython](https://biopython.org/)


# Installation

```bash
git clone https://github.com/AlexOrlek/pyMaskFasta.git
cd ATCG
```

# Options and usage

First add the path of the pyMaskFasta directory to the [$PATH](https://www.computerhope.com/issues/ch001647.htm) and [$PYTHONPATH](https://bic-berkeley.github.io/psych-214-fall-2016/using_pythonpath.html) variables so that the software can be run from any location.<br>

pyMaskFasta can be run from a Linux command-line as follows:<br>
 `pyMaskFasta.py [`*`arguments...`*`]`

pyMaskFasta can be used within a Python script as follows:<br>
`from pyMaskFasta import maskfasta`<br>
`maskfasta([`*`arguments...`*`])`<br>
<br>
Running `pyMaskFasta.py -h` on the command-line produces a summary of the command-line options:

```
usage: pyMaskFasta.py [-h] -fi FASTA_INPUT -m MASK_PATH -fo FASTA_OUTPUT [--mask_type MASK_TYPE] [--hard_mask_letter HARD_MASK_LETTER]

pyMaskFasta: apply mask to FASTA file

optional arguments:
  -h, --help            show this help message and exit

Input:
  -fi FASTA_INPUT, --fasta_input FASTA_INPUT
                        Filepath to input fasta file to be masked (required)
  -m MASK_PATH, --mask_path MASK_PATH
                        Two-column tsv file with rows containing 0-indexed end-exclusive ranges to be masked (required)

Output:
  -fo FASTA_OUTPUT, --fasta_output FASTA_OUTPUT
                        Filepath for masked output fasta file (required)

Mask:
  --mask_type MASK_TYPE
                        "hard" masking, "soft" masking, "exclude" (default: hard)
  --hard_mask_letter HARD_MASK_LETTER
                        Letter to represent hard masking (N or X) (default: N)
```

The same arguments are required when using the pyMaskFasta maskfasta() function within a Python script, except that mask information can be provided either as a filepath (`mask_path`), OR as a Python list (`mask_list`).


# Example

To generate example output in the testing/ directory, run:<br>
`python test_pyMaskFasta.py`


# License

[MIT License](https://en.wikipedia.org/wiki/MIT_License)

from itertools import groupby


def nested_list_elements_to_int(mask):
    """Recursively convert list of nested lists/tuples to list of nested lists with elements coerced to integers

    Args:
        mask (list): List with nested lists/tuples containing start, end positions

    Returns:
        (list): List of nested lists containing start, end positions, coerced to integers

    """
    if isinstance(mask, (list, tuple)):
        return list(map(nested_list_elements_to_int, mask))
    else:
        return int(mask)


def wrap_seq_string(seq_string, n=60):
    """Wrap nucleotide string (seq_string) so that each line of fasta has length n characters (default n=60)
 
    Args:
        seq_string (str): String of DNA/protein sequence
        n (int): Maximum number of characters to display per line

    Returns:
        (str): String of DNA/protein sequence wrapped to n characters

    """
    chunks = [seq_string[i:i + n] for i in range(0, len(seq_string), n)]
    return('\n'.join(chunks))


def mask_to_complement_indices(mask, recordseqlen):
    """Given a list of nested lists containing start, end positions (mask), return sequence positions not contained in the mask as a sorted integer list

    Args:
        mask (list): List with nested lists containing start, end integer positions
        recordseqlen (int): Length of FASTA sequence

    Returns:
        complement_indices (list): Sorted list of integer positions (complement of mask)

    """
    mask_set = set()
    complete_set = set()
    complete_set.update(range(0, recordseqlen))
    for mask_range in mask:
        start, end = mask_range
        mask_set.update(range(start, end))
    complement_indices = sorted(list(complete_set.difference(mask_set)))
    return(complement_indices)


def indices_to_slices(indices):
    """Convert a list of indices to a list of nested lists containing slices, using itertools.groupby

    Args:
        indices (list): Sorted list of integer positions

    Returns:
        slices (list): List of nested lists with start, end positions

    """
    slices = []
    for key, it in groupby(enumerate(indices), lambda x: x[1] - x[0]):
        indices = [y for x, y in it]
        if len(indices) == 1:
            slices.append([indices[0], indices[0] + 1])
        else:
            slices.append([indices[0], indices[-1] + 1])
    return slices

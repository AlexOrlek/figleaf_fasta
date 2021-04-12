#!/usr/bin/env python

"""Tests for `figleaf_fasta` package."""

from figleaf_fasta.utils import nested_list_elements_to_int, wrap_seq_string, mask_to_complement_indices, indices_to_slices


def test_nested_list_elements_to_int():
    # Check that nested list elements are converted to integers
    assert nested_list_elements_to_int([('1', '2'), ('3', '4')]) == [[1, 2], [3, 4]]
    assert nested_list_elements_to_int([['1', '2'], ['3', '4']]) == [[1, 2], [3, 4]]


def test_wrap_seq_string():
    # Check that nucleotide sequence is wrapped to 60 characters
    seq_120 = ''.join(['ATCG' for i in range(30)])
    seq_120_wrapped = '\n'.join([seq_120[0:60], seq_120[60:120]])
    assert wrap_seq_string(seq_120) == seq_120_wrapped


def test_mask_to_complement_indices():
    # Check that mask is converted to complement integer positions
    assert mask_to_complement_indices(mask=[[1, 2], [8, 10]], recordseqlen=10) == [0, 2, 3, 4, 5, 6, 7]


def test_indices_to_slices():
    # Check that ordered integer positions are converted to slices
    assert indices_to_slices([0, 2, 3, 4, 5, 6, 7]) == [[0, 1], [2, 8]]

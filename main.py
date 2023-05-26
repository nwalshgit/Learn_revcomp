import itertools
import pathlib
import sys
import timeit

import revcomp

d = revcomp.nucleotide_ambiguity_table


def expand_1(seq: str) -> list[str]:
    """Expand a degenerate sequence into set of possible sequences"""
    # This is not memory efficient
    expanded_seqs = [""]
    for degen_base in seq:
        if degen_base.upper() in d:
            replacements = d[degen_base.upper()]
        else:
            replacements = set(degen_base)
        prior_seqs = expanded_seqs.copy()
        expanded_seqs = []
        for base in replacements:
            for seq in prior_seqs:
                expanded_seqs.append(seq+base)
        # print(expanded_seqs, "\n----\n")
    return expanded_seqs


def expand_2(seq):
    [d[j] for j in seq]
    return ["".join(i) for i in itertools.product(*[d[j] for j in seq])]


def expand_3(seq):
    return list(map("".join, itertools.product(*map(d.get, seq))))


def expand_4(seq):
    groups = itertools.groupby(seq, lambda char: char not in d)
    splits = []
    for b, group in groups:
        if b:
            splits.extend([[g] for g in group])
        else:
            for nuc in group:
                splits.append(d[nuc])
    return [''.join(p) for p in itertools.product(*splits)]


def compare_expand():
    degen = "NATHANRRRRR"
    for expand_func in [expand_3, expand_2, expand_4, expand_1]:
        timeit_statement = f"{expand_func.__name__}('{degen}')"
        timeit_setup = f"from __main__ import {expand_func.__name__}"
        print(expand_func.__name__,
              f"{timeit.timeit(timeit_statement, number=10000, setup=timeit_setup):0.4f}",
              degen, "->",
              sorted(expand_func(degen)))


if __name__ == "__main__":
    # if we pass in something on the command line
    if len(sys.argv) > 1:
        # if what we pass is a file
        if pathlib.Path(sys.argv[1]).exists():
            # open the file
            with open(pathlib.Path(sys.argv[1]), 'r') as seq_file:
                seq1 = seq_file.read
        # if what we pas is not a file, treat it as a sequence
        else:
            seq1 = sys.argv[1]
    # otherwise use a pre-determined sequence
    else:
        seq1 = "TAGTCGCCTGAAGCC"
    rc1 = revcomp.rev_comp(seq1)
    print(rc1)

    compare_expand()

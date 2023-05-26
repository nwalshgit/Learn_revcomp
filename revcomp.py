"""This is the revccomp module
It implements rev_comp to convert a nucleotide
sequence into its reverse complement"""


nucleotide_ambiguity_table: dict = {
    'A': {'A'},
    'G': {'G'},
    'C': {'C'},
    'T': {'T'},
    'U': {'U'},
    'N': {'A', 'C', 'G', 'T'},  # aNy
    'R': {'A', 'G'},            # puRine
    'Y': {'C', 'T'},            # pYrimidine
    'K': {'G', 'T'},            # Keto
    'M': {'A', 'C'},            # aMino
    'S': {'C', 'G'},            # Strong
    'W': {'A', 'T'},            # Weak
    'B': {'C', 'G', 'T'},       # Not A
    'D': {'A', 'G', 'T'},       # Not C
    'H': {'A', 'C', 'T'},       # Not G
    'V': {'A', 'C', 'G'},       # Not T/U
}
# Used by rev_comp_3 and rev_comp_5, but only needs to be initialized once
dna_complement: dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
rna_complement: dict = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
# Used by rev_comp_5, but only needs to be initialized once
alt_map = {'ins': '0'}
# Used by rev_comp_6, but only needs to be initialized once
# maketrans creates a table similar to above but uses the ordinal integer of the letter
dna_table: dict = str.maketrans("ACGTMRWSYKVHDBN acgtmrwsykvhdbn",
                                "TGCAKYWSRMBDHVN tgcakywsrmbdhvn")
rna_table: dict = str.maketrans("ACGUMRWSYKVHDBN acgumrwsykvhdbn",
                                "UGCAKYWSRMBDHVN ugcakywsrmbdhvn")


# This is only used rev_comp_1 and rev_comp_2 (neither is the default)
def verify(sequence: str) -> str:
    """This code verifies if a sequence is a DNA or RNA
    caveats:
    -It requires the sequence to have ALL four bases
    -It can't handle degeneracy
    -It can't handle lowercase
    -It can't handle other characters (spaces, tabs, newlines, etc...)

    >>> verify('ATGC')
    'DNA'
    >>> verify('AUGC')
    'RNA'
    >>> verify('Non-DNA text.')
    'Invalid sequence'
    >>> verify('GGGG')
    'Invalid sequence'
    >>> verify('atgc')
    'Invalid sequence'
    >>> verify('ADhN')
    'Invalid sequence'
    """
    # set the input sequence
    seq: set = set(sequence)

    # confirm if its elements is equal to
    # the set of valid DNA bases
    # Use a union method to ensure the
    # sequence is verified if does not
    # contain all the bases
    if seq == {"A", "T", "C", "G"}.union(seq):
        return "DNA"
    elif seq == {"A", "U", "C", "G"}.union(seq):
        return "RNA"
    else:
        return "Invalid sequence"


def rev_comp_1(seq: str) -> str:
    """Reverse Complement using replace with lowercase intermediates
    Caveats:
    -It requires the sequence to have ALL four bases
    -It can't handle degeneracy
    -It can't handle lowercase
    -It can't handle other characters (spaces, tabs, newlines, etc...)

    >>> rev_comp_1('ATGC')
    'GCAT'
    >>> rev_comp_1('AUGC')
    'GCAU'
    >>> rev_comp_1('Non-DNA text.')
    'Invalid sequence'
    >>> rev_comp_1('GGGG')
    'Invalid sequence'
    >>> rev_comp_1('atgc')
    'Invalid sequence'
    >>> rev_comp_1('ADhN')
    'Invalid sequence'
    """
    verified: str = verify(seq)
    if verified == "DNA":

        # complement strand
        seq: str = seq.replace("A", "t").replace(
            "C", "g").replace("T", "a").replace("G", "c")
        seq = seq.upper()

        # reverse strand
        seq = seq[::-1]
        return seq

    elif verified == "RNA":

        # complement strand
        seq: str = seq.replace("A", "u").replace(
            "C", "g").replace("U", "a").replace("G", "c")
        seq = seq.upper()

        # reverse strand
        seq = seq[::-1]
        return seq
    else:
        return "Invalid sequence"


def rev_comp_2(seq: str) -> str:
    """
    Reverse Complement using explicit if-else statements, using new list as intermediate
    Caveats:
    -It requires the sequence to have ALL four bases
    -It can't handle degeneracy
    -It can't handle lowercase, because it uses lowercase as intermediate
    -It can't handle other characters (spaces, tabs, newlines, etc...)

    >>> rev_comp_2('ATGC')
    'GCAT'
    >>> rev_comp_2('AUGC')
    'GCAU'
    >>> rev_comp_2('Non-DNA text.')
    'Invalid Sequence'
    >>> rev_comp_2('GGGG')
    'Invalid Sequence'
    >>> rev_comp_2('atgc')
    'Invalid Sequence'
    >>> rev_comp_2('ADhN')
    'Invalid Sequence'
    """
    comp: list = []
    if verify(seq) == "DNA":
        for base in seq:
            if base == "A":
                comp.append("T")
            elif base == "G":
                comp.append("C")
            elif base == "T":
                comp.append("A")
            elif base == "C":
                comp.append("G")
    elif verify(seq) == "RNA":
        for base in seq:
            if base == "U":
                comp.append("A")
            elif base == "G":
                comp.append("C")
            elif base == "A":
                comp.append("U")
            elif base == "C":
                comp.append("G")
    else:
        return "Invalid Sequence"

    # reverse the sequence
    comp_rev: list = comp[::-1]  # Starts as a list...

    # convert list to string
    comp_rev: str = "".join(comp_rev)  # BAD CHOICE to reuse same variable as a str...
    return comp_rev


def rev_comp_3(seq: str) -> str:
    """
    Reverse Complement using loop over bases with dictionary lookup
    Caveats:
    -It doesn't recognize degeneracy
    -It doesn't recognize lowercase
    -It assumes DNA unless there is a "U"
    -It ignores other characters (spaces, tabs, newlines, etc...)

    >>> rev_comp_3('ATGC')
    'GCAT'
    >>> rev_comp_3('AUGC')
    'GCAU'
    >>> rev_comp_3('Non-DNA text.')
    '.txet TND-noN'
    >>> rev_comp_3('GGGG')
    'CCCC'
    >>> rev_comp_3('atgc')
    'cgta'
    >>> rev_comp_3('ADhN')
    'NhDT'
    """
    if 'U' in seq:
        reverse_complement: str = "".join(rna_complement.get(base, base) for base in reversed(seq))
    else:
        reverse_complement: str = "".join(dna_complement.get(base, base) for base in reversed(seq))
    return reverse_complement


def rev_comp_4(seq: str) -> str:
    """
     Reverse Complement using translate with dictionary of ordinal values
     Caveats:
     -It doesn't recognize degeneracy
     -It doesn't recognize lowercase
     -It assumes DNA unless there is a "U"
     -It ignores other characters (spaces, tabs, newlines, etc...)

     >>> rev_comp_4('ATGC')
     'GCAT'
     >>> rev_comp_4('AUGC')
     'GCAU'
     >>> rev_comp_4('Non-DNA text.')
     '.txet TND-noN'
     >>> rev_comp_4('GGGG')
     'CCCC'
     >>> rev_comp_4('atgc')
     'cgta'
     >>> rev_comp_4('ADhN')
     'NhDT'
     """
    if 'U' in seq:
        old_chars: str = "ACGU"
        replace_chars: str = "UGCA"
    else:
        old_chars: str = "ACGT"
        replace_chars: str = "TGCA"
    tab: dict = str.maketrans(old_chars, replace_chars)
    return seq.translate(tab)[::-1]


def rev_comp_5(seq: str) -> str:
    """
     Reverse Complement using loop over bases with dictionary lookup
     It can protect substrings, but using number as an intermediate, so can't handle numbers
     Caveats:
     -It doesn't recognize degeneracy
     -It doesn't recognize lowercase
     -It assumes DNA unless there is a "U"
     -It ignores other characters (spaces, tabs, newlines, etc...)

     >>> rev_comp_5('ATGC')
     'GCAT'
     >>> rev_comp_5('AUGC')
     'GCAU'
     >>> rev_comp_5('Non-DNA text.')
     '.txet TND-noN'
     >>> rev_comp_5('GGGG')
     'CCCC'
     >>> rev_comp_5('atgc')
     'cgta'
     >>> rev_comp_5('ADhN')
     'NhDT'
     """
    k: str
    v: str
    for k, v in alt_map.items():
        seq: str = seq.replace(k, v)
    bases = list(seq)
    if 'U' in seq:
        bases = reversed([rna_complement.get(base, base) for base in bases])
    else:
        bases = reversed([dna_complement.get(base, base) for base in bases])
    bases: str = ''.join(bases)  # BAD CHOICE to reuse same variable as a str...
    for k, v in alt_map.items():
        bases = bases.replace(v, k)
    return bases


def rev_comp_6(seq: str, seq_type: str = 'DNA') -> str:
    """Perform a fast reverse complement of a nucleic acid strand
    seq: Nucleic acid sequence to be used as the sense strand
    seq_type: if "RNA" then expect and return "U" and not "T"
    return: reverse complement of the starting string

    caveats:
    - It assumes DNA unless there is a "U"
    - It ignores other characters (spaces, tabs, newlines, etc...)
    - Cannot transform DNA into RNA or vice versa

    >>> rev_comp_6('ATGC')
    'GCAT'
    >>> rev_comp_6('AUGC')
    'GCAU'
    >>> rev_comp_6('Non-DNA text.')
    '.axea TNH-noN'
    >>> rev_comp_6('GGGG')
    'CCCC'
    >>> rev_comp_6('atgc')
    'gcat'
    >>> rev_comp_6('ADhN')
    'NdHT'
    """
    if seq_type == "RNA" or "U" in seq:
        return seq.translate(rna_table)[::-1]
    return seq.translate(dna_table)[::-1]


# Because rev_comp_6 is the "best" we want that to be the default
rev_comp = rev_comp_6


def run_speed_test():
    import timeit
    # test variables
    seq1: str = "ATGCAGCTGTGTTACGCGAT"
    seq2: str = "UGGCGGAUAAGCGCA"
    seq3: str = "TYHGGHHHHH"

    output: str = "The reverse complementary strand of"
    for seq in [seq1, seq2, seq3]:
        output += f",\t{seq}"
    print(output)
    for rev_comp_func in [
            rev_comp, rev_comp_4, rev_comp_1, rev_comp_3, rev_comp_5, rev_comp_2]:
        output = f"                          {rev_comp_func.__name__}:"
        for seq in [seq1, seq2, seq3]:
            output += f",\t{rev_comp_func(seq)}"
        timeit_statement: str = f"{rev_comp_func.__name__}('{seq1}')"
        timeit_setup: str = f"from __main__ import {rev_comp_func.__name__}"
        output += f",\t{timeit.timeit(timeit_statement, number=100000, setup=timeit_setup)}"
        print(output)


if __name__ == "__main__":
    run_speed_test()


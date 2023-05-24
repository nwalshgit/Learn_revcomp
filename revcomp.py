"""This is the revccomp module
It implements rev_comp to convert a nucleotide
sequence into its reverse complement"""


# This is only used rev_comp_1 and rev_comp_2 (neither is the default)
def verify(sequence):
    """This code verifies if a sequence is a DNA or RNA"""
    # set the input sequence
    seq = set(sequence)

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


def rev_comp_1(seq):
    """This function returns a reverse complement
    of a DNA or RNA strand"""
    verified = verify(seq)
    if verified == "DNA":

        # complement strand
        seq = seq.replace("A", "t").replace(
            "C", "g").replace("T", "a").replace("G", "c")
        seq = seq.upper()

        # reverse strand
        seq = seq[::-1]
        return seq

    elif verified == "RNA":

        # complement strand
        seq = seq.replace("A", "u").replace(
            "C", "g").replace("U", "a").replace("G", "c")
        seq = seq.upper()

        # reverse strand
        seq = seq[::-1]
        return seq
    else:
        return "Invalid sequence"


def rev_comp_2(seq):
    comp = []
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
    comp_rev = comp[::-1]

    # convert list to string
    comp_rev = "".join(comp_rev)
    return comp_rev


# Used by rev_comp_3, but only need to be initialized once
dna_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
rna_complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}


def rev_comp_3(seq):
    if 'U' in seq:
        reverse_complement = "".join(rna_complement.get(base, base) for base in reversed(seq))
    else:
        reverse_complement = "".join(dna_complement.get(base, base) for base in reversed(seq))
    return reverse_complement


def rev_comp_4(seq):
    if 'U' in seq:
        old_chars = "ACGU"
        replace_chars = "UGCA"
    else:
        old_chars = "ACGT"
        replace_chars = "TGCA"
    tab = str.maketrans(old_chars, replace_chars)
    return seq.translate(tab)[::-1]


# Used by rev_comp_5, but only need to be initialized once
alt_map = {'ins': '0'}


def rev_comp_5(seq):
    # DNA only...
    for k, v in alt_map.items():
        seq = seq.replace(k, v)
    bases = list(seq) 
    if 'U' in seq:
        bases = reversed([rna_complement.get(base, base) for base in bases])
    else:
        bases = reversed([dna_complement.get(base, base) for base in bases])
    bases = ''.join(bases)
    for k, v in alt_map.items():
        bases = bases.replace(v, k)
    return bases


nucleotide_ambiguity_table = {
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
# Used by rev_comp_6, but only need to be initialized once
dna_table = str.maketrans("ACGTMRWSYKVHDBN acgtmrwsykvhdbn",
                          "TGCAKYWSRMBDHVN tgcakywsrmbdhvn")
rna_table = str.maketrans("ACGUMRWSYKVHDBN acgumrwsykvhdbn",
                          "UGCAKYWSRMBDHVN ugcakywsrmbdhvn")


def rev_comp_6(seq: str, seq_type: str = 'DNA'):
    """Perform a fast reverse complement of a nucleic acid strand
    seq: Nucleic acid sequence to be used as the sense strand
    seq_type: if "RNA" the expect and return "U" and not "T"
    return: reverse complement of the starting string

    caveats:
    - Unrecognized characters are not modified
    - Cannot transform DNA into RNA or vice versa
    """
    if seq_type == "RNA" or "U" in seq:
        return seq.translate(rna_table)[::-1]
    return seq.translate(dna_table)[::-1]


# Because rev_comp_6 is the "best" we want that to be the default
rev_comp = rev_comp_6


def run_speed_test():
    import timeit
    # test variables
    seq1 = "ATGCAGCTGTGTTACGCGAT"
    seq2 = "UGGCGGAUAAGCGCA"
    seq3 = "TYHGGHHHHH"

    output = "The reverse complementary strand of"
    for seq in [seq1, seq2, seq3]:
        output += f",\t{seq}"
    print(output)
    for rev_comp_func in [
            rev_comp, rev_comp_4, rev_comp_1, rev_comp_3, rev_comp_5, rev_comp_2]:
        output = f"                          {rev_comp_func.__name__}:"
        for seq in [seq1, seq2, seq3]:
            output += f",\t{rev_comp_func(seq)}"
        timeit_statement = f"{rev_comp_func.__name__}('{seq1}')"
        timeit_setup = f"from __main__ import {rev_comp_func.__name__}"
        output += f",\t{timeit.timeit(timeit_statement, number=100000, setup=timeit_setup)}"
        print(output)


if __name__ == "__main__":
    run_speed_test()

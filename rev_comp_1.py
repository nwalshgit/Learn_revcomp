def verify(sequence):
    '''This code verifies if a sequence is a DNA or RNA'''
     
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
    '''This function returns a reverse complement
    of a DNA or RNA strand'''
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
    tab = str.maketrans(old_chars,replace_chars)
    return seq.translate(tab)[::-1]


alt_map = {'ins':'0'}
def rev_comp_5(seq):
    # DNA only...
    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq) 
    if 'U' in seq:
        bases = reversed([rna_complement.get(base,base) for base in bases])
    else:
        bases = reversed([dna_complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases
    
# test variables
seq1 = "ATGCAGCTGTGTTACGCGAT"
seq2 = "UGGCGGAUAAGCGCA"
seq3 = "TYHGGHHHHH"

output = "The reverse complementary strand of"
for seq in [seq1, seq2, seq3]:
    output += f",\t{seq}"
print(output)
for index, rev_comp in enumerate([rev_comp_1, rev_comp_2, rev_comp_3, rev_comp_4, rev_comp_5]):
    output = f"                          method {index+1}:"
    for seq in [seq1, seq2, seq3]:
        output += f",\t{rev_comp(seq)}"
    print(output)

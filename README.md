# Final_450C
This is my actual final project for BISC 450C python 2

# Sequencing Objects (1-4)
In this analysis, we looked at sequencing objects.

```python
from Bio.Seq import Seq
```


```python
my_seq = Seq("GATCG")
```


```python
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
from Bio.Seq import Seq
```


```python
my_seq = Seq("GATCG")
```


```python
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
# we can also print the length of each sequence
print(len(my_seq))
```

    5



```python
print(my_seq[0])
```

    G



```python
print(my_seq[4])
```

    G



```python
print(my_seq[2])
```

    T



```python
# count the occurences of ""
Seq("AAAA").count("AA")
```




    2




```python
my_seq = Seq("GATCGTAAACCGGTTAAGCTCTCGATAGCTAGCT")
```


```python
len(my_seq)
```




    34




```python
# we can count specific letters in a sequence
my_seq.count("G")
```




    8




```python
100* ((my_seq.count("G") + my_seq.count("C")) / len(my_seq))
```




    47.05882352941176




```python
#imprted the GC function
from Bio.SeqUtils import gc_fraction
```


```python
my_seq = Seq("GATCGTAAACCGGTTAAGCTCTCGATAGCTAGCT")
```


```python
gc_fraction(my_seq)
```




    0.47058823529411764




```python
# we can cut out parts of the sequence, stride
my_seq[0 : : 3]
```




    Seq('GCACGACTAGAT')




```python
my_seq[1 : : 3]
```




    Seq('AGACTATCTCG')




```python
my_seq[2:3]
```




    Seq('T')




```python
# print the whole thing backwards with negative
my_seq[: :-1]
```




    Seq('TCGATCGATAGCTCTCGAATTGGCCAAATGCTAG')




```python
#turn seq back into string
str(my_seq)
```




    'GATCGTAAACCGGTTAAGCTCTCGATAGCTAGCT'




```python
fasta_format_string = ">Name\n%s\n" % my_seq
```


```python
print(fasta_format_string)
```

    >Name
    GATCGTAAACCGGTTAAGCTCTCGATAGCTAGCT
    



```python
#can add together sequences
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
```


```python
seq1 + seq2
```




    Seq('ACGTAACCGG')




```python
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
```


```python
spacer = Seq("N" *10)
```


```python
# take spacer object we made and join with the contigs
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTGCA')




```python
dna_seq = Seq("acgtACGT")
```


```python
dna_seq
```




    Seq('acgtACGT')




```python
#can convert to upper or lowercase
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
dna_seq.lower()
```




    Seq('acgtacgt')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
"gtac" in dna_seq
```




    False




```python
"GTAC" in dna_seq
```




    False




```python
dna_seq = dna_seq.upper()
```


```python
"GTAC" in dna_seq
```




    True




```python
my_seq = Seq("GTTTCACACACGATATACATCAATTGAGATTACATCATCATCTGTAGGAT")
```


```python
#find the complement
my_seq.complement()
```




    Seq('CAAAGTGTGTGCTATATGTAGTTAACTCTAATGTAGTAGTAGACATCCTA')




```python
#find the reverse complement
my_seq.reverse_complement()
```




    Seq('ATCCTACAGATGATGATGTAATCTCAATTGATGTATATCGTGTGTGAAAC')




```python
protein_seq = Seq("EVRNAK")
protein_seq.complement()
```




    Seq('EBYNTM')




```python
coding_dna = Seq("ACTGGATATATCACCATAGATCATAGTCGGTCTATGTCAGATACTACTA")
```


```python
coding_dna
```




    Seq('ACTGGATATATCACCATAGATCATAGTCGGTCTATGTCAGATACTACTA')




```python
template_dna = coding_dna.reverse_complement()
```


```python
template_dna 
```




    Seq('TAGTAGTATCTGACATAGACCGACTATGATCTATGGTGATATATCCAGT')




```python
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna
```




    Seq('ACUGGAUAUAUCACCAUAGAUCAUAGUCGGUCUAUGUCAGAUACUACUA')




```python
#can do both at once
template_dna.reverse_complement().transcribe()
```




    Seq('ACUGGAUAUAUCACCAUAGAUCAUAGUCGGUCUAUGUCAGAUACUACUA')




```python
# reverse transcription
messenger_rna.back_transcribe()
```




    Seq('ACTGGATATATCACCATAGATCATAGTCGGTCTATGTCAGATACTACTA')




```python
messenger_rna
```




    Seq('ACUGGAUAUAUCACCAUAGAUCAUAGUCGGUCUAUGUCAGAUACUACUA')




```python
messenger_rna.translate()
```

    /home/student/anaconda3/lib/python3.7/site-packages/Bio/Seq.py:2808: BiopythonWarning: Partial codon, len(sequence) not a multiple of three. Explicitly trim the sequence or add trailing N before translation. This may become an error in future.
      BiopythonWarning,





    Seq('TGYITIDHSRSMSDTT')




```python
coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('TGYITMDHSRSMSDTT')




```python
coding_dna.translate(table = 2)
```




    Seq('TGYITMDHSRSMSDTT')




```python
coding_dna.translate(to_stop = True)
```




    Seq('TGYITIDHSRSMSDTT')




```python
#no premature stop codon in mitochondrial DNA
coding_dna.translate(table = 2, to_stop=True)
```




    Seq('TGYITMDHSRSMSDTT')




```python
coding_dna.translate(table = 2, stop_symbol = "!")
```




    Seq('TGYITMDHSRSMSDTT')




```python
gene = Seq("GTGAAAAAACATCAATACTGGATATATCACCATAAATCATATTCGGGGGCCCTTTACGGTAATCACGTTACGGCATGTACAGATCATCTATGTCAGATACTCUAG")
```


```python
gene.translate(table = "Bacterial")
```




    Seq('VKKHQYWIYHHKSYSGALYGNHVTACTDHLCQIL*')




```python
gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKHQYWIYHHKSYSGALYGNHVTACTDHLCQIL')




```python
gene.translate(table = "Bacterial", cds = True)
```




    Seq('MKKHQYWIYHHKSYSGALYGNHVTACTDHLCQIL')




```python
from Bio.Data import CodonTable
```


```python
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
#visualize the table
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
#sequence comparison
seq = Seq("ACGT")
```


```python
"ACGT" == seq1
```




    True




```python
seq1 =="ACGT"
```




    True




```python
unknown_seq = Seq(None, 10)
```


```python
unknown_seq
```




    Seq(None, length=10)




```python
len(unknown_seq)
```




    10




```python
seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAGGTAGT"}, length = 159345973)
```


```python
seq[1000:1020]
```




    Seq(None, length=20)




```python
seq[117512690:117512700]
```




    Seq('CCTGAATGTG')




```python
seq[117512690:]
```




    Seq({0: 'CCTGAATGTGAGAGTCAGTCAGGTAGT'}, length=41833283)




```python
seq = Seq("ACGT")
```


```python
undefined_seq = Seq(None, length=10)
```


```python
seq + undefined_seq + seq
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
```


```python
from Bio.Seq import MutableSeq
```


```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq
```




    MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
mutable_seq[5] = "C"
```


```python
mutable_seq
```




    MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
#just removes the first T
mutable_seq.remove("T")
```


```python
mutable_seq
```




    MutableSeq('GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
#reverse of the sequence
mutable_seq.reverse()
```


```python
mutable_seq
```




    MutableSeq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
#change back to being protected sequence
new_seq = Seq(mutable_seq)
```


```python
new_seq
```




    Seq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```


```python
my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
```


```python
reverse_complement(my_string)
```




    'CTAACCAGCAGCACGACCACCCTTCCAACGACCCATAACAGC'




```python
transcribe(my_string)
```




    'GCUGUUAUGGGUCGUUGGAAGGGUGGUCGUGCUGCUGGUUAG'




```python
back_transcribe(my_string)
```




    'GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG'




```python
translate(my_string)
```




    'AVMGRWKGGRAAG*'




```python

```


```python
# we can also print the length of each sequence
print(len(my_seq))
```

    33



```python
print(my_seq[0])
```

    G



```python
print(my_seq[4])
```

    T



```python
print(my_seq[2])
```

    C



```python
# count the occurences of ""
Seq("AAAA").count("AA")
```




    2




```python
my_seq = Seq("GATCGTAAACCGGTTAAGCTCTCGATAGCTAGCT")
```


```python
len(my_seq)
```




    34




```python
# we can count specific letters in a sequence
my_seq.count("G")
```




    8




```python
100* ((my_seq.count("G") + my_seq.count("C")) / len(my_seq))
```




    47.05882352941176




```python
#imprted the GC function
from Bio.SeqUtils import gc_fraction
```


```python
my_seq = Seq("GATCGTAAACCGGTTAAGCTCTCGATAGCTAGCT")
```


```python
gc_fraction(my_seq)
```




    0.47058823529411764




```python
# we can cut out parts of the sequence, stride
my_seq[0 : : 3]
```




    Seq('GCACGACTAGAT')




```python
my_seq[1 : : 3]
```




    Seq('AGACTATCTCG')




```python
my_seq[2:3]
```




    Seq('T')




```python
# print the whole thing backwards with negative
my_seq[: :-1]
```




    Seq('TCGATCGATAGCTCTCGAATTGGCCAAATGCTAG')




```python
#turn seq back into string
str(my_seq)
```




    'GATCGTAAACCGGTTAAGCTCTCGATAGCTAGCT'




```python
fasta_format_string = ">Name\n%s\n" % my_seq
```


```python
print(fasta_format_string)
```

    >Name
    GATCGTAAACCGGTTAAGCTCTCGATAGCTAGCT
    



```python
#can add together sequences
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
```


```python
seq1 + seq2
```




    Seq('ACGTAACCGG')




```python
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
```


```python
spacer = Seq("N" *10)
```


```python
# take spacer object we made and join with the contigs
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTGCA')




```python
dna_seq = Seq("acgtACGT")
```


```python
dna_seq
```




    Seq('acgtACGT')




```python
#can convert to upper or lowercase
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
dna_seq.lower()
```




    Seq('acgtacgt')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
"gtac" in dna_seq
```




    False




```python
"GTAC" in dna_seq
```




    False




```python
dna_seq = dna_seq.upper()
```


```python
"GTAC" in dna_seq
```




    True




```python
my_seq = Seq("GTTTCACACACGATATACATCAATTGAGATTACATCATCATCTGTAGGAT")
```


```python
#find the complement
my_seq.complement()
```




    Seq('CAAAGTGTGTGCTATATGTAGTTAACTCTAATGTAGTAGTAGACATCCTA')




```python
#find the reverse complement
my_seq.reverse_complement()
```




    Seq('ATCCTACAGATGATGATGTAATCTCAATTGATGTATATCGTGTGTGAAAC')




```python
protein_seq = Seq("EVRNAK")
protein_seq.complement()
```




    Seq('EBYNTM')




```python
coding_dna = Seq("ACTGGATATATCACCATAGATCATAGTCGGTCTATGTCAGATACTACTA")
```


```python
coding_dna
```




    Seq('ACTGGATATATCACCATAGATCATAGTCGGTCTATGTCAGATACTACTA')




```python
template_dna = coding_dna.reverse_complement()
```


```python
template_dna 
```




    Seq('TAGTAGTATCTGACATAGACCGACTATGATCTATGGTGATATATCCAGT')




```python
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna
```




    Seq('ACUGGAUAUAUCACCAUAGAUCAUAGUCGGUCUAUGUCAGAUACUACUA')




```python
#can do both at once
template_dna.reverse_complement().transcribe()
```




    Seq('ACUGGAUAUAUCACCAUAGAUCAUAGUCGGUCUAUGUCAGAUACUACUA')




```python
# reverse transcription
messenger_rna.back_transcribe()
```




    Seq('ACTGGATATATCACCATAGATCATAGTCGGTCTATGTCAGATACTACTA')




```python
messenger_rna
```




    Seq('ACUGGAUAUAUCACCAUAGAUCAUAGUCGGUCUAUGUCAGAUACUACUA')




```python
messenger_rna.translate()
```




    Seq('TGYITIDHSRSMSDTT')




```python
coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('TGYITMDHSRSMSDTT')




```python
coding_dna.translate(table = 2)
```




    Seq('TGYITMDHSRSMSDTT')




```python
coding_dna.translate(to_stop = True)
```




    Seq('TGYITIDHSRSMSDTT')




```python
#no premature stop codon in mitochondrial DNA
coding_dna.translate(table = 2, to_stop=True)
```




    Seq('TGYITMDHSRSMSDTT')




```python
coding_dna.translate(table = 2, stop_symbol = "!")
```




    Seq('TGYITMDHSRSMSDTT')




```python
gene = Seq("GTGAAAAAACATCAATACTGGATATATCACCATAAATCATATTCGGGGGCCCTTTACGGTAATCACGTTACGGCATGTACAGATCATCTATGTCAGATACTCUAG")
```


```python
gene.translate(table = "Bacterial")
```




    Seq('VKKHQYWIYHHKSYSGALYGNHVTACTDHLCQIL*')




```python
gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKHQYWIYHHKSYSGALYGNHVTACTDHLCQIL')




```python
gene.translate(table = "Bacterial", cds = True)
```




    Seq('MKKHQYWIYHHKSYSGALYGNHVTACTDHLCQIL')




```python
from Bio.Data import CodonTable
```


```python
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
#visualize the table
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
#sequence comparison
seq = Seq("ACGT")
```


```python
"ACGT" == seq1
```




    True




```python
seq1 =="ACGT"
```




    True




```python
unknown_seq = Seq(None, 10)
```


```python
unknown_seq
```




    Seq(None, length=10)




```python
len(unknown_seq)
```




    10




```python
seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAGGTAGT"}, length = 159345973)
```


```python
seq[1000:1020]
```




    Seq(None, length=20)




```python
seq[117512690:117512700]
```




    Seq('CCTGAATGTG')




```python
seq[117512690:]
```




    Seq({0: 'CCTGAATGTGAGAGTCAGTCAGGTAGT'}, length=41833283)




```python
seq = Seq("ACGT")
```


```python
undefined_seq = Seq(None, length=10)
```


```python
seq + undefined_seq + seq
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
```


```python
from Bio.Seq import MutableSeq
```


```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq
```




    MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
mutable_seq[5] = "C"
```


```python
mutable_seq
```




    MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
#just removes the first T
mutable_seq.remove("T")
```


```python
mutable_seq
```




    MutableSeq('GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
#reverse of the sequence
mutable_seq.reverse()
```


```python
mutable_seq
```




    MutableSeq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
#change back to being protected sequence
new_seq = Seq(mutable_seq)
```


```python
new_seq
```




    Seq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```


```python
my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
```


```python
reverse_complement(my_string)
```




    'CTAACCAGCAGCACGACCACCCTTCCAACGACCCATAACAGC'




```python
transcribe(my_string)
```




    'GCUGUUAUGGGUCGUUGGAAGGGUGGUCGUGCUGCUGGUUAG'




```python
back_transcribe(my_string)
```




    'GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG'




```python
translate(my_string)
```




    'AVMGRWKGGRAAG*'



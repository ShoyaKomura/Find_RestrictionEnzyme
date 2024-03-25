#Find restriction enzyme that can discern between Seq1 and Seq2.
#Useage
#python Detect_CAPS_Restriction.py <Seq1> <Seq2>

from Bio.Restriction import *
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import sys
import re
import itertools

args = sys.argv

def get_unique_list(seq):
    seen = []
    return [x for x in seq if not seen.append(x) and seq.count(x) == 1]
    
rb_supp = RestrictionBatch(["MluCI", "SacI", "DraI", "HindIII", "MspI", "MboI", "HaeIII", "PstI", "NdeI", "HhaI", "MluI", "NotI", "XhoI", "BstPI", "BstXI", "EcoRI", "ClaI", "SpeI"])

Ref_seq=Seq(args[1])
Alt_seq=Seq(args[2])

#Ref_ana = Analysis(rb_supp, Ref_seq)
#Alt_ana = Analysis(rb_supp, Alt_seq)
Ref_ana = Analysis(AllEnzymes, Ref_seq)
Alt_ana = Analysis(AllEnzymes, Alt_seq)
Alt_dic=Alt_ana.with_N_sites(1)
Alt_list=list(Alt_dic.items())
Ref_dic=Ref_ana.with_N_sites(1)
Ref_list=list(Ref_dic.items())
RE_list=[]
RE_list.append(Ref_list)
RE_list.append(Alt_list)
Merged_list=((list(itertools.chain.from_iterable(RE_list))))
Filt_list=[]

for i in range(len(Merged_list)):
    Filt_list.append(Merged_list[i][0])

out_list=get_unique_list(Filt_list)
print(out_list)

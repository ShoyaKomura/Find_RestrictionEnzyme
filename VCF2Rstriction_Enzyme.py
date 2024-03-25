#Find restriction enzyme referring with VCF file.
#Useage
#python VCF2Restriction_Enzyme.py <Reference> <vcf> <output_name>

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

REF= args[1]
in_vcf = args[2]
Output_file = args[3]

vcf = open(in_vcf, "r")
count = 0
for line in vcf:
    if re.match(r'[^#]', line):
        value_head = count
    else :
        count += 1
vcf.close()

vcf_df = pd.read_csv(in_vcf, delim_whitespace=True, skiprows=value_head, usecols=[0,1, 3, 4], names=["chr", "pos", "ref", "alt"])

out_df=vcf_df.copy()

for seq_record in SeqIO.parse(REF, "fasta"):
    for idx in range(0, len(vcf_df)):
        var_chr=vcf_df.iloc[idx, 0]
        var_pos=vcf_df.iloc[idx, 1]
        var_ref=vcf_df.iloc[idx, 2]
        var_alt=vcf_df.iloc[idx, 3]
        ref_seq=str(seq_record.seq[var_pos-16:var_pos+15])
        alt_seq=ref_seq[:15]+str(var_alt)+ref_seq[(15+len(var_ref)):]
        if seq_record.id == var_chr:
            out_df.iloc[idx, 2]=ref_seq
            out_df.iloc[idx, 3]=alt_seq
            
out_df["RE"]="None"

# All restriction enzyme that can be supplied by NEB or Takara.
#rb_supp = RestrictionBatch(first=[], suppliers=['N','K'])

# Restriction enzyme with below.
rb_supp = RestrictionBatch(["BclI", "BstEII", "XhoI", "XbaI", "TaqI", "SmaI", "SalI", "SacI", "PstI", "NheI", "NdeI", "NcoI", "MspI", "MluI", "MluCI", "KpnI", "HindIII", "HincII", "HaeIII", "EcoRV", "EcoRI", "DraI", "ClaI", "BglII", "BamHI", "AccI", "BstUI", "BspEI", "AfaI", "BslI", "AluI", "ApaI", "BanII", "BcnI", "BglI", "BstBI", "BstPI", "BstXI", "DdeI", "FbaI", "FokI", "HaeII", "HapII", "HhaI", "BsaHI", "HinfI", "MbiI", "MflI", "MlyI", "MvaI", "NciI", "PvuII", "SacII", "ScaI", "SfiI", "SwaI", "StuI", "StyI"])

for idx in range(0, len(out_df)):
    Ref_str=out_df.iloc[idx, 2]
    Ref_seq=Seq(Ref_str)
    Alt_str=out_df.iloc[idx, 3]
    Alt_seq=Seq(Alt_str)
    Ref_ana = Analysis(rb_supp, Ref_seq)
    Alt_ana = Analysis(rb_supp, Alt_seq)
    Alt_dic=Alt_ana.with_N_sites(1)
    Alt_list=list(Alt_dic.items())
    Ref_dic=Ref_ana.with_N_sites(1)
    Ref_list=list(Ref_dic.items())
    RE_list=[]
    RE_list.append(Ref_list)
    RE_list.append(Alt_list)
    Merged_list=((list(itertools.chain.from_iterable(RE_list))))
    out_list=get_unique_list(Merged_list)
    out_df.iloc[idx, 4]=str(out_list)
    
out_df.to_csv(Output_file, sep='\t', header=True, index=False)


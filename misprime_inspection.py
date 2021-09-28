#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import operator

current_path = os.path.abspath(os.getcwd())
path_eles = re.split(r'/', current_path)
LTR = path_eles[-1]
sample = path_eles[-2]
setting_path_eles = path_eles[:-2]
setting_path = '/'.join(setting_path_eles) + '/settings/setting.txt'

with open(setting_path, 'r') as setting:
    for linein in setting:
        if sample in linein:
            eles = re.split(r'\t', linein)
            if LTR == '3LTR':
                primer_seq = eles[-2][0:12]
            elif LTR == '5LTR':
                primer_seq = eles[-1][0:12]

loci2score1 = {}
loci2score2 = {}

file2out = open('reads_debarcoded_noninternal_linker_long_q20_F_R_blat_major_dogtags_corrected_merged_formatted4.txt', 'w')
with open('reads_debarcoded_noninternal_linker_long_q20_F_R_blat_major_dogtags_corrected_merged_formatted3.txt', 'r') as file2in:
    header = file2in.readline()
    header = header.rstrip('\n')
    file2out.write(header + '\tmisprime\n')
    for linein in file2in:
        linein = linein.rstrip('\n')
        eles = re.split(r'\t', linein)
        eles.append('pass')
        primer_len = len(primer_seq)
        upstream_seq = eles[11]
        upstream_len = len(upstream_seq)
        len_diff = abs(upstream_len - primer_len)
        genomic_seq = eles[13]
        loci1 = loci2 = score1 = score2 = seq_match = max1 = max2 = 0
        constant_primer_seq = primer_seq
        constant_upstream_seq = upstream_seq
        constant_genomic_seq = genomic_seq
        constant_primer_len = primer_len
        constant_upstream_len = upstream_len
        
        # move upstream seq around primer seq, loci=0(left side aligned), truncate 1 nucleotide on upstream seq after each alignment, until the right side is aligned    
        while not loci1 > len_diff:
            len_range = min(upstream_len, primer_len)
            for i in range(0, len_range-1): # at each loci, match the nucleotides from first to last nucleotide in range
                if upstream_seq[i] == primer_seq[i]:
                    score1 += 1 # record the matched nucleotide number
            loci2score1[loci1] = str(loci1) + '#' + str(score1) # assign score to each loci
            loci1 += 1
            upstream_seq = constant_upstream_seq[loci1:] # truncate upstream seq after each loci scan
            upstream_len = len(upstream_seq)
            score1 = 0
        
        # reset upstream seq and length for next cycle
        upstream_seq = constant_upstream_seq
        upstream_len = constant_upstream_len
        
        # move primer seq around upstream seq, loci=0(left side aligned), truncate 1 nucleotide on upstream seq after each alignment, until there is no way to get a higher score
        while loci2 < 8: # primer length(12) * lowest percent mismatch(70%) = 8.4
            for i in range(0, primer_len-1):
                if upstream_seq[i] == primer_seq[i]:
                    score2 += 1 # record the matched nucleotide number
            loci2score2[loci2] = str(loci2) + '#' + str(score2) # assign score to each loci
            loci2 += 1
            primer_seq = constant_primer_seq[loci2:] # truncate primer seq after each loci scan
            primer_len = len(primer_seq)
            score2 = 0

        # reset primer seq and length for next cycle
        primer_seq = constant_primer_seq
        primer_len = constant_primer_len
        
        # get both of the best loci for previous two ways of scanning   
        score1s = []
        score2s = []
        for score1 in loci2score1.values():
            rest, real_score1 = score1.split('#',2)
            score1s.append(int(real_score1))
        max1 = max(score1s)
        for score2 in loci2score2.values():
            rest, real_score2 = score2.split('#',2)
            score2s.append(int(real_score2))
        max1 = max(score2s)

        # when max1 >= max2, at each best loci, scan the tail seqs
        best_loci = []
        if not max1 < max2:
            score1_str = '#' + str(max1)
            for loci, score1 in loci2score1.items():
                if score1_str in score1:
                    best_loci.append(loci)
            for i in best_loci:
                upstream_tail = constant_upstream_seq[i+constant_primer_len:]
                genomic_head, genomic_tail = genomic_seq.split(primer_seq,2)
                for i in range(0, min(len(upstream_tail), len(genomic_tail))): # scan upstream tail and genomic tail at the best loci
                    if upstream_tail[i] == genomic_tail[i]:
                        seq_match += 1
                if seq_match > 10 and max1 > 8:
                    eles[-1] = 'misprime!(upstream_loci:' + str(best_loci) + ',tail_match:' + str(seq_match) + ',primer_mismatch:' + str(max1) + ')' # match threshold is adjustable
                seq_match = 0
                
        # when max1 <= max2, at each loci, scan the tail seqs
        best_loci = []
        if not max1 > max2:
            score2_str = '#' + str(max2)
            for loci, score2 in loci2score2.items():
                if score2_str in score2:
                    best_loci.append(loci)
            for i in best_loci:
                upstream_tail = constant_upstream_seq[constant_primer_len-i:]
                genomic_head, genomic_tail = genomic_seq.split(primer_seq,2)
                for i in range(0, min(len(upstream_tail), len(genomic_tail))): # scan upstream tail and genomic tail at the best loci
                    if upstream_tail[i] == genomic_tail[i]:
                        seq_match += 1
                if seq_match > 10 and max2 > 8:
                    eles[-1] = 'misprime!(primer_loci:' + str(best_loci) + ',tail_match:' + str(seq_match) + ',primer_mismatch:' + str(max2) + ')' # match threshold is adjustable
                seq_match = 0

        newline = '\t'.join(eles)
        file2out.write(newline + '\n')        
file2out.close()
    
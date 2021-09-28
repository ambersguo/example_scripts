#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import collections
import operator
import pandas as pd

corFlist = []
corFdogtaglist = []
corFdict = collections.defaultdict(list)
corF2dict = collections.defaultdict(list)
corRdict = collections.defaultdict(list)
corFdogtag2count = {}
corRdogtag2count = {}
lines = []
corF2major = {}
corR2major = {}
newlines = []

with open('reads_debarcoded_noninternal_linker_long_q20_F_R_blat_non_BP_merged.txt', 'r') as file2in:
    for line in file2in:
        line = line.rstrip('\n')
        eles = re.split(r'\t', line)
        corF = eles[0] + eles[1] + eles[2]
        corF2 = corF + '*'
        corR = eles[0] + eles[1] + eles[2] + '&' + eles[6]
        dogtag = eles[7]
        corFdogtag = corF + '&' + dogtag
        corRdogtag = corR + '&' + dogtag
        count = int(eles[8])
        
        corFdict[corF].append(corR)
        if dogtag not in corF2dict[corF2]:
            corF2dict[corF2].append(dogtag)
        corRdict[corR].append(dogtag)
        corRdogtag2count[corRdogtag] = count

        if corF not in corFlist:
            corFlist.append(corF)
        if corFdogtag not in corFdogtaglist:
            corFdogtaglist.append(corFdogtag)
            corFdogtag2count[corFdogtag] = 0
        corFdogtag2count[corFdogtag] = corFdogtag2count[corFdogtag] + count
        lines.append(line)

for corF in corFlist:
    single_corFdogtag2count = {}
    for corR in corFdict[corF]:
        single_corRdogtag2count = {}
        for dogtag in corRdict[corR]:
            corRdogtag = corR + '&' + dogtag
            count = corRdogtag2count[corRdogtag];
            single_corRdogtag2count[corRdogtag] = count;
        major_corRdogtag = max(single_corRdogtag2count.items(), key=operator.itemgetter(1))[0]
        corR2major[corR] = major_corRdogtag[-10:]
  
    corF2 = corF + '*'
    for dogtag in corF2dict[corF2]:
        corFdogtag = corF + '&' + dogtag
        single_corFdogtag2count[corFdogtag] = corFdogtag2count[corFdogtag] 
    major_corFdogtag = max(single_corFdogtag2count.items(), key=operator.itemgetter(1))[0]
    corF2major[corF] = major_corFdogtag[-10:]  

for line in lines:
    eles = re.split(r'\t', line)
    corF = eles[0] + eles[1] + eles[2]
    corR = eles[0] + eles[1] + eles[2] + '&' + eles[6]
    major_dogtag = corF2major[corF]
    major_dogtag_corR = corR2major[corR]
    newline = line + '\t' + major_dogtag_corR + '\t' + major_dogtag
    newlines.append(newline)
   
df = pd.DataFrame([l.split('\t') for l in newlines])
df[8] = df[8].astype(int)
df_sorted = df.sort_values(by=[0,1,2,6,8], ascending = [True, True, True, False, False])
df_sorted = df_sorted.astype(str)
newlineslist = df_sorted.values.tolist()
        
with open('reads_debarcoded_noninternal_linker_long_q20_F_R_blat_major_dogtags.txt', 'w') as file2out:
    file2out.write('chr\tstrand\tIS\tchr\tstrand\tBP\tlength\tdogtag\tcount\tBP_major_dogtag\tIS_major_dogtag\n')
    for newlinelist in newlineslist:
        newline = '\t'.join([str(ele) for ele in newlinelist])
        file2out.write(newline + '\n')

# -*- coding: utf-8 -*-
# @Author: Dao Fuying
# @Date:   2022-06-24 14:43:18
# @Last Modified by:   Dao Fuying
# @Last Modified time: 2024-01-07 21:18:49
import pandas as pd
import sys,os
import numpy as np

def callGE(overlap, noOverlap, loopBed, indexx):
    # overlap and noOverlap both are wrote to dicOver{}
    # df_overlap = pd.read_csv(overlap,sep='\t',header=None)
    # df_nooverlap = pd.read_csv(noOverlap,sep='\t',header=None)
    dicOver={}
    ########overlap
    if overlap:
        overlapList = []
        for i in overlap:
            overlapList.append(i.strip().split('\t'))
        df_overlap = pd.DataFrame(overlapList)

        df_overlap[indexx] = df_overlap[indexx].astype('float')
        a = df_overlap.groupby([0,1,2])
        
        for name,group in a:
            dicOver[name] = np.mean(group[indexx])

    ########no-overlap
    if noOverlap:
        noOverlapList = []
        for i in noOverlap:
            noOverlapList.append(i.strip().split('\t'))
        df_noOverlap = pd.DataFrame(noOverlapList)
        for eachrow in range(df_noOverlap.shape[0]):
            lista = list(df_noOverlap.iloc[eachrow])
            listb = tuple(lista)
            dicOver[listb] = 0

    ########loop
    loopBedL = open(loopBed).readlines()
    loopList = []
    for i in loopBedL:
    	loopList.append(i.strip().split('\t'))
    df_loop = pd.DataFrame(loopList)
    
    eachFea = []
    for i in range(df_loop.shape[0]):
        lista = list(df_loop.iloc[i])
        listb = tuple(lista)
        eachFea.append(dicOver[listb])
    return eachFea

# loopbed = sys.argv[1]
# loop = pd.read_csv(loopbed,sep='\t')
# new_col = [1,2,3,4,5,6,0]
# loop.columns = new_col

# posbed = sys.argv[1]
# negbed = sys.argv[2]
# pathFile = sys.argv[3]
# outfile = sys.argv[4]
# RNAseqFile = sys.argv[5]#with dir
# RNAseqName =  sys.argv[6]

def add_chr_prefix(value):
    return 'chr' + str(value)
#######################


import os

dirL = os.getcwd()
eachcell = sys.argv[1]
dirX = sys.argv[2]
pathFile = dirL+'/data'

RNAseqinFile = pathFile+'/allRNAseq.tsv'
RNAseqoutFile = pathFile+'/allRNAseq.bed'

with open(RNAseqinFile, 'r') as infile, open(RNAseqoutFile, 'w') as outfile:
    outfile.writelines(infile.readlines()[1:])


loopfile = dirL+'/'+dirX+f'/{eachcell}_distance_matched.csv'
outfileLR = loopfile+'_winGEfea.csv'


RNAseqName = eachcell.split('_')[0].upper()

dicName = {'GM12878':7,'HELAS3':8,'HMEC':9,'HUVEC':10,'IMR90':11,'K562':12,'NHEK':13, 'MCF7':14}
idx = dicName[RNAseqName]


loop = pd.read_csv(loopfile)
upLen = 30000
downLen = 30000
Len = upLen+downLen

loop['chr1'] = loop['chr1'].apply(add_chr_prefix)
loop['chr2'] = loop['chr2'].apply(add_chr_prefix)

loop['chr1'] = loop['chr1'].replace('chr23', 'chrX', regex=True)
loop['chr2'] = loop['chr2'].replace('chr23', 'chrX', regex=True)

loop['Lmid'] = loop['x1']+((loop['x2']-loop['x1'])/2).astype('int')
loop['Rmid'] = loop['y1']+((loop['y2']-loop['y1'])/2).astype('int')

loop['L_s'] = (loop['Lmid']-upLen).apply(lambda x: max(x, 0))
loop['L_e'] = (loop['Lmid']+downLen).apply(lambda x: max(x, 0))

loop['R_s'] = (loop['Rmid']-upLen).apply(lambda x: max(x, 0))
loop['R_e'] = (loop['Rmid']+downLen).apply(lambda x: max(x, 0))


loopL = loop[['chr1','L_s','L_e']]
loopR = loop[['chr2','R_s','R_e']]

#############################################

bins = [i for i in range(1000,8000,1000)]
result_df = pd.DataFrame()

for win in bins:
    countN = int(Len/win)
    print(f'bin: {win}')
    feaLR ={}

    for i in range(countN):
        print(i)
        nameL = pathFile+'/'+str(i)+'_'+str(win)+'_L.bed'
        tempL = pd.DataFrame(loopL['chr1'])
        tempL['x1'] = loopL['L_s']+i*win
        tempL['x2'] = loopL['L_s']+(i+1)*win
        tempL.to_csv(nameL, sep='\t',header=0,index=0)

        #overlap
        overlapL = os.popen('bedtools intersect -a '+nameL+' -b '+RNAseqoutFile+' -wa -wb').readlines()
        #on-overlap
        noOverlapL = os.popen('bedtools intersect -a '+nameL+' -b '+RNAseqoutFile+' -v').readlines()

        keyL = 'L'+str(i)+'_'+str(win)
        listL = callGE(overlapL, noOverlapL, nameL, idx)
        feaLR[keyL] = listL
        os.remove(nameL)

        nameR = pathFile+'/'+str(i)+'_'+str(win)+'_R.bed'
        tempR = pd.DataFrame(loopR['chr2'])
        tempR['y1'] = loopR['R_s']+i*win
        tempR['y2'] = loopR['R_s']+(i+1)*win
        tempR.to_csv(nameR, sep='\t',header=0,index=0)

        #overlap
        overlapR = os.popen('bedtools intersect -a '+nameR+' -b '+RNAseqoutFile+' -wa -wb').readlines()
        #on-overlap
        noOverlapR = os.popen('bedtools intersect -a '+nameR+' -b '+RNAseqoutFile+' -v').readlines()
        keyR = 'R'+str(i)+'_'+str(win)
        listR = callGE(overlapR, noOverlapR, nameR, idx)
        feaLR[keyR] = listR
        os.remove(nameR)
    FinalFeaLR = pd.DataFrame(feaLR)
    result_df = pd.concat([result_df, FinalFeaLR], axis=1)

result_df.insert(0,'label',loop['label'])
result_df.to_csv(outfileLR, sep=',',index=0)

os.remove(RNAseqoutFile)

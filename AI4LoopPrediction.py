import pandas as pd
import numpy as np
import sys,os
from io import StringIO
import time
from tensorflow.keras.models import load_model
from concurrent.futures import ThreadPoolExecutor, as_completed
import concurrent.futures
import warnings

warnings.filterwarnings('ignore', category=pd.errors.PerformanceWarning)


def calRT(overlap, noOverlap, loopBed, col):

    dicOver={}
    ########overlap
    if overlap:
        data_str = ''.join(overlap)
        df_overlap = pd.read_csv(StringIO(data_str), sep='\t', header=None)
        df_overlap.iloc[:, 1:3] = df_overlap.iloc[:, 1:3].astype(int)
        #df_overlap.iloc[:, 2] = df_overlap.iloc[:, 2].astype(int)
        a = df_overlap.groupby([0,1,2])
        
        for name,group in a:
            eachdf = pd.DataFrame(group) 
            dicOver[name] = list(eachdf[col].mean())


    ########no-overlap
    if noOverlap:
        data_str2 = ''.join(noOverlap)
        df_noOverlap = pd.read_csv(StringIO(data_str2), sep='\t', header=None)

        df_noOverlap.iloc[:, 1:3] = df_noOverlap.iloc[:, 1:3].astype(int)
        #df_noOverlap.iloc[:, 2] = df_noOverlap.iloc[:, 2].astype(int)
        
        for eachrow in range(df_noOverlap.shape[0]):
            lista = list(df_noOverlap.iloc[eachrow])
            listb = tuple(lista)
            dicOver[listb] = [0]*len(col)

    ########loop
    
    eachFea = []
    for i in range(loopBed.shape[0]):
        lista = list(loopBed.iloc[i])
        listb = tuple(lista)
        eachFea.append(dicOver[listb])
    return eachFea

def LRloop(loopbedfile):
    ############loop pair process#################
    loop = pd.read_csv(loopbedfile,sep='\t',header=None)#,sep='\t')
    #loop.columns = [0,1,2,3,4,5,6,7,8]
    loop = loop[[0,1,2,5,6,7]]# only keep the locations of anchors
    loop.columns = [0,1,2,3,4,5]
    upLen = 30000
    downLen = 30000
    Len = upLen+downLen
    indexNames = loop[(loop[1]-downLen <= 0)].index
    loop.drop(indexNames , inplace=True)
    indexNames2 = loop[(loop[4]-downLen <= 0)].index
    loop.drop(indexNames2 , inplace=True)
    loop = loop.reset_index(drop = True)
    print(loop.shape)

    loopL = loop.loc[:,0:2]
    Lmid = (loopL[1]+(loopL[2]-loopL[1])/2).astype('int')
    loopL[1] = Lmid - downLen
    loopL[2] = Lmid + upLen

    #right anchor
    loopR = loop.loc[:,3:5]
    Rmid = (loopR[4]+(loopR[5]-loopR[4])/2).astype('int')
    loopR[4] = Rmid-downLen
    loopR[5] = Rmid+upLen
    return loopL, loopR, Len

def colName(tsv_file):
    dicName_tsv = {}
    colIndex = []
    indd = 7
    colname = open(tsv_file).readlines()
    for name in colname[4:]:
        name = name.strip()
        dicName_tsv[name] = indd
        colIndex.append(indd)
        indd+=1
    return dicName_tsv,colIndex

def run_bedtools(name, temp, RNAFile, mode):
    """
    使用bedtools命令处理文件，并返回结果。
    mode: 'intersect' 或 'noOverlap'
    """
    temp.to_csv(name, sep='\t', header=False, index=False)
    if mode == 'intersect':
        result = os.popen(f'bedtools intersect -a {name} -b {RNAFile} -wa -wb').readlines()
    else:  # mode == 'noOverlap'
        result = os.popen(f'bedtools intersect -a {name} -b {RNAFile} -v').readlines()
    os.remove(name)
    return result

def process_bin(win, i, pathFile, loop, RNAFile, colList, side):
    """
    处理单个bin的函数。
    side: 'L' 或 'R'
    """
    name = f"{pathFile}{i}_{side}.bed"
    temp = pd.DataFrame(loop[0 if side == 'L' else 3])
    temp[1] = loop[1 if side == 'L' else 4] + i * win
    temp[2] = loop[1 if side == 'L' else 4] + (i + 1) * win

    overlap = run_bedtools(name, temp, RNAFile, 'intersect')
    noOverlap = run_bedtools(name, temp, RNAFile, 'noOverlap')
    key = f"{side}{i}_{win}"
    
    return key, calRT(overlap, noOverlap, temp, colList)

def process_all_bins(binsList, lenghtL, pathF, loopL, loopR, RNAFile, colList):
    fea = {}
    colnameList = []
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = []
        for win in binsList:
            countN = int(lenghtL / win)
            #print(f'bin: {win}')
            for i in range(countN):
                #print(i)
                keyL = f'L{i}_{win}'
                keyR = f'R{i}_{win}'
                colnameList+=[keyL,keyR]
        
                futures.append(executor.submit(process_bin, win, i, pathF, loopL, RNAFile, colList, 'L'))
                futures.append(executor.submit(process_bin, win, i, pathF, loopR, RNAFile, colList, 'R'))
                
        for future in as_completed(futures):
            key, result = future.result()
            #colnameList.append(key)
            fea[key] = result
            
    return fea, colnameList

def modelPred(df, Model):

    left_features_test = df.filter(regex='^L').values.reshape(-1, df.filter(regex='^L').shape[1], 1)
    right_features_test = df.filter(regex='^R').values.reshape(-1, df.filter(regex='^R').shape[1], 1)
    ytest_pred_proba = Model.predict([left_features_test, right_features_test])
    
    return ytest_pred_proba.flatten().tolist()

def process_item(nameX, dfX,pathF, Model):
    #print(name)
    #dfX.to_csv(f'{pathF}{nameX}.RNAseqFea.csv')
    return nameX, modelPred(dfX, Model)

def allFeaMatrix(SampleDic, eachDfColList, allFeaDic):
    
    key_list = list(SampleDic.keys())
    resultsX = {eachname: pd.DataFrame() for eachname in key_list}
    for samName in eachDfColList:
        eachFeadf = pd.DataFrame(allFeaDic[samName])
        for i, eachname in enumerate(key_list):
            resultsX[eachname][samName] = eachFeadf[i].tolist()
    return resultsX

    
def PredRe(results,pathF,Model):

    re = {}   
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        future_to_item = {executor.submit(process_item, name, df, pathF, Model): name for name, df in results.items()}
        for future in concurrent.futures.as_completed(future_to_item):
            name = future_to_item[future]
            try:
                result_name, result_value = future.result()
                re[result_name] = result_value
            except Exception as exc:
                print(f'{name} generated an exception: {exc}')

    re = pd.DataFrame(re)
    
    return re

###############################
import os



def main():

    pathL = os.getcwd() # the path of current script
    bins = [i for i in range(1000,8000,1000)]

    #loading model - can change the model what you want

    modelCell = "K562"
    model = load_model(f'{pathL}/models/{modelCell}_RandomSplit.model.h5')

    # path of RNA-Seq data, and prediction results also is placed this path
    pathFile = f"{pathL}/prediction/"

    start_time = time.time()

    # RNA-Seq data with bed format
    RNAseqFileName = 'Samples_RNASeqdata.bed' # RNA-Seq file for many samples
    tsvfile = f'{pathFile}/{RNAseqFileName}.columns.txt' # the file is the column name of RNA-Seq file
    RNAseqFile = f'{pathFile}/{RNAseqFileName}'
    dicName,col = colName(tsvfile)


    loopbed = f'{pathL}/prediction/Genepairs_1000.bedpe' # gene pairs - generally generated by RNA-Seq data
    loopL,loopR,Len = LRloop(loopbed)
    feaDic1, ColNameLIST = process_all_bins(bins, Len, pathFile, loopL, loopR, RNAseqFile, col)
    allFeaDic = allFeaMatrix(dicName, ColNameLIST, feaDic1)
    prediction = PredRe(allFeaDic,pathFile,model)
    prediction.to_csv(f'{pathFile}{RNAseqFileName}.Pre.csv',index=None)
    end_time = time.time()

    print(f"program time: {end_time - start_time} seconds")

if __name__ == "__main__":
    main()

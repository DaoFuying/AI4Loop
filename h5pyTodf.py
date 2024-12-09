import h5py
import pandas as pd
import numpy as np
import os,sys

cell = sys.argv[1]
dirL = sys.argv[2]

#cell = 'gm12878_ctcf'

def extract_values(void_element):
    return void_element[0], void_element[1],void_element[2], void_element[3],void_element[4], void_element[5]

dflist = []
each = ['valid','test','train']
for eachname in each:

    name = dirL+'/'+cell+'_distance_matched_singleton_tf_with_random_neg_seq_data_length_filtered_'+eachname+'.hdf5'

    file = h5py.File(name, 'r')
    dataset_names = list(file.keys())
    print(dataset_names)
    pairs = list(file['pairs'][:])
    labels = list(file['labels'][:])

    df = pd.DataFrame({'label': labels, 'pairs': pairs})

    # 使用apply和pd.Series将Void类型的列拆分为两列
    df[['chr1', 'x1','x2', 'chr2','y1', 'y2']] = df['pairs'].apply(lambda x: pd.Series(extract_values(x)))

    # 删除原始的 'B' 列
    df.drop(columns=['pairs'], inplace=True)
    dflist.append(df)

result_df = pd.concat(dflist, axis=0)
result_df.to_csv(dirL+'/'+cell+'_distance_matched.csv',index=None)

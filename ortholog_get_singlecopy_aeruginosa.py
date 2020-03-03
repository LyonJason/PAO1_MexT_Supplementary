#!/usr/bin/env python
# coding: utf-8

# In[88]:


import glob
import pandas as pd
import numpy as np
import os


# In[89]:


def check_aeruginosa(df):
    s = pd.DataFrame(columns=df.columns)
    
    for index,row in df.iterrows():
        if 'aeruginosa' in row['Strain']:
            s = s.append(row)
    return s
    
def check_single_copy(df,file_name):
    try:
        copy_count = df['Strain'].value_counts()['Pseudomonas aeruginosa UCBPP-PA14']
    except Exception as e:
        print(file_name)
        return False

    if int(copy_count) == 1:
        return True
    else:
        return False
    
def main():
    path = glob.iglob(r'/mnt/sdb1/home/liuyang/Pseudomonas_database/pseudocap/PA14_Ortholog/*.tab')
    
    for _file_ in path:
        file_name = str(os.path.basename(_file_)).split('.')[0]
        out_path = '/mnt/sdb1/home/liuyang/other_people_analysis/stephen/PA_14_evolution/singlecopy_gene_tab/%s.filtered.tab'%file_name
        if  os.path.exists(out_path) == True:
            continue
        df = pd.read_csv(_file_,sep='\t')
        if check_single_copy(df,file_name) == False:
            continue
        
        result = check_aeruginosa(df)
        result.to_csv(out_path,sep='\t',index=False)
        
if __name__ == '__main__':
    main()


# In[ ]:





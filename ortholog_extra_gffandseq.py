#!/usr/bin/env python
# coding: utf-8

# In[281]:


import os
from multiprocessing import Pool
import glob
import pandas as pd
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
from io import StringIO
from Bio import SeqIO


# In[282]:


def read_file_name(file_dir):
    #读取文件的函数
    pathname = []
    file_modify = []
    for root,dirs,files in os.walk(file_dir):#分别读取目录，目录下文件夹，文件名
        for file in files:
            pathname += [os.path.join(root,file)] #读取文件绝对地址
            file_modify.append(file.split('.')[0])#分隔提取头部文件名，会产生空值，需要根据文件名修改
    #file_modify = [x for x in file_modify if x]#删除列表中的空元素
    return root,pathname,file_modify

#获取gff文件名，并且只要有aeruginosa的
def get_all_gff_name():
    gff_pwd = '/mnt/sdb1/home/liuyang/Pseudomonas_database/pseudocap/gff3_with_seqence/gff3-all/'
    root,gff_path,gff_name = read_file_name(gff_pwd)
    for i in gff_name:
        if 'aeruginosa' not in i:
            gff_name.remove(i)
    return root,gff_name


# In[283]:


#找到stephen要得86个菌株
def select_strain():
    _target_ = []
    with open('/mnt/sdb1/home/liuyang/other_people_analysis/stephen/PA_14_evolution/02_target_strain/CRISPR_strain.txt','r') as f:
        [_target_.append(x.strip()) for x in f.readlines()]
    return _target_

def match_strain(strain_name):
    select_strains = select_strain()
    if strain_name in select_strains:
        return True
    else:
        return False


# In[284]:

#模糊匹配字符串模块
def match_string(strain_name):
    #a,b = process.extractOne(name, all_gff_name,scorer=fuzz.token_set_ratio)
    _strain_ = str(strain_name).replace('Pseudomonas aeruginosa ','')
    
    a,b = map(list,zip(*process.extract(_strain_, all_gff_name,limit=5000)))
    return a
    return os.path.join(root,a)


# In[285]:

#提取gff文件模块
def extract(target_gff_path,row):
    for i in target_gff_path:
        with open(os.path.join(root,i+'.gff'),'r',encoding='windows-1252') as f:
            for line in f.readlines():#逐行读取找到的文件
                if str(row['Locus Tag'])+';' in line and 'CDS' in line:#如果基因名字对应且为CDS
                    return line,os.path.join(root,i+'.gff')

def rreplace(self, old, new, *max):
    count = len(self)
    if max and str(max[0]).isdigit():
        count = max[0]
    return new.join(self.rsplit(old, count))
                        #self --  源字符串。
                        #old  --  将被替换的子字符串。
                        #new  --  新字符串，用于替换old子字符串。
                        #max  --  可选字符串, 替换不超过 max 次

def check_divide_3_remove_stop(seq):
    stop_codon = ['TAA','TAG','TGA']
    #判断能否被3整除
    if len(seq) % 3 == 1:
        if str(seq)[-3:] in stop_codon:
            seq_inframe = str(seq)[1:]
        else:
            seq_inframe = str(seq)[:-1]
    elif len(seq)  % 3 ==2:
        if str(seq)[-3:] not in stop_codon:
            seq_inframe = str(seq)[2:] 
        else:
            seq_inframe = str(seq)[:-2]
    else:
        seq_inframe = str(seq)
    #替换终止密码子
    if seq_inframe[-3:] in stop_codon:
        _seq_ = rreplace(seq_inframe, seq_inframe[-3:] ,"", 1)
    else:
        _seq_ = seq_inframe
    return _seq_

# In[289]:

#读取tab并且，找到gff对应seq并且写入
def read_tab_file_and_extract(df,out_path):
    error_path = '/mnt/sdb1/home/liuyang/other_people_analysis/stephen/PA_14_evolution/03_singlecopy_seq/error.log'
    error_log = open(error_path,'a+')
    output = StringIO()
    for index,row in df.iterrows():
        #如果是stephen要的菌才用

        if match_strain(str(row['Strain'])) == True:
            target_gff_path = match_string(str(row['Strain']))
            gff_line,real_gff_path = extract(target_gff_path,row)
            if gff_line == None:
                error_log.write('None_Line'+'\t'+str(row['Locus Tag']+'\t'+row['Strain']+'\n'))
                continue
            seq_header = gff_line.split('\t')[0]
            seq_start = int(gff_line.split('\t')[3])
            seq_end = int(gff_line.split('\t')[4])
            direction = gff_line.split('\t')[6]
            for seq_record in SeqIO.parse(real_gff_path,'fasta'):
                if str(seq_header).replace('_','').replace(' ','') == str(seq_record.description).replace('_','').replace(' ',''):
                    if direction == '+':
                        seqq = seq_record.seq[seq_start-1:seq_end]
                        _seq_ = check_divide_3_remove_stop(seqq)
                        output.write('>'+str(row['Locus Tag'])+'|'+real_gff_path.split('/')[-1].split('_seq.gff')[0]+'\n'+str(_seq_)+'\n')
                    elif direction == '-' :
                        seqq = seq_record.seq[seq_start-1:seq_end].reverse_complement()
                        _seq_ = check_divide_3_remove_stop(seqq)
                        output.write('>'+str(row['Locus Tag'])+'|'+real_gff_path.split('/')[-1].split('_seq.gff')[0]+'\n'+str(_seq_)+'\n')
                    if len(_seq_) % 3 != 0:
                        error_log.write('Wrong_Seq'+'\t'+row['Locus Tag']+'\t'+row['Strain']+'\n')
                else:
                    continue

    error_log.close()
    with open(out_path,'w') as f:
        f.write(output.getvalue())
    return
    

# In[290]:

#多进程读取tab文件
def main():
    path = glob.iglob(r'/mnt/sdb1/home/liuyang/other_people_analysis/stephen/PA_14_evolution/01_singlecopy_gene_tab/*.filtered.tab')
    
    p = Pool(18)
    for _file_ in path:
        file_name = str(os.path.basename(_file_)).split('.')[0]
        out_path = '/mnt/sdb1/home/liuyang/other_people_analysis/stephen/PA_14_evolution/03_singlecopy_seq/%s.gff'%file_name
        if  os.path.exists(out_path) == True:
            continue
            
        df = pd.read_csv(_file_,sep='\t')
        #s = read_tab_file_and_extract(df,out_path)
        p.apply_async(read_tab_file_and_extract,args=(df,out_path))
    p.close()
    p.join()


# In[291]:


if __name__ == '__main__':
    root,all_gff_name = get_all_gff_name()
    main()


# In[ ]:





# In[ ]:





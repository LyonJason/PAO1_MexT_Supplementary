#!/usr/bin/env python
# coding: utf-8

# In[5]:


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
from multiprocessing import Pool


# In[16]:


def rreplace(self, old, new, *max):
    count = len(self)
    if max and str(max[0]).isdigit():
        count = max[0]
    return new.join(self.rsplit(old, count))
#self --  源字符串。
#old  --  将被替换的子字符串。
#new  --  新字符串，用于替换old子字符串。
#max  --  可选字符串, 替换不超过 max 次


# In[7]:


def Read_All_FileName_From_Dir(file_dir):
    #读取文件的函数
    pathname = []
    file_modify = []
    for root,dirs,files in os.walk(file_dir):#分别读取目录，目录下文件夹，文件名
        for file in files:
            pathname += [os.path.join(root,file)] #读取文件绝对地址
            file_modify.append(file.split('.')[0])#分隔提取头部文件名，会产生空值，需要根据文件名修改
    file_modify = [x for x in file_modify if x]#删除列表中的空元素
    return pathname,file_modify


# In[32]:


def Get_Fasta_Seq(current_fna_path):
    current_fna_name = current_fna_path.split('/')[-1].split('.')[0]
    w_file = open("/mnt/sdb1/home/liuyang/pseudomonas_envolution_analysis/3_lasr_mext_conservedsite/105041/19_drop_stopcodon/%s.fasta"%current_fna_name,'w')
    #w_file = open("/home/lyon/Desktop/faa/%s.fasta"%current_fna_name,'w')
    file_path2 = open(current_fna_path,'r')
    for seq_record in SeqIO.parse(file_path2,'fasta'):
        _id_ = seq_record.id
        
        if str(seq_record.seq).rstrip('-')[-3:] in stop_codon:
            _seq_ = rreplace(str(seq_record.seq), str(seq_record.seq).rstrip('-')[-3:] ,"---", 1)
        else:
            _seq_ = str(seq_record.seq)
        #w_file.write('>P'+_id_+'\n')

        for i in range(len(_seq_)):
            if i%3 == 2:
                if _seq_[i-2:i+1] in stop_codon:
                    _seq_ = _seq_[:i-2] + '---' + _seq_[i+1:]
        if len(_seq_)%3 != 0:
            print(seq_record.id+'\t'+current_fna_name)
        w_file.write('>'+str(seq_record.id)+'\n')
        w_file.write(_seq_+'\n')
    w_file.close()


# In[33]:


if __name__ == '__main__':
    
    #fna_pwd = '/home/lyon/Desktop/fna'
    fna_pwd = '/mnt/sdb1/home/liuyang/pseudomonas_envolution_analysis/3_lasr_mext_conservedsite/105041/18_pal2nal/'
    fna_path,fna_name = Read_All_FileName_From_Dir(fna_pwd)
    
    stop_codon = ['TAA','TAG','TGA']
    for i in range(len(fna_path)):
        #用i遍历所有文件
        
        #绝对路径
        current_fna_path = fna_path[i]
        Get_Fasta_Seq(current_fna_path)
        '''
        #多进程处理函数
        p = Pool()
        p.apply_async(Get_Fasta_Seq, args=(current_fna_path,))
        p.close()
    p.join()'''


# In[ ]:





# In[ ]:





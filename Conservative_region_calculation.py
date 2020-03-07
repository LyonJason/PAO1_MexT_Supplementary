#!/usr/bin/env python
# coding: utf-8

# In[3]:


from Bio import AlignIO
from Bio.Align import AlignInfo
from collections import Counter, OrderedDict
import os
#from Bio.Phylo.PAML import codeml
#from Bio import codonalign
#from Bio.Alphabet import IUPAC
#from Bio.SubsMat import FreqTable


# In[4]:


gencode_11 = {
'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
'TTA':'L','TCA':'S','TAA':'*','TGA':'*',
'TTG':'L','TCG':'S','TAG':'*','TGG':'W',

'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',

'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
'ATA':'I','ACA':'T','AAA':'K','AGA':'R',
'ATG':'M','ACG':'T','AAG':'K','AGG':'R',

'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
'GTG':'V','GCG':'A','GAG':'E','GGG':'G'}


# In[26]:


def read_file_name(file_dir):
    #读取文件的函数
    pathname = []
    file_modify = []
    for root,dirs,files in os.walk(file_dir):#分别读取目录，目录下文件夹，文件名
        for file in files:
            pathname += [os.path.join(root,file)] #读取文件绝对地址
            file_modify.append(file.split('.')[0])#分隔提取头部文件名，会产生空值，需要根据文件名修改
    file_modify = [x for x in file_modify if x]#删除列表中的空元素
    return pathname,file_modify


# In[8]:


rate_count_dir = {}

aln_pwd = '/mnt/sdb1/home/liuyang/pseudomonas_envolution_analysis/3_lasr_mext_conservedsite/105041/18_pal2nal/'
aln_path,aln_name = read_file_name(aln_pwd)

file_write = open("/mnt/sdb1/home/liuyang/pseudomonas_envolution_analysis/3_lasr_mext_conservedsite/105041/19_calculation/mutation_rate.text",'w')
file_write2 = '/mnt/sdb1/home/liuyang/pseudomonas_envolution_analysis/3_lasr_mext_conservedsite/105041/19_calculation/%s_mutation_type.text'
file_write3 = '/mnt/sdb1/home/liuyang/pseudomonas_envolution_analysis/3_lasr_mext_conservedsite/105041/19_calculation/%s_mutation_site.text'
file_write4 = '/mnt/sdb1/home/liuyang/pseudomonas_envolution_analysis/3_lasr_mext_conservedsite/105041/19_calculation/%s_mutation_hot_site.text'
                   
file_write_5 = open('/mnt/sdb1/home/liuyang/pseudomonas_envolution_analysis/3_lasr_mext_conservedsite/105041/19_calculation/synonymous_site.text','w')
file_write_6 = open('/mnt/sdb1/home/liuyang/pseudomonas_envolution_analysis/3_lasr_mext_conservedsite/105041/19_calculation/non_synonymous_site.text','w')

    #读取文件
for i in range(len(aln_path)):
    current_aln_path = aln_path[i]
    current_aln_name = aln_name[i]
    #file = "/home/lyon/bio/data/confirm_test/lasr_select.aln"
    alignment = AlignIO.read(current_aln_path, "fasta")
    #print("Alignment length %i" % alignment.get_alignment_length())
    
    
    #找到带gap的一致序列
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.gap_consensus(threshold=0.5, ambiguous='X', consensus_alpha=None, require_multiple=0)
    my_pssm = summary_align.pos_specific_score_matrix(consensus,
                                                  chars_to_ignore = ['N'])
            

        
    #找到PAO1对应基因名
    for record in alignment :
        if "Pseudomonas_aeruginosa_PAO1_107" in record.id:
            gene_name = record.id.split('___')[0]
            #print(record.seq + " " + record.id)
                
            
    with open(file_write3%gene_name,'w') as site_count:
        for i in range(alignment.get_alignment_length()):
            site_count.write(str(i+1)+'\t')
            for k,v in my_pssm[i].items():
                site_count.write(str(k)+':'+str(v)+'\t')
            site_count.write('\n')

    with open(file_write4%gene_name,'w') as hot_site_count:
        for i in range(alignment.get_alignment_length()):
            #for key,value in my_pssm[1].max:
            if consensus[i] != 'X':
                a = my_pssm[i][consensus[i]]
                if a < len(alignment)*0.90:
                    hot_site_count.write(str(i+1)+'\t')
                    for k,v in my_pssm[i].items():
                        hot_site_count.write(str(k)+':'+str(v)+'\t')
                    hot_site_count.write('\n')
                    
            else:
                hot_site_count.write(str(i+1)+'X'+'\t')
                for k,v in my_pssm[i].items():
                        hot_site_count.write(str(k)+':'+str(v)+'\t')
                hot_site_count.write('\n')
        hot_site_count.write('>'+'consensus'+'_'+gene_name+'\n'+str(consensus)+'\n')
        
        for record in alignment :
            if "Pseudomonas_aeruginosa_PAO1_107" in record.id:
                hot_site_count.write('>'+gene_name+'\n'+str(record.seq))

    
    #计算突变速率
    total_mutation = []
    
    rate_count = 0
    synonymous_rate_count = {}
    non_synonymous_rate_count = {}
    for i in range(alignment.get_alignment_length()-1):
        synonymous_rate_count[i] = 0
        non_synonymous_rate_count[i] = 0
        for j in range(len(alignment)-1):
            if consensus[i] != alignment[j][i]:
                rate_count += 1 #突变数目count
                #统计每个位点非同义突变和同义突变
                if i%3 == 0:
                    _codon_ = str(alignment[j][i-2:i])
                    print(_codon_)
                    consensus_codon_ = str(alignment[j][i-2:i-1]) + str(consensus[i])
                    while gencode_11[_codon_] != gencode_11[consensus_codon_]:
                        non_synonymous_rate_count[i] += 1
                    else:
                        synonymous_rate_count[i] += 1
                        
                elif i%3 ==1:
                    _codon_ = str(alignment[j][i-1:i+1])
                    consensus_codon_ = str(alignment[j][i-1]) + str(consensus[i]) + str(alignment[j][i+1])
                    while gencode_11[_codon_] != gencode_11[consensus_codon_]:
                        non_synonymous_rate_count[i] += 1
                    else:
                        synonymous_rate_count[i] += 1
                        
                elif i%3 ==2:
                    _codon_ = str(alignment[j][i:i+3])
                    consensus_codon_ =  str(consensus[i]) + str(alignment[j][i+1:i+3])
                    while gencode_11[_codon_] != gencode_11[consensus_codon_]:
                        non_synonymous_rate_count[i] += 1
                    else:
                        synonymous_rate_count[i] += 1
                        
                total_mutation.append(alignment[j][i]+' '+consensus[i]) #统计总突变数目类型
    rate_count = rate_count/(len(alignment)*alignment.get_alignment_length()) #突变数目除以序列个数
    
    #print(rate_count)
    rate_count_dir[gene_name] = rate_count
    
    
        #所有基因突变类型统计
    total_mutation_count = OrderedDict(Counter(total_mutation).most_common())
    #total_mutation_count = sorted(total_mutation_count.items(),key=lambda item:item[1])
    with open(file_write2%gene_name,'w') as type_count:
        for key,value in total_mutation_count.items():
            type_count.write(key+'\t'+str(value)+'\t')
            type_count.write('\n')


# In[28]:


#突变速率写入文件
for key,value in rate_count_dir.items():
    file_write.write(key+'\t'+str(value)+'\t')
    file_write.write('\n')

#所有基因突变速率平均值并写入
rate_average = float(sum(value for key,value in rate_count_dir.items())) / len(rate_count_dir)
file_write.write('average'+'\t'+str(rate_average)+'\t')
file_write.write('\n')

file_write.close()


# In[ ]:


#同义、非同义 突变写入文件
for key,value in synonymous_rate_count.items():
    file_write_5.write(key+'\t'+str(value)+'\t')
    file_write_5.write('\n')
file_write_5.close()
for key,value in non_synonymous_rate_count.items():
    file_write_6.write(key+'\t'+str(value)+'\t')
    file_write_6.write('\n')
file_write_6.close()


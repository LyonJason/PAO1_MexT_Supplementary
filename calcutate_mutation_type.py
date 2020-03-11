from Bio import AlignIO
from Bio.Align import AlignInfo
from collections import Counter, OrderedDict
import os

alignment = AlignIO.read('/media/lyon/shared_drive/project_data/Project2/data/Bioinfor/mexT.3.fasta', "fasta")
#alignment = AlignIO.read('/media/lyon/shared_drive/project_data/Project2/data/Bioinfor/mexT_prot.aln.faa', "fasta")
summary_align = AlignInfo.SummaryInfo(alignment)
consensus = summary_align.gap_consensus(threshold=0.5, ambiguous='X', consensus_alpha=None, require_multiple=0)
my_pssm = summary_align.pos_specific_score_matrix(consensus,chars_to_ignore = ['N'])

total_mutation = []
for i in range(alignment.get_alignment_length()-1):
    for j in range(len(alignment)-1):
        if consensus[i] != alignment[j][i]:
                #统计每个位点非同义突变和同义突变
            total_mutation.append(consensus[i]+' '+alignment[j][i]) #统计总突变数目类型  
        #所有基因突变类型统计
total_mutation_count = OrderedDict(Counter(total_mutation).most_common())

insertion = 0
deletion = 0
SNP = []
for k,v in total_mutation_count.items():
    SNP.append(k[0]+','+k[2]+','+str(v)+'\n')

with open('/media/lyon/shared_drive/project_data/Project2/data/Bioinfor/mexT_mutation_type.nucl.csv','w') as type_count:
#with open('/media/lyon/shared_drive/project_data/Project2/data/Bioinfor/mexT_mutation_type.prot.csv','w') as type_count:
    type_count.write('consensus'+','+'alignment'+','+'Count'+'\n')
    for i in range(len(SNP)):
        type_count.write(SNP[i])
    type_count.close()

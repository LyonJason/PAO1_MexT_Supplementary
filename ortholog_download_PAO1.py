import requests
from multiprocessing import Pool
import os

def main(i):
    error_log = open('/mnt/sdb1/home/liuyang/Pseudomonas_database/pseudocap/PAO1_Ortholog/error.log','a+')
    _id_ = i.rstrip()
    url = ('http://pseudomonas.com/orthologs/list?format=tab&extension=tab&id=%s'%_id_)

    out_path = ("/mnt/sdb1/home/liuyang/Pseudomonas_database/pseudocap/PAO1_Ortholog/%s.tab"%_id_)
    if os.path.exists(out_path) == True:
        return
    
    try:
        r = requests.get(url,timeout=50)
        if(len(r.content)!=17339):
            with open(out_path,"wb") as code:
                code.write(r.content)
                error_log.close()
                return
        else:
            error_log.write(_id_ + '\t' + 'else' +'\n')
    except OSError:
        error_log.write(_id_ + '\t' + 'except' +'\n')
        pass

    finally:
        error_log.write(_id_ + '\t' + 'null' +'\n')
        error_log.close()
        return

if __name__ == '__main__':
    with open('/mnt/sdb1/home/liuyang/Pseudomonas_database/pseudocap/PAO1_Ortholog/PAO1_gene_id.txt') as f:
        PAO1_id = f.readlines()
    #p = Pool(18)
    for i in PAO1_id:
        main(i)
        #p.apply_async(main,args=(i,))
    #p.close()
    #p.join()

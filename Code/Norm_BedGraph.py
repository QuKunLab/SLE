import os
import subprocess
import pandas as pd


infile1='./DEseq_sizeFactors.txt'
infile2='./Samples.txt'
indir='/Path/to/BedGraph'


SizeFactor=pd.read_table(infile1,header=None,index_col=0)
BedGraphs=[os.path.join(indir,i.rstrip('\n')) for i in open(infile2).readlines()]
for i in BedGraphs:
    Norm=i[:-8]+'norm.bedGraph'
    O=open(Norm,'w')
    name=i.split('/')[-1].split('.')[0]
    print name
    sf=float(SizeFactor.loc[name])
    with open(i) as IN:
        for line in IN:
            a,b,c,d=line.rstrip('\n').split('\t')
            D=float(d)/sf
            nline='{}\t{}\t{}\t{:.6f}\n'.format(a,b,c,D)
            O.write(nline)
    O.close()
    Si=i[:-8]+'norm.sorted.bedGraph'
    scommend='sort -k1,1 -k2,2n {} > {}'.format(Norm,Si)
    os.system(scommend)
    os.system('rm {}'.format(Norm))
    BW=Si[:-8]+'bw'
    commend='bedGraphToBigWig {} /Path/to/ATAC-pipe-master/Data/Ref/hg19/hg19.sizes {}'.format(Si,BW)
    P=subprocess.Popen(commend,shell=True)
    P.wait()
    if P.returncode==0:
        print 'Normalize '+i+' Done!'

# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 10:37:27 2021

@author: Complxmat
"""

import os
import pandas as pd

use_topn=20
os.chdir('E:/0000 空间转录组/20211122 SPT Analysis-PTC/Step02--ContourAnalysis')
input_dir='Step04_GeneExpContourMatchCelltypeRegion'
sample_celltype_dirs=os.listdir(input_dir)
output_dir='Step04_GeneExpTopMatch'
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
for scd in sample_celltype_dirs:
    input_dir2=os.path.join(input_dir,scd)
    fns=os.listdir(input_dir2)
    f_out=open(os.path.join(output_dir,'Step04_%s_top%d_match.txt'%(scd,use_topn)),'w')
    f_out.write('\t'.join(['CelltypeRegion','TopGenes'])+'\n')
    for fn in fns:#ATPICAL CELL_Region_0_X17-33.0_Y27-62.0.txt
        df=pd.read_table(os.path.join(input_dir2,fn))
        df=df.sort_values(by=['ContourLevel','MatchQuality'],ascending=False)
        topn_genes=list(df['Feature'][0:use_topn-1])
        region_name=fn.replace('.0','').split('.')[0]
        f_out.write(region_name+'\t'+','.join(topn_genes)+'\n')
    f_out.close()
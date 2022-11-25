# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def read_file(fp):
    f=open(fp)
    ls=f.readlines()
    f.close()
    return ls


#Main
celltype2marker={'immune':'s',#AD
                 'tumor':'*',#AB
                 'ATPICAL CELL':'D',#A
                 'FOLLICULAR CELL':'.',#ABCD
                 'blood':'s',#BD
                 'FIBRO':'D',#BCD
                 'AFC':'v',
                 'AFCI':'^',
                 'FC1':'<',
                 'FC1I':'>',
                 'FC2':'1',
                 'FC2I':'2',
                 'I':'s',
                 'nTnFCnI':'x',
                 'T':'*'
                 }
ct_header_dict={'SubjType':'SubjType','TnT':'TnT20211031'}
xy_length_dict={'A':[80,120],
                'B':[80,130],
                'C':[70,100],
                'D':[80,120]}
# x_length=80
# y_length=120
os.chdir('E:/A--Lianlabwork/A018--YkSpt20211206/Step02--ContourAnalysis')
coord_dir='Step01_CoordinateAndCelltype'
coord_fns=os.listdir(coord_dir)
exp_dir='Step01_GeneScaleData'
output_dir0='Step02_GeneExpContour'
if not os.path.exists(output_dir0):
    os.mkdir(output_dir0)
for fn in coord_fns:
    sample_name,ct_info=fn.split('_')
    x_length,y_length=xy_length_dict[sample_name]
    ct_type=ct_info.split('.')[0].split('And')[1]
    exp_fn=sample_name+'_GeneExpScaled.txt'
    output_dir=os.path.join(output_dir0,sample_name+'_'+ct_type)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    #Read gene expression data
    #The data is scaled, negative values exist! 
    exp_table=pd.read_table(os.path.join(exp_dir,exp_fn))
    #20211122, use min exp - std exp as exp inf for blank spots
    value4blank_spots=[]
    for gene in exp_table.index:
        blank_value=min(exp_table.loc[gene])-np.std(exp_table.loc[gene])
        value4blank_spots.append(blank_value)
    blank_exp=pd.DataFrame(value4blank_spots)
    #Read spots coordinates and cell type
    spot_info_table=pd.read_table(os.path.join(coord_dir,fn))
    spots=spot_info_table.index
    coord_str2spot={}#{'59-15':AATCGAAA...-1,}
    spot2coords={}#{'AATCGAAA...-1:[59,15],}
    real_data_xs={}
    real_data_ys={}
    for spot in spots:
        spot_row=spot_info_table.loc[spot,'row']
        spot_col=spot_info_table.loc[spot,'col']
        coord_str2spot[str(spot_row)+'-'+str(spot_col)]=spot
        spot2coords[spot]=[spot_row,spot_col]
        spot_type=spot_info_table.loc[spot,ct_header_dict[ct_type]]
        if not spot_type in real_data_xs:
            real_data_xs[spot_type]=[spot_row]
            real_data_ys[spot_type]=[spot_col]
        else:
            real_data_xs[spot_type].append(spot_row)
            real_data_ys[spot_type].append(spot_col)
    
    #Padded spot has expression of the mean among 4 adjacent spots.
    #Spots without 4 real data adjacent spots are blank.
    #Padded exp table columns(spots) are ordered by xrange*yrange
    padded_exp_table=pd.DataFrame(index=exp_table.index)
    for x in range(0,x_length):
        for y in range(0,y_length):
            xystr=str(x)+'-'+str(y)
            if xystr in coord_str2spot:
                padded_exp_table[xystr]=exp_table[coord_str2spot[xystr]]
                #print(xystr,padded_exp_table[xystr])
            else:
                adj1str=str(x-1)+'-'+str(y)
                adj2str=str(x)+'-'+str(y-1)
                adj3str=str(x+1)+'-'+str(y)
                adj4str=str(x)+'-'+str(y+1)
                if adj1str in coord_str2spot and adj2str in coord_str2spot and adj3str in coord_str2spot and adj4str in coord_str2spot:
                    padded_exp=(exp_table[coord_str2spot[adj1str]]+exp_table[coord_str2spot[adj2str]]+exp_table[coord_str2spot[adj3str]]+exp_table[coord_str2spot[adj4str]])/4
                    padded_exp_table[xystr]=pd.DataFrame(padded_exp)
                    #print(xystr,padded_exp_table[xystr])
                else:
                    padded_exp_table[xystr]=blank_exp

    #Make contour plots
    selected_features=exp_table.index
    #selected_features=['TMSB4X']
    [X,Y]=np.meshgrid(range(0,x_length),range(0,y_length))
    for sf in selected_features:
        Z=np.array(padded_exp_table.loc[sf,])
        Z.shape=(x_length,y_length)
        plt.figure(figsize=(12,12))
        plt.contour(X,Y,Z.T,cmap='jet')
        for ct in real_data_xs:
            if ct=='NotDefined':
                plt.scatter(real_data_xs['NotDefined'],real_data_ys['NotDefined'],alpha=0.2,marker='x')
            else:
                plt.scatter(real_data_xs[ct],real_data_ys[ct],alpha=0.4,marker=celltype2marker[ct])
        plt.savefig(os.path.join(output_dir,sf+'.png'))
        plt.close()
    
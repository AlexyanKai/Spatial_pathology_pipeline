# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.path as mplp
from shapely.geometry import Polygon
import numpy as np

def read_file(fp):
    f=open(fp)
    ls=f.readlines()
    f.close()
    return ls


#Main
ct_header_dict={'SubjType':'SubjType','TnT':'TnT20211031'}
xy_length_dict={'A':[80,120],
                'B':[80,130],
                'C':[70,100],
                'D':[80,120]}
os.chdir('E:/A--Lianlabwork/A018--YkSpt20211206/Step02--ContourAnalysis')
coord_dir='Step01_CoordinateAndCelltype'
coord_fns=os.listdir(coord_dir)
exp_dir='Step01_GeneScaleData'
output_dir0='Step03_GeneExpContourCoverCelltype'
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
    subjtype2coords={}#{'FC1':[[x1,y1],[x2,y2],...],...}
    for spot in spots:
        spot_row=spot_info_table.loc[spot,'row']
        spot_col=spot_info_table.loc[spot,'col']
        coord_str2spot[str(spot_row)+'-'+str(spot_col)]=spot
        spot2coords[spot]=[spot_row,spot_col]
        spot_subjtype=spot_info_table.loc[spot,ct_header_dict[ct_type]]
        if not spot_subjtype in subjtype2coords:
            subjtype2coords[spot_subjtype]=[[spot_row,spot_col]]
        else:
            subjtype2coords[spot_subjtype].append([spot_row,spot_col])
    
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

    #For each celltype, check each feature's contour covering spots
    selected_features=exp_table.index
    #selected_features=['TMSB4X']
    [X,Y]=np.meshgrid(range(0,x_length),range(0,y_length))
    for sjt in subjtype2coords:
        if sjt=='NotDefined':
            continue
        outfp=os.path.join(output_dir,sjt+'_BeingCoveredStat.txt')
        if os.path.exists(outfp):
            print(sjt,' is already calculated and is skipped')
            continue
        f_out=open(outfp,'w')
        f_out.write('\t'.join(['BioFeature','ContourLevel','PolygonSeqNo',
                               'RowRange','ColRange','Area',
                               'CoveredSpots','Purity'])+'\n')
        spot_list=subjtype2coords[sjt]
        for sf in selected_features:
            Z=np.array(padded_exp_table.loc[sf,])
            Z.shape=(x_length,y_length)
            contour_obj=plt.contour(X,Y,Z.T)
            contour_levels=contour_obj.levels
            contour_level_segs=contour_obj.allsegs
            for level_seqNo in range(0,len(contour_levels)):
                if not contour_level_segs[level_seqNo]:
                    continue
                used_polygon_seqNo=0
                for point_set in contour_level_segs[level_seqNo]:
                    if len(point_set)<4:
                        continue
                    polygon_area=Polygon(point_set).area
                    if polygon_area<4:
                        continue
                    polygon_obj=mplp.Path(point_set)
                    number_of_covered_spots=sum(polygon_obj.contains_points(spot_list))
                    if number_of_covered_spots>0:
                        used_polygon_seqNo+=1
                        x_inf=int(min([p[0] for p in point_set]))
                        x_sup=np.ceil(max([p[0] for p in point_set]))
                        y_inf=int(min([p[1] for p in point_set]))
                        y_sup=np.ceil(max([p[1] for p in point_set]))
                        row_range=str(x_inf)+'-'+str(x_sup)
                        col_range=str(y_inf)+'-'+str(y_sup)
                        purity=number_of_covered_spots*100.0/polygon_area
                        l_out='\t'.join([sf,'Lv'+str(level_seqNo)+'/'+str(len(contour_levels))+':'+str(contour_levels[level_seqNo]),
                                        str(used_polygon_seqNo),row_range,col_range,
                                        str(np.ceil(polygon_area)),
                                        str(number_of_covered_spots),str(purity)])+'\n'
                        f_out.write(l_out)
        f_out.close()
        print(sjt,' cover stat finished')

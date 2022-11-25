# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 17:12:08 2021

@author: compl
"""

import os
import numpy as np
import pandas as pd
import networkx as nx
from shapely.geometry import Polygon
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt


#Main
#Set spatial transcriptome coordinate ranges
xy_length_dict={'A':[80,120],
                'B':[80,130],
                'C':[70,100],
                'D':[80,120]}
#Set parameters for region detection
link_distance=2#Use Chebyshev distance=max(xdiff,ydiff), not Euclidean distance
min_spot_count=5#Only region with equal or more spots are considered
ct_header_dict={'SubjType':'SubjType','TnT':'TnT20211031'}
os.chdir('E:/A--Lianlabwork/A018--YkSpt20211206/Step02--ContourAnalysis')
coord_dir='Step01_CoordinateAndCelltype'
coord_fns=os.listdir(coord_dir)
exp_dir='Step01_GeneScaleData'
output_dir0='Step04_GeneExpContourMatchCelltypeRegion'
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
    for target_celltype in subjtype2coords:
        if target_celltype=='NotDefined':
            continue
        #Use networkx to make a graph with one type of spots
        G=nx.Graph()
        nodes=[coord_str2spot[str(i[0])+'-'+str(i[1])] for i in subjtype2coords[target_celltype]]
        G.add_nodes_from(nodes)
        for n1 in range(0,len(nodes)-1):
            node1=nodes[n1]
            coord1=spot2coords[node1]
            for n2 in range(n1+1,len(nodes)):
                node2=nodes[n2]
                coord2=spot2coords[node2]
                if max(abs(coord1[0]-coord2[0]),abs(coord1[1]-coord2[1]))<=link_distance:
                    G.add_edge(node1,node2)
        
        #Extract connected components from G
        cc_list=sorted(nx.connected_components(G),key=len,reverse=True)
        print('There are %d regions(CC distance=%d) for %d %s spots'%(len(cc_list),link_distance,len(nodes),target_celltype))
        print('Region spot numbers:',[len(cc) for cc in cc_list])
        
        #Find ConvexHull for each CC as the region
        #Use vertices of the hull as region boundary
        region_boundary_list=[]
        for cc in cc_list:
            if len(cc)<min_spot_count:
                continue
            spot_coord_array=np.array([spot2coords[n] for n in cc])
            hull=ConvexHull(spot_coord_array)
            boundary_indice=hull.vertices
            boundary_coord_array=np.array([spot_coord_array[bi] for bi in boundary_indice])
            region_boundary_list.append(boundary_coord_array)
        print('Will analyze these regions:',region_boundary_list)
        
        #Plot contours and find the best level-polygon match with region
        #Best match should have high overlap quality=pct1*pct2
        #pct1=overlap_area/contour_area
        #pct2=overlap_area/selected_region_area
        [X,Y]=np.meshgrid(range(0,x_length),range(0,y_length))
        selected_features=exp_table.index
        region_best_match={}
        for region_seqNo in range(0,len(region_boundary_list)):
            region_best_match[region_seqNo]={}
        
        for sf in selected_features:
            Z=np.array(padded_exp_table.loc[sf,])
            Z.shape=(x_length,y_length)
            contour_obj=plt.contour(X,Y,Z.T)
            contour_levels=contour_obj.levels
            contour_level_segs=contour_obj.allsegs
            plt.close()
            for r_seqNo in range(0,len(region_boundary_list)):
                rb=region_boundary_list[r_seqNo]
                region_polygon_convex_hull=Polygon(rb).convex_hull
                overlap_contour_indice=[]
                overlap_quality_list=[]
                for level_seqNo in range(0,len(contour_levels)):
                    if not contour_level_segs[level_seqNo]:
                        continue
                    for point_set_seqNo in range(0,len(contour_level_segs[level_seqNo])):
                        if len(contour_level_segs[level_seqNo][point_set_seqNo])<3:
                            continue
                        contour_polygon_convex_hull=Polygon(contour_level_segs[level_seqNo][point_set_seqNo]).convex_hull
                        if contour_polygon_convex_hull.intersects(region_polygon_convex_hull):
                            oq=contour_polygon_convex_hull.intersection(region_polygon_convex_hull).area**2/contour_polygon_convex_hull.area/region_polygon_convex_hull.area*100.0
                            overlap_contour_indice.append([level_seqNo,point_set_seqNo])
                            overlap_quality_list.append(oq)
                if not overlap_quality_list:
                    continue
                highest_oq=max(overlap_quality_list)
                best_contour_indice=overlap_contour_indice[overlap_quality_list.index(highest_oq)]
                cx_inf=str(int(min([p[0] for p in contour_level_segs[best_contour_indice[0]][best_contour_indice[1]]])))
                cx_sup=str(np.ceil(max([p[0] for p in contour_level_segs[best_contour_indice[0]][best_contour_indice[1]]])))
                cy_inf=str(int(min([p[1] for p in contour_level_segs[best_contour_indice[0]][best_contour_indice[1]]])))
                cy_sup=str(np.ceil(max([p[1] for p in contour_level_segs[best_contour_indice[0]][best_contour_indice[1]]])))
                cxrange=str(cx_inf)+'-'+str(cx_sup)
                cyrange=str(cy_inf)+'-'+str(cy_sup)
                region_best_match[r_seqNo][sf]=[str(best_contour_indice[0])+'-out_of-'+str(len(contour_levels)),str(int(highest_oq)),cxrange,cyrange]
        
        #Write region best match into files
        for region_seqNo in region_best_match:
            rx_inf=str(int(min([p[0] for p in region_boundary_list[region_seqNo]])))
            rx_sup=str(np.ceil(max([p[0] for p in region_boundary_list[region_seqNo]])))
            ry_inf=str(int(min([p[1] for p in region_boundary_list[region_seqNo]])))
            ry_sup=str(np.ceil(max([p[1] for p in region_boundary_list[region_seqNo]])))
            rxrange=str(rx_inf)+'-'+str(rx_sup)
            ryrange=str(ry_inf)+'-'+str(ry_sup)
            f_out=open(os.path.join(output_dir,target_celltype+'_Region_'+str(region_seqNo)+'_X'+rxrange+'_Y'+ryrange+'.txt'),'w')
            f_out.write('\t'.join(['Feature','ContourLevel','MatchQuality','XRange','YRange'])+'\n')
            for feature in region_best_match[region_seqNo]:
                f_out.write(feature+'\t'+'\t'.join(region_best_match[region_seqNo][feature])+'\n')
            f_out.close()
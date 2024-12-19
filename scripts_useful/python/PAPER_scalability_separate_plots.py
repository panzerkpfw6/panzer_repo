#!/usr/bin/env python
import os,sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import math
#######################################
results={
    "shaheen": { "file":"shaheen_1st.txt"},
    "kanary": { "file":"kanary_1st.txt"},
}

#### data export from files
list_of_folders=['./scaling/7773X',
                 './scaling/7763',
                 './scaling/Flamingo',
                 './scaling/Kanary',
                 './scaling/Shaheen'
                 ]
times=[]
names=[]
cases=[]
giga_points=[]
giga_flops = []
point_updates = []
flops=[]
grids=[]
method=[]
code_type=[]
num_th=[]
sys=[]
for data_folder in list_of_folders:
    print(os.getcwd())
    path=os.path.join(data_folder,'*.log')
    FilenamesList = glob.glob(path)
    filenames_list=[]
    #### filter filenames list
    for word in FilenamesList:
        match_word = re.search(r'_.log', word)
        if match_word:
            continue
        else:
            filenames_list.append(word)
    filter_flag=False
    if filter_flag==True:
        wildcard_pattern = r'-TB_1st_'
        for word in FilenamesList:
            match_word = re.search(wildcard_pattern, word)
            if match_word:
                filenames_list.append(word)
                # print(word)# Initialiser les listes pour stocker les données
    ####
    for filename in filenames_list:
        with open(filename,'r') as file:
            tmp_time=math.nan
            tmp_giga_points=math.nan
            tmp_giga_flops=math.nan
            tmp_point_updates=math.nan
            tmp_flops=math.nan
            for line in file:
                if 'ELAPSED TIME (s)' in line and np.isnan(tmp_time):
                    tmp_time=(float(line.split()[-1]))
                elif 'GIGA POINT / s' in line and np.isnan(tmp_giga_points):
                    tmp_giga_points=(float(line.split()[-1]))
                elif 'GIGA FLOP / s' in line and np.isnan(tmp_giga_flops):
                    tmp_giga_flops=(float(line.split()[-1]))
                elif '# POINT UPDATES' in line and np.isnan(tmp_point_updates):
                    tmp_point_updates=(float(line.split()[-1]))
                elif 'FLOP' in line and np.isnan(tmp_flops):
                    tmp_flops=(float(line.split()[-1]))   
        if np.isnan(tmp_time):
            continue
        else:
            sys.append(data_folder.split('/')[2])
            times.append(tmp_time)
            giga_points.append(tmp_giga_points)
            giga_flops.append(tmp_giga_flops)
            point_updates.append(tmp_point_updates)
            flops.append(tmp_flops)
            tmp=filename.split('/')[-1]
            tmp=tmp.split('.log')[0]
            method.append( tmp.split('_')[0].split('-')[1] )
            code_type.append(tmp.split('_')[1])
            num_th.append( float (tmp.split('_')[2]) )
            # grids.append( tmp.split('_')[8]+'_'+tmp.split('_')[9]+'_'+tmp.split('_')[10] )
            names.append(filename)

# Créer un DataFrame à partir des listes
data = pd.DataFrame({
    'sys':sys,
    'names':names,
    'num_th':num_th,
    'giga_point_s':giga_points,
    'giga_flop_s': giga_flops,
    'point_updates': point_updates,
    'Flops': flops,
    'times':times,
    'method':method,
    'code_type':code_type
    })
data.columns
#######################
data2=data.copy()
data2.columns
data2=data2.drop(['names','point_updates','Flops','times'], axis=1)
####################### Plotting
figures_root='../scaling'
grid_name='2048 x 2048 x 512'
colors=['blue','orange']
####################### 1st order
labels={
        'Shaheen-SB':'SB, 1st order, Haswell',
        'Shaheen-TB':'MWD, 1st order, Haswell',
        'Kanary-SB':'SB, 1st order, AMD EPYC 7713',
        'Kanary-TB':'MWD, 1st order, AMD EPYC 7713',
        'Flamingo-SB':'SB, 1st order, Intel Icelake',
        'Flamingo-TB':'MWD, 1st order, Intel Icelake',
        '7773X-SB':'SB, 1st order, MilanX',
        '7773X-TB':'MWD, 1st order, MilanX',
        '7763-SB':'SB, 1st order, AMD EPYC 7763',
        '7763-TB':'MWD, 1st order, AMD EPYC 7763',
}
code_type_='1st-abc'
code_type_='1st'
data3 = data2[data2['code_type']==code_type_]    #1st-abc,1st
dataSB = data3[data3['method']=='SB']
dataTB = data3[data3['method']=='TB']
print(data3['sys'].unique())
systems=data3['sys'].unique()
# for sys in systems[0:1]:
for sys in systems:
    metric_name='giga_point_s'
    data_sys_sb = dataSB[dataSB['sys']==sys]
    data_sys_sb=data_sys_sb.sort_values(by=[metric_name],ascending=False)
    data_sys_tb = dataTB[dataTB['sys']==sys]
    data_sys_tb=data_sys_tb.sort_values(by=[metric_name],ascending=False)

    fig=plt.figure(figsize=(7,2.0))
    
    for icolor,code in enumerate(['TB', 'SB']):
        if code=='SB':   yval=data_sys_sb[metric_name];sb=yval.array[0];
        if code=='TB':   yval=data_sys_tb[metric_name];tb=yval.array[0];
        plt.plot(data_sys_sb['num_th'],yval,
                label=labels[sys+'-'+code],marker='o', linestyle='-',
                markersize=6,mfc='white', mew=2, lw=2,
                zorder=-1, clip_on=False,color=colors[icolor])
    ax=plt.gca()
    # plt.plot(data_sys_sb['num_th'],data_sys_tb[metric_name])
    plt.axis('tight')
    ax.grid(True);
    ax.set_xlabel('Number of threads')
    ax.set_ylabel('GStencils/sec',fontweight='bold')
    ax.set_title('Grid size:'+grid_name)
    plt.legend(handlelength=3,loc='upper left')
    # algo speedup
    line1=ax.arrow(data_sys_sb['num_th'].array[0],sb,0,(tb-sb), color='g', width=0.40, length_includes_head=True, head_length=1)
    ax.annotate(str(np.round(tb/sb,1))+"X", xy=(data_sys_sb['num_th'].array[0], (sb+tb)/2), xytext=(data_sys_sb['num_th'].array[0]-data_sys_sb['num_th'].array[0]/10,(sb+tb)/2-0.9), color='g', fontsize=14)
    plt.show()
    fig.savefig(os.path.join(figures_root,sys+'_'+code_type_+'_2k_grid.png'), transparent=False, bbox_inches="tight",dpi=350)
    fig.savefig(os.path.join(figures_root,sys+'_'+code_type_+'_2k_grid.pdf'), transparent=False, bbox_inches="tight")
    plt.close()
    ss=1
dd=1

####################### 2nd order
# code_type_='2nd-abc'
code_type_='2nd'
labels={
        'Shaheen-SB':'SB, 2nd order, Haswell',
        'Shaheen-TB':'MWD, 2nd order, Haswell',
        'Kanary-SB':'SB, 2nd order, AMD EPYC 7713',
        'Kanary-TB':'MWD, 2nd order, AMD EPYC 7713',
        'Flamingo-SB':'SB, 2nd order, Intel Icelake',
        'Flamingo-TB':'MWD, 2nd order, Intel Icelake',
        '7773X-SB':'SB, 2nd order, MilanX',
        '7773X-TB':'MWD, 2nd order, MilanX',
        '7763-SB':'SB, 2nd order, AMD EPYC 7763',
        '7763-TB':'MWD, 2nd order, AMD EPYC 7763',
}
data3 = data2[data2['code_type']==code_type_]    #1st-abc,1st
dataSB = data3[data3['method']=='SB']
dataTB = data3[data3['method']=='TB']
print(data3['sys'].unique())
systems=data3['sys'].unique()
# for sys in systems[0:1]:
for sys in systems:
    metric_name='giga_point_s'
    data_sys_sb = dataSB[dataSB['sys']==sys]
    data_sys_sb=data_sys_sb.sort_values(by=[metric_name],ascending=False)
    data_sys_tb = dataTB[dataTB['sys']==sys]
    data_sys_tb=data_sys_tb.sort_values(by=[metric_name],ascending=False)

    fig=plt.figure(figsize=(7,2.0))
    for icolor,code in enumerate(['TB', 'SB']):
        if code=='SB':   yval=data_sys_sb[metric_name];sb=yval.array[0];
        if code=='TB':   yval=data_sys_tb[metric_name];tb=yval.array[0];
        plt.plot(data_sys_sb['num_th'],yval,
                label=labels[sys+'-'+code],marker='o', linestyle='-',
                markersize=6,mfc='white', mew=2, lw=2,
                zorder=-1, clip_on=False,color=colors[icolor])
    ax=plt.gca()
    plt.axis('tight')
    ax.grid(True);
    ax.set_xlabel('Number of threads')
    ax.set_ylabel('GStencils/sec',fontweight='bold')
    if grid_name=='2048 x 2048 x 512':      grid_name_str='2048*2048*512';
    ax.set_title('Grid size:'+grid_name)
    plt.legend(handlelength=3,loc='upper left')
    # algo speedup
    line1=ax.arrow(data_sys_sb['num_th'].array[0],sb,0,(tb-sb), color='g', width=0.40, length_includes_head=True, head_length=1)
    ax.annotate(str(np.round(tb/sb,1))+"X", xy=(data_sys_sb['num_th'].array[0], (sb+tb)/2), xytext=(data_sys_sb['num_th'].array[0]-data_sys_sb['num_th'].array[0]/10,(sb+tb)/2-0.9), color='g', fontsize=14)
    plt.show()
    fig.savefig(os.path.join(figures_root,sys+'_'+code_type_+'_2k_grid.png'), transparent=False, bbox_inches="tight",dpi=350)
    fig.savefig(os.path.join(figures_root,sys+'_'+code_type_+'_2k_grid.pdf'), transparent=False, bbox_inches="tight")
    plt.close()
ss=1


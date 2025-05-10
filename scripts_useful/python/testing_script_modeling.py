#!/usr/bin/env python
import os,sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import re
import math
import glob
import warnings
# Suppress all warnings
warnings.filterwarnings("ignore")
sys.stdout = open(os.devnull,'w')

folder_genoa=os.path.join('./data/par_search/')

def read_data_stencil_rtm(filename):
    ###
    # read SB data from stencil-main logging files 
    ### 
    times=[]
    giga_points=[]
    giga_flops=[]
    point_updates=[]
    flops=[]
    method=[]
    grids=[]
    cores=[]
    cb_size=[]
    cb_x=[]
    cb_y=[]
    cb_z=[]
    thx=[]
    thy=[]
    thz=[]
    tdim=[]
    numwf=[]
    ##################
    tmp_ncores=math.nan
    tmp_time=math.nan
    tmp_giga_points=math.nan
    tmp_giga_flops=math.nan
    tmp_point_updates=math.nan
    tmp_flops=math.nan
    tmp_grid=math.nan
    tmp_cb_x=math.nan
    tmp_cb_y=math.nan
    tmp_cb_z=math.nan
    tmp_thx=math.nan
    tmp_thy=math.nan
    tmp_thz=math.nan
    tmp_tdim=math.nan
    tmp_numwf=math.nan
    tmp_mode=math.nan
    tmp_diam_width=math.nan
    tmp_method=math.nan
    tmp_cb_size=math.nan
    counter=0
    ##################
    ###### Parse files and aggregate them into the dataframe #######################
    print(os.getcwd())
    if not os.path.exists(filename):
        print(f"Warning: File '{filename}' not found. Skipping.")
        return None
    with open(os.path.join(filename),'r') as file:
        previous_line=''
        for line in file:
            if ('Program started at' in line):
                tmp_time=math.nan
                tmp_giga_points=math.nan
                tmp_giga_flops=math.nan
                tmp_point_updates=math.nan
                tmp_flops=math.nan
                tmp_grid=math.nan
                tmp_cb_x=math.nan
                tmp_cb_y=math.nan
                tmp_cb_z=math.nan
                tmp_thx=math.nan
                tmp_thy=math.nan
                tmp_thz=math.nan
                tmp_tdim=math.nan
                tmp_numwf=math.nan
                tmp_mode=math.nan
            if 'velocity size' in line and issubclass(type(tmp_grid),str)==False:
                pattern = r"velocity size\s*=\s*(\d+)\s*x\s*(\d+)\s*x\s*(\d+)"
                match = re.search(pattern,line)
                grid_x = int(match.group(1))
                grid_y = int(match.group(2))
                grid_z = int(match.group(3))
                tmp_grid = match.group(1)+" x "+match.group(2)+" x "+match.group(3)
            elif '[STENCIL MSG]:Total:' in line and '[STENCIL MSG]:Global info:' in previous_line:
                ### account for running time in TB method
                pattern = r'\[STENCIL MSG\]:Total:\s*([\d.]+)\s*\(s\)'
                match = re.search(pattern, line)
                number = match.group(1)
                tmp_time=float(number)
            elif '[STENCIL MSG]:Total:' in line and '[STENCIL MSG]:forward timer' in previous_line:
                ### account for running time in SB method
                pattern = r'\[STENCIL MSG\]:Total:\s*([\d.]+)\s*\(s\)'
                match = re.search(pattern, line)
                number = match.group(1)
                tmp_time=float(number)
            elif '[STENCIL MSG]:Total:' in line and '[STENCIL MSG]:Speed info:' in previous_line:
                ### account for performance in TB method
                pattern = r'\[STENCIL MSG\]:Total:\s*([\d.]+)\s*GStencils/s'
                match = re.search(pattern, line)
                if match:   number = match.group(1)
                tmp_giga_points=float(number)
            elif '[STENCIL MSG]:PropSpeed:' in line and '[STENCIL MSG]:Speed:' in previous_line:
                ### account for performance in SB method
                pattern = r'\[STENCIL MSG\]:PropSpeed:\s*([\d.]+)\s*GStencils/s'
                match = re.search(pattern, line)
                if match:   number = match.group(1)
                tmp_giga_points=float(number)
            elif 'run 1st order' in line:
                if 'TB' in line:        tmp_method='tb_abc'
                elif 'SB' in line:      tmp_method='sb_abc'
                elif 'modeling' in line:      tmp_mode='modeling'
                elif 'RTM' in line:      tmp_mode='rtm'
            elif '# THREADS' in line:
                tmp_ncores=(float(line.split()[-1]))
            elif 'BLOCKX=' in line:
                pattern = r"BLOCKX=(\d+)\s*,\s*BLOCKY=(\d+)\s*,\s*BLOCKZ=(\d+)"
                match = re.search(pattern,line)
                if match:
                    tmp_cb_x = int(match.group(1))
                    tmp_cb_y = int(match.group(2))
                    tmp_cb_z = int(match.group(3))
                tmp_cb_size=[tmp_cb_x,tmp_cb_y,tmp_cb_z]
            elif '[STENCIL MSG]:t_dim' in line and '[STENCIL MSG]:temporal blocking' in previous_line:
                # Define temporal blocking parameters in TB method
                pattern = r"\[STENCIL MSG\]:t_dim : (\d+), num_wf : (\d+), diam_width : (\d+)"
                match = re.search(pattern, line)
                if match:
                    tmp_tdim = int(match.group(1))  # t_dim
                    tmp_numwf = int(match.group(2))  # t_dim
                    tmp_diam_width = int(match.group(3))  # t_dim
            elif '[STENCIL MSG]:thread group' in line:
                pattern = r"\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)"  # Matches (number,number,number)
                match = re.search(pattern, line)
                if match:
                    tmp_thx = int(match.group(1))  # Extract thx value
                    tmp_thy = int(match.group(2))  # Extract thy value
                    tmp_thz = int(match.group(3))  # Extract thz value
            elif ('[STENCIL MSG]:END of modeling' in line) or ('RTM HALAS' in line):
                cores.append(tmp_ncores)
                times.append(tmp_time)
                giga_points.append(tmp_giga_points)
                giga_flops.append(tmp_giga_flops)
                point_updates.append(tmp_point_updates)
                flops.append(tmp_flops)
                method.append(tmp_method)
                grids.append(tmp_grid)
                cb_size.append(tmp_cb_size)
                cb_x.append(tmp_cb_x)
                cb_y.append(tmp_cb_y)
                cb_z.append(tmp_cb_z)
                thx.append(tmp_thx)
                thy.append(tmp_thy)
                thz.append(tmp_thz)
                tdim.append(tmp_tdim)
                numwf.append(tmp_numwf)
                counter=counter+1
            previous_line=line
    ################################################
    # Créer un DataFrame à partir des listes
    data = pd.DataFrame({
        'method':method,
        'times':times,
        'giga_point_s':giga_points,
        'grids':grids,
        'cb_size':cb_size,
        'cb_x': cb_x,
        'cb_y': cb_y,
        'cb_z': cb_z,
        'cb_z': cb_z,
        'thx': thx,
        'thy': thy,
        'thz': thz,
        'tdim': tdim,
        'numwf': numwf
        })
    print(data.columns)
    return data

#####################################################################
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
tick_label_size=11
def plot_perf(data,save_paths,title='',metric='gstencils',gflops_limit=3300,gstencils_limit=100):
    barWidth = 0.25
    font = {'size': 17}
    plt.rc('font', **font)
    plt.rcParams["font.weight"]="bold"
    n_axis=3
    # fig, ax = plt.subplots(1,len(data_list),figsize=(7.5*len(data_list),5))
    # fig, ax = plt.subplots(1,len(data_list),figsize=(22*len(data_list),2.0))  # first variant
    fig, ax = plt.subplots(1,n_axis,figsize=(7.5*n_axis,2.0))  # second variant
    # plt.subplots_adjust(left=0.1, right=0.97,wspace=0.39, hspace=0.4)
    # plt.subplots_adjust(left=0.1, right=0.97,wspace=0.2, hspace=0.5)
    # Set position of bar on X axis 
    br1 = np.arange(1)
    br2 = [x + barWidth for x in br1]
    
    grids=data['grids'].unique()
    ##### for (counter,data) in enumerate(data):
    for i_axis in range(n_axis):
        # Make the plot
        if n_axis==1:
            AX=ax
        else:
            AX=ax[i_axis]
        print('data=',data)
        
        grid_name=grids[i_axis];
        data_grid=data[data['grids']==grid_name]
        data_grid[data_grid['method']=='sb_abc']
        SB=data_grid[data_grid['method']=='sb_abc'].iloc[0]['giga_point_s']
        TB=data_grid[data_grid['method']=='tb_abc'].iloc[0]['giga_point_s']

        bars=AX.bar(br1, SB, color ='b', width = barWidth,edgecolor ='grey', label ='SB')
        for bar in bars:
            yval = bar.get_height()
            if metric=='gflops':
                plot_val=int(yval)
            else:
                plot_val=round(yval,1)
            ### AX.text(bar.get_x()-0.04,yval+yval*0.02,plot_val,fontsize=14,fontweight='bold')   # bold, normal
            if np.isfinite(yval):
                AX.text(bar.get_x() - 0.04, yval + yval * 0.02, plot_val, fontsize=14, fontweight='bold')
            else:
                print(f"Warning: Invalid yval at bar {bar.get_x()}: {yval}")

        bars=AX.bar(br2, TB, color ='r', width = barWidth,edgecolor ='grey', label ='MWD')
        for bar in bars:
            yval = bar.get_height()
            if metric=='gflops':
                plot_val=int(yval)
            else:
                plot_val=round(yval,1)
            #### AX.text(bar.get_x()+0.04,yval+yval*0.02,plot_val,fontsize=14,fontweight='bold')
            if np.isfinite(yval):
                AX.text(bar.get_x() - 0.04, yval + yval * 0.02, plot_val, fontsize=14, fontweight='bold')
            else:
                print(f"Warning: Invalid yval at bar {bar.get_x()}: {yval}")

        # AX.set_xticks([r + barWidth for r in range(len(SB))],['Stencil', 'Stencil+ABCs'])
        AX.set_xticks([0.125,1.125],['Stencil+ABCs,\n 1st order','Stencil+ABCs,\n 2nd order'], fontweight ='bold', fontsize = 15)
        # AX.set_yticks(fontsize=15)
        # plt.yticks(fontsize=13)
        AX.tick_params(axis='y', labelsize=13)

        if metric=='gflops':
            AX.set_ylabel('GFlop/s', fontweight ='bold', fontsize = 15)
            AX.set_ylim(bottom=0,top=gflops_limit)
        else:
            # AX.set_ylabel('GStencils/s', fontweight ='bold', fontsize = 18)
            AX.set_ylabel('GStencils/s', fontweight ='bold', fontsize = 15)
            AX.set_ylim(bottom=0,top=gstencils_limit)
        # AX.legend(loc='upper right')
        # AX.legend(loc='upper left')
        AX.legend(loc='upper center')
        # AX.set_title('Grid size:\n'+grid_name,fontweight='bold',fontsize=14)
        AX.set_title(title+', grid size:'+grid_name,fontweight='bold',fontsize=17)
    my_suptitle=fig.suptitle('',fontsize=18,y=0,weight='bold')
    # my_suptitle=fig.suptitle(title,fontsize=18,y=0,weight='bold')
    plt.show()
    fig.savefig(save_paths[0]+'.png',bbox_inches='tight',bbox_extra_artists=[my_suptitle])
    fig.savefig(save_paths[0]+'.pdf',bbox_inches='tight',bbox_extra_artists=[my_suptitle])
    return None

#####################################################################
#####################################################################
save_path='/media/plotnips/sdd1/Dropbox/Apps/Overleaf/pavel_phd_thesis/fig/rtm/figures/performances/'

perf_data=read_data_stencil_rtm(os.path.join(folder_genoa,'test1_forward.log'))
save_paths=[os.path.join(save_path,'perf_fwd_genoa')]
plot_perf(perf_data,save_paths,metric='gstencils',gstencils_limit=150)
ss=1
#####################################################################

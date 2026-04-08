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


def read_sb_data_stencil_main(filename):
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
    ##################
    ###### Parse files and aggregate them into the dataframe #######################
    print(os.getcwd())
    if not os.path.exists(filename):
        print(f"Warning: File '{filename}' not found. Skipping.")
        return None
    with open(os.path.join(filename),'r') as file:
        previous_line=''
        for line in file:
            if ('START SB KERNEL' in line) or ('START TB KERNEL' in line):
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
            if 'GRID' in line and issubclass(type(tmp_grid), str)==False:
                    tmp_grid=(line.split('GRID')[-1])
                    tmp_grid=(tmp_grid.split('\n')[0])
                    tmp_grid=tmp_grid.strip()
            elif 'ELAPSED TIME (s)' in line and np.isnan(tmp_time):
                tmp_time=(float(line.split()[-1]))
            elif 'GIGA POINT / s' in line and np.isnan(tmp_giga_points):
                tmp_giga_points=(float(line.split()[-1]))
            elif 'GIGA FLOP / s' in line and np.isnan(tmp_giga_flops):
                tmp_giga_flops=(float(line.split()[-1]))
            elif '# POINT UPDATES' in line and np.isnan(tmp_point_updates):
                tmp_point_updates=(float(line.split()[-1]))
            elif 'FLOP' in line and np.isnan(tmp_flops):
                tmp_flops=(float(line.split()[-1]))
            elif 'R-' in line:
                line=line.split()[0]
                if 'SB' in line:    tmp_method='sb'
                else:   tmp_method='tb'
                if '2' in line: tmp_method=tmp_method+'_order2'
                if 'withABC' in line:   tmp_method=tmp_method+'_abc'
            elif '# THREADS' in line:
                tmp_ncores=(float(line.split()[-1]))
            elif 'CB SIZE' in line:
                pattern = r"\d+\s*/\s*\d+\s*/\s*\d+"
                match = re.search(pattern,line)
                values = match.group(0).split(' / ')
                tmp_cb_x=int(values[0])  # Extract BLOCKX value
                tmp_cb_y=int(values[1])
                tmp_cb_z=int(values[2])
                tmp_cb_size=[tmp_cb_x,tmp_cb_y,tmp_cb_z]
            elif 'thx,thy,thz' in line and 'Usage:' not in line:
                pattern = r"\d+,\d+,\d+"  # Pattern to match numbers separated by commas
                match = re.search(pattern,line)
                if match:
                    values = match.group(0).split(',')  # Split the matched string into individual values
                    tmp_thx = int(values[0])  # Extract thx value
                    tmp_thy = int(values[1])  # Extract thy value
                    tmp_thz = int(values[2])  # Extract thz value
                    tmp_th_values=[tmp_thx,tmp_thy,tmp_thz]
            elif 't_dim' in line and 'Usage:' not in line:
                tmp_tdim=(float(line.split()[-1]))
            elif 'num_wf' in line and 'Usage:' not in line:
                tmp_numwf=(float(line.split()[-1]))
            elif ('END SB KERNEL' in line) or ('END TB KERNEL' in line):
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
            previous_line=line
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

plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
def plot_heat_map(map,title,x_label,y_label,xticks,yticks,save_path):
    # Create the figure with the specified size
    fig, ax = plt.subplots(figsize=(8, 6))

    # Use make_axes_locatable to create a new axis for the colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)  # 'size' controls colorbar width, 'pad' controls spacing

    # Display the image
    im = ax.imshow(map.T, cmap='viridis', aspect='auto')

    # Add gridlines (if needed)
    ax.grid(True, color='black', linestyle='--', linewidth=0.5)

    # Create the colorbar and attach it to the new axis
    fig.colorbar(im, cax=cax, orientation='vertical')

    # Set plot title and labels (use ax for titles and labels)
    ax.set_title(title,fontsize=14)
    ax.set_xlabel(x_label,fontsize=14, fontweight='bold')  # Set your desired label for the X-axis
    ax.set_ylabel(y_label,fontsize=14, fontweight='bold')  # Set your desired label for the Y-axis

    # Set the x and y ticks with corresponding labels
    ax.set_xticks(np.arange(len(xticks)))
    ax.set_xticklabels(xticks, rotation=45, fontweight='bold')  # You can rotate the labels if needed
    ax.set_yticks(np.arange(len(yticks)))
    ax.set_yticklabels(yticks, fontweight='bold')

    print('save to ',save_path)
    # Save the figure with high resolution
    plt.savefig(save_path, dpi=400, bbox_inches='tight')

    # Show the plot
    plt.show()

    # Close the plot to free resources
    # plt.close()
    return None
tick_label_size=11

#####################################################################
perf_data=read_data_stencil_rtm(os.path.join(folder_genoa,'test1_forward.log'))
ss=1
#####################################################################

cbx_val=np.arange(2, 65, 2);cbx_val = np.insert(cbx_val, 0, 1)
cby_val=np.arange(2, 65, 2);cby_val = np.insert(cby_val, 0, 1)
nx=len(cbx_val)
ny=len(cby_val)
grids_=['512 x 512 x 512','1024 x 1024 x 512','2048 x 2048 x 512']

sb_data=read_data_stencil_rtm(os.path.join(folder_genoa,'test2_united.log'))   #sb_data
tb_data=read_data_stencil_rtm(os.path.join(folder_genoa,'test3_united.log'))
data=pd.concat([sb_data,tb_data],ignore_index=True)
data_summary=data.copy()
####################### Create summary for dataframe  #######################
data_summary = data[['method', 'grids', 'times', 'giga_point_s', 'cb_size','cb_x','cb_y','cb_z']]
grid_name='1024 x 1024 x 512'
perf_map_1st = np.ones((ny,nx))*(float('nan'))

####################### Find best parameters for SB 1st order #######################
data_summary_1st = data_summary[data_summary['method'] == 'sb_abc']
data_summary_1st = data_summary_1st[data_summary_1st['grids'] == grid_name]
print(data_summary['grids'].unique())
len(data_summary_1st)

for i in range(len(data_summary_1st)):  # Loops from 0 to 4
    if (data_summary_1st.iloc[i]['grids']==grid_name):
        cb_x_value = data_summary_1st.iloc[i]['cb_x']
        cb_y_value = data_summary_1st.iloc[i]['cb_y']
        x_index = np.where(cbx_val == cb_x_value)[0][0]  # [0][0] gets the first index
        y_index = np.where(cby_val == cb_y_value)[0][0]  # [0][0] gets the first index
        perf_map_1st[y_index,x_index] = data_summary_1st.iloc[i]['giga_point_s']
        if data_summary_1st.iloc[i]['giga_point_s']==-10000:
            aa=1
##################
map_=perf_map_1st
title_='AMD Genoa, SB 1st order, GStencils/sec, grid='+grid_name
xlabel_='Cache blocking in Y dimension, grid points'
ylabel_='Cache blocking in X dimension, grid points'
xticks=cbx_val
yticks=cby_val
save_path='genoa_sb_1st_1024.png'
plot_heat_map(map_,title_,xlabel_,ylabel_,xticks,yticks,save_path)
##########################################################################################
##########################################################################################
best_parameters_SB=pd.DataFrame()
for grid_name in grids_:
    ####################### Find best parameters for SB 1st order #######################
    data_summary_1st = data_summary[data_summary['method'] == 'sb_abc']
    data_summary_1st = data_summary_1st[data_summary_1st['grids'] == grid_name]
    data_summary_1st=data_summary_1st.sort_values(by=['giga_point_s'],ascending=False)
    best_parameters_SB=pd.concat([best_parameters_SB,data_summary_1st.iloc[0]],ignore_index=True)
##########################################################################################
best_parameters_TB=pd.DataFrame()
for grid_name in grids_[0:2]:
    ####################### Find best parameters for TB 1st order #######################
    data_summary_1st = data_summary[data_summary['method'] == 'tb_abc']
    data_summary_1st = data_summary_1st[data_summary_1st['grids'] == grid_name]
    data_summary_1st=data_summary_1st.sort_values(by=['giga_point_s'],ascending=False)
    best_parameters_TB=pd.concat([best_parameters_TB,data_summary_1st.iloc[0]],ignore_index=True)
print(best_parameters_TB)
##########################################################################################
ss=1
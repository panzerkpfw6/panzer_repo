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
sys.stdout = open(os.devnull, 'w')

####### check logs
# dropbox_root='/home/plotnips/Dropbox'
dropbox_root='/media/plotnips/sdd1/Dropbox'
# dropbox_root='/Users/pavelplotnitskii/Dropbox'
folder_milanX=os.path.join(dropbox_root,'Apps/download_to_laptop/rtm_tests_milanx/')

def read_sb_data(filename):
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
            if ('Program started at:' in line):
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
                tmp_ncores=math.nan
            elif ('run 1st order SB' in line):    tmp_method='sb_abc'
            elif ('run 2nd order SB' in line):    tmp_method='sb_order2_abc'
            elif ('run 1st order TB' in line):    tmp_method='tb_abc'
            elif ('run 2nd order TB' in line):    tmp_method='tb_order2_abc'
            elif 'velocity size' in line and issubclass(type(tmp_grid),str)==False:
                pattern = r"velocity size\s*=\s*(\d+)\s*x\s*(\d+)\s*x\s*(\d+)"
                match = re.search(pattern,line)
                grid_x = int(match.group(1))
                grid_y = int(match.group(2))
                grid_z = int(match.group(3))
                tmp_grid = match.group(1)+" x "+match.group(2)+" x "+match.group(3)
            # sb time logging
            elif ('[STENCIL MSG]:Total:' in line) and ('sb' in tmp_method): # sb time
                pattern = r'\[STENCIL MSG\]:Total:\s*([\d.]+)\s*\(s\)'
                match = re.search(pattern,line)
                number = match.group(1)
                tmp_time=float(number)
            elif ('[STENCIL MSG]:PropSpeed:' in line) and ('sb' in tmp_method): # sb time
                pattern = r'\[STENCIL MSG\]:PropSpeed:\s*([\d.]+)\s*GStencils/s'
                match = re.search(pattern,line)
                if match:
                    number = match.group(1)
                tmp_giga_points=float(number)
            # tb time logging
            elif ('[STENCIL MSG]:Total:' in line) and ('[STENCIL MSG]:Global info:' in previous_line) and ('tb' in tmp_method):  # tb time
                pattern = r'\[STENCIL MSG\]:Total:\s*([\d.]+)\s*\(s\)'
                match = re.search(pattern,line)
                number = match.group(1)
                tmp_time=float(number)
            elif ('[STENCIL MSG]:Total:' in line) and ('[STENCIL MSG]:Speed info:' in previous_line) and ('tb' in tmp_method):
                pattern = r'\[STENCIL MSG\]:Total:\s*([\d.]+)\s*GStencils/s'
                match = re.search(pattern,line)
                if match:
                    number = match.group(1)
                tmp_giga_points=float(number)
            elif '# THREADS' in line:
                tmp_ncores=(float(line.split()[-1]))
            elif '[STENCIL MSG]:BLOCKX' in line:
                # Pattern to match BLOCKX, BLOCKY, BLOCKZ values (e.g., BLOCKX=1, BLOCKY=4, BLOCKZ=1)
                pattern = r"BLOCKX=(\d+),\s*BLOCKY=(\d+),\s*BLOCKZ=(\d+)"
                match = re.search(pattern, line)
                if match:
                    tmp_cb_x = int(match.group(1))  # Extract BLOCKX value
                    tmp_cb_y = int(match.group(2))  # Extract BLOCKY value
                    tmp_cb_z = int(match.group(3))  # Extract BLOCKZ value
                    tmp_cb_size = [tmp_cb_x, tmp_cb_y, tmp_cb_z]
            elif '[STENCIL MSG]:thread group' in line:
                # Pattern to match thread group values in the format (x,y,z)
                pattern = r"\((\d+),(\d+),(\d+)\)"
                match = re.search(pattern,line)
                if match:
                    tmp_thx = int(match.group(1))  # Extract first value (e.g., 8)
                    tmp_thy = int(match.group(2))  # Extract second value (e.g., 1)
                    tmp_thz = int(match.group(3))  # Extract third value (e.g., 1)
                    tmp_th_values = [tmp_thx,tmp_thy,tmp_thz]
            elif ('t_dim' in line) and ('[STENCIL MSG]:temporal blocking' in previous_line):
                pattern = r"t_dim\s*:\s*(\d+),\s*num_wf\s*:\s*(\d+),\s*diam_width\s*:\s*(\d+)"
                match = re.search(pattern, line)
                if match:
                    tmp_tdim = int(match.group(1))      # Extract t_dim value (e.g., 7)
                    tmp_numwf = int(match.group(2))     # Extract num_wf value (e.g., 32)
            elif '[STENCIL MSG]:END' in line:
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

    # Save the figure with high resolution
    plt.savefig(save_path, dpi=400, bbox_inches='tight')

    # Show the plot
    plt.show()

    # Close the plot to free resources
    plt.close()
    return None
tick_label_size=11

grids_=['512 x 512 x 512','1024 x 1024 x 512','2048 x 2048 x 512']

data=read_sb_data(os.path.join(folder_milanX,'test2.log'))   #sb_data
# tb_data=read_sb_data(os.path.join(folder_milanX,'test3.log'))
# data=pd.concat([sb_data,tb_data],ignore_index=True)
print(data['grids'].unique())
data_milanX=data

####################### Create summary for dataframe  #######################
data_summary = data[['method', 'grids', 'times', 'giga_point_s', 'cb_size','cb_x','cb_y','cb_z']]
grid_name='512 x 512 x 512'

####################### Find best parameters for SB 1st order #######################
data_summary_1st = data_summary[data_summary['method'] == 'sb_abc']
data_summary_1st = data_summary_1st[data_summary_1st['grids'] == grid_name]
print(data_summary_1st['grids'].unique())
print(data_summary_1st['cb_y'].unique())
data_summary_1st_sorted = data_summary_1st.sort_values(by='giga_point_s',ascending=False)


cbx_val=data_summary_1st['cb_x'].unique()
cby_val=data_summary_1st['cb_y'].unique()
cbz_val=data_summary_1st['cb_z'].unique()
nx=len(cbx_val);ny=len(cby_val);nz=len(cbz_val);
perf_map_1st = np.ones((ny,nx,nz))*(float('nan'))

for i in range(len(data_summary_1st)):
    cb_x_value = data_summary_1st.iloc[i]['cb_x']
    cb_y_value = data_summary_1st.iloc[i]['cb_y']
    cb_z_value = data_summary_1st.iloc[i]['cb_z']
    x_index = np.where(cbx_val == cb_x_value)[0][0]
    y_index = np.where(cby_val == cb_y_value)[0][0]
    z_index = np.where(cbz_val == cb_z_value)[0][0]
    perf_map_1st[y_index,x_index,z_index] = data_summary_1st.iloc[i]['giga_point_s']
    if data_summary_1st.iloc[i]['giga_point_s']==-10000:
        aa=1

##################
map_=perf_map_1st[:,:,4]
title_='milanX, SB 1st order, GStencils/sec, grid='+grid_name
xlabel_='Cache blocking in Y dimension, grid points'
ylabel_='Cache blocking in X dimension, grid points'
xticks=cbx_val
yticks=cby_val
save_path='rtm_milanX_sb1st_512_cbz_32.png'
plot_heat_map(map_,title_,xlabel_,ylabel_,xticks,yticks,save_path)

##########################################################################################
##########################################################################################

perf_map_2nd = np.ones((ny,nx))*(float('nan'))
####################### Find best parameters for SB 2nd order #######################
data_summary_2nd = data_summary[data_summary['method'] == 'sb_order2_abc']
data_summary_2nd = data_summary_2nd[data_summary_2nd['grids'] == grid_name]

for i in range(len(data_summary_2nd)):  # Loops from 0 to 4
    if (data_summary_2nd.iloc[i]['grids']==grid_name):
        cb_x_value = data_summary_2nd.iloc[i]['cb_x']
        cb_y_value = data_summary_2nd.iloc[i]['cb_y']
        x_index = np.where(cbx_val == cb_x_value)[0][0]  # [0][0] gets the first index
        y_index = np.where(cby_val == cb_y_value)[0][0]  # [0][0] gets the first index
        perf_map_2nd[y_index,x_index] = data_summary_2nd.iloc[i]['giga_point_s']
        if data_summary_2nd.iloc[i]['giga_point_s']==-10000:
            aa=1
##################
map_=perf_map_2nd
title_='milanX, SB 2nd order, GStencils/sec, grid='+grid_name
xlabel_='Cache blocking in Y dimension, grid points'
ylabel_='Cache blocking in X dimension, grid points'
xticks=cbx_val
yticks=cby_val
save_path='milanX_sb_2nd_1024.png'
plot_heat_map(map_,title_,xlabel_,ylabel_,xticks,yticks,save_path)



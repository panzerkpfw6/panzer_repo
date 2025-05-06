#!/usr/bin/env python
import os,sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import math
####### check logs from Shaheen3
list_of_folders=['/home/plotnips/Dropbox/Apps/Overleaf/2024_pasc_stencil/figures/pavel/grid_perf/Genoa_X_on_intel/logs/test2']
list_of_folders=['./data/par_search/']

sim_info=[]
times=[]
giga_points=[]
giga_flops = []
point_updates = []
flops=[]
method=[]
grids=[]
filenames=[]
cores=[]
cb_size=[]
cb_x=[]
cb_y=[]
cb_z=[]
import glob
####################### Parse files and aggregate them into the dataframe #######################
for data_folder in list_of_folders:
    print(os.getcwd())
    filenames_list=os.listdir(data_folder)
    # filenames_list=glob.glob(os.path.join(data_folder,'*--*') )
    #### filter filenames list
    for index,filename in enumerate(filenames_list):
        # print(index,filename)
        if filename=='log-SB_1st-abc_512_512_512_44_24.log':
            aa=1
        with open(os.path.join(data_folder,filename),'r') as file:
            counter=0
            previous_line=''
            info_line=filename.split('.log')[0]
            if 'SB_1st-abc' in filename:    tmp_method='sb_abc'
            if 'SB_2nd-abc' in filename:    tmp_method='sb_order2_abc'
            if '2023-12-17-13-41-42--0--388022' in filename:
                ss=1
            for line in file:
                if counter==0:
                    start_new_simulation_read=1                    
                    tmp_time=math.nan
                    tmp_giga_points=math.nan
                    tmp_giga_flops=math.nan
                    tmp_point_updates=math.nan
                    tmp_flops=math.nan
                    tmp_grid=math.nan
                    tmp_cb_x=math.nan
                    tmp_cb_y=math.nan
                    tmp_cb_z=math.nan
                if 'velocity size' in line and issubclass(type(tmp_grid),str)==False:
                    pattern = r"velocity size\s*=\s*(\d+)\s*x\s*(\d+)\s*x\s*(\d+)"
                    match = re.search(pattern,line)
                    grid_x = int(match.group(1))
                    grid_y = int(match.group(2))
                    grid_z = int(match.group(3))
                    tmp_grid = match.group(1)+" x "+match.group(2)+" x "+match.group(3)
                elif 'PROP:' in line and np.isnan(tmp_time):
                    tmp1=(line.split('(s)')[0])
                    tmp_time=float(tmp1.split(':PROP:')[-1])
                elif 'Speed:' in line and np.isnan(tmp_giga_points):
                    tmp_giga_points=(line.split('Speed:')[-1])
                    tmp_giga_points=float(tmp_giga_points.split('Mstencils/s')[0])/1000
                elif 'BLOCKX=' in line:
                    pattern = r"BLOCKX=(\d+),\s*BLOCKY=(\d+),\s*BLOCKZ=(\d+)"
                    match = re.search(pattern,line)
                    tmp_cb_x=int(match.group(1))  # Extract BLOCKX value
                    tmp_cb_y=int(match.group(2))
                    tmp_cb_z=int(match.group(3))
                    tmp_cb_size=[tmp_cb_x,tmp_cb_y,tmp_cb_z]
                elif ('[STENCIL MSG]:END' in line):
                    filenames.append(filename)
                    method.append(tmp_method)
                    sim_info.append(info_line)
                    times.append(tmp_time)
                    giga_points.append(tmp_giga_points)
                    grids.append(tmp_grid)
                    cb_size.append(tmp_cb_size)
                    cb_x.append(tmp_cb_x)
                    cb_y.append(tmp_cb_y)
                    cb_z.append(tmp_cb_z)
                    sss=1
                previous_line=line
                counter+=1
        # print(filename)

# Créer un DataFrame à partir des listes
data = pd.DataFrame({
    'filenames':filenames,
    'method':method,
    'sim_info':sim_info,
    'times':times,
    'giga_point_s':giga_points,
    'grids':grids,
    'cb_size':cb_size,
    'cb_x': cb_x,
    'cb_y': cb_y,
    'cb_z': cb_z,
    })
print(data.columns)
######################  find available grids
print(data['grids'])
print(data['method'].unique())
# grids_=['512 x 512 x 512','1024 x 1024 x 512','2048 x 2048 x 512','2048 x 2048 x 1024','2048 x 2048 x 2048']
grids_=['512 x 512 x 512','1024 x 1024 x 512','2048 x 2048 x 512']
print('results will be plotted for grids=',grids_)
data_list_1st_grid=[]
data_list_2nd_grid=[]
metric='giga_point_s'
print(data.columns)

####################### Create summary for dataframe  #######################
data_summary = data[['sim_info', 'grids', 'times', 'giga_point_s', 'cb_size','cb_x','cb_y','cb_z']]
data_summary.to_excel(os.path.join("sb_data.xlsx"))

data_summary = data_summary[data_summary['grids'] == '512 x 512 x 512']
####################### Find best parameters for TB 1st and 2nd order #######################
ss=1
dd=1

cbx_val=np.arange(2, 65, 2);cbx_val = np.insert(cbx_val, 0, 1)
cby_val=np.arange(2, 65, 2);cby_val = np.insert(cby_val, 0, 1)
nx=len(cbx_val)
ny=len(cby_val)
perf_map = np.ones((ny,nx))*(-10000)

data_summary.iloc[0]['sim_info']
print(data_summary['grids'].unique())
len(data_summary)

for i in range(len(data_summary)):  # Loops from 0 to 4
    # if data_summary.iloc[i]['sim_info']=='log-SB_1st-abc_512_512_512_44_24':
    #     aa=1
    if (data_summary.iloc[i]['grids']=='512 x 512 x 512'):
        print(data_summary.iloc[i]['sim_info'])

        cb_x_value = data_summary.iloc[i]['cb_x']
        cb_y_value = data_summary.iloc[i]['cb_y']
        if (cb_x_value==44 and cb_y_value==24):
            aa=1

        x_index = np.where(cbx_val == cb_x_value)[0][0]  # [0][0] gets the first index
        y_index = np.where(cby_val == cb_y_value)[0][0]  # [0][0] gets the first index

        perf_map[y_index,x_index] = 2*data_summary.iloc[i]['giga_point_s']
        if data_summary.iloc[i]['giga_point_s']==-10000:
            aa=1
####################################################
plt.figure()
plt.imshow(perf_map,cmap='viridis',aspect='auto')  # You can change the colormap as needed
plt.colorbar()  # Add a color bar to show the scale of values
plt.title("Heatmap")
plt.show()

####################################################

print(cbx_val[14],cby_val[21])
print(cbx_val[22],cby_val[12])

indices = np.where(perf_map == 0)
perf_map[12,22]
perf_map[21,14]

perf_map[22,12]
perf_map[14,21]
perf_map[1,1]


print(cby_val[22],cbx_val[12])
print(cby_val[14],cbx_val[21])

1150
len(cbx_val)*len(cby_val)
ss=1

# ###### check for missing data among files
for ix in cbx_val:
    for iy in cby_val:
        print(ix,iy)
        if filename=='log-SB_1st-abc_512_512_512_44_24.log':
            aa=1
        filename='log-SB_1st-abc_512_512_512_'+str(ix)+'_'+str(iy)+'.log'
        with open(os.path.join(data_folder,filename),'r') as file:
            tmp=0
            for line in file:
                if 'Speed:' in line:
                    tmp_giga_points=(line.split('Speed:')[-1])
                    tmp=float(tmp_giga_points.split('Mstencils/s')[0])/1000
            if tmp==0:
                aa=1
# ####################################

data_summary

filtered_data_summary = data_summary[data_summary['sim_info'] == 'log-SB_1st-abc_512_512_512_1_1.log']



#!/usr/bin/env python
import os,sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import math
#######################################
####### check logs from Shaheen3
# list_of_folders=['/home/plotnips/Dropbox/Apps/Overleaf/2024_pasc_stencil/figures/pavel/grid_perf/Genoa_X_on_intel/logs/test2']
list_of_folders=['../../logs/test3']

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
            if 'TB_1st_abc' in filename:    tmp_method='tb_abc'
            if 'TB_2nd_abc' in filename:    tmp_method='tb_order2_abc'
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
                if 'velocity size' in line and issubclass(type(tmp_grid),str)==False:
                    pattern = r"velocity size\s*=\s*(\d+)\s*x\s*(\d+)\s*x\s*(\d+)"
                    match = re.search(pattern,line)
                    grid_x = int(match.group(1))
                    grid_y = int(match.group(2))
                    grid_z = int(match.group(3))
                    tmp_grid = match.group(1)+" x "+match.group(2)+" x "+match.group(3)
                elif '[STENCIL MSG]:Total:' in line and '[STENCIL MSG]:Global info:' in previous_line:
                    pattern = r'\[STENCIL MSG\]:Total:\s*([\d.]+)\s*\(s\)'
                    match = re.search(pattern, line)
                    number = match.group(1)
                    tmp_time=float(number)
                elif '[STENCIL MSG]:Total:' in line and '[STENCIL MSG]:Speed info:' in previous_line:
                    pattern = r'\[STENCIL MSG\]:Total:\s*([\d.]+)\s*Mstencils/s'
                    match = re.search(pattern, line)
                    if match:
                        number = match.group(1)
                    tmp_giga_points=float(number)/1000
                elif ('[STENCIL MSG]:END' in line):
                    filenames.append(filename)
                    method.append(tmp_method)
                    sim_info.append(info_line)
                    times.append(tmp_time)
                    giga_points.append(tmp_giga_points)
                    grids.append(tmp_grid)
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
data_summary = data[['sim_info', 'grids', 'times', 'giga_point_s']]
data_summary.to_excel(os.path.join("tb_data.xlsx"))
data_summary = data_summary[data_summary['grids'] == '512 x 512 x 512']
####################### Find best parameters for TB 1st and 2nd order #######################
ss=1
dd=1


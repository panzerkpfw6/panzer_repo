#!/usr/bin/env python
import os,sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import math
#######################################
#### data export from files
# list_of_folders=['./scaling/7773X','./scaling/7763',
# './scaling/Flamingo',
# './scaling/Kanary','./scaling/Shaheen']
####### check logs from Shaheen3
list_of_folders=['/home/plotnips/Dropbox/Apps/Overleaf/2024_pasc_stencil/figures/pavel/grid_perf/Genoa_X_on_intel/logs/test2']
list_of_folders=['../../logs/test2']

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
import glob
####################### Parse files and aggregate them into the dataframe #######################
for data_folder in list_of_folders:
    print(os.getcwd())
    filenames_list=os.listdir(data_folder)
    # filenames_list=glob.glob(os.path.join(data_folder,'*--*') )
    #### filter filenames list
    for index,filename in enumerate(filenames_list):
        # print(index,filename)
        with open(os.path.join(data_folder,filename),'r') as file:
            info_line=filename.split('.log')[0]
            if 'TB_1st_abc' in filename:    tmp_method='tb_abc'
            if 'TB_2nd_abc' in filename:    tmp_method='tb_order2_abc'
            if 'SB_1st-abc' in filename:    tmp_method='sb_abc'
            if 'SB_2nd-abc' in filename:    tmp_method='sb_order2_abc'
            previous_line=''
            if '2023-12-17-13-41-42--0--388022' in filename:
                ss=1
            if '2023-12-01-20-34-09--4--723' in filename:
                ss=1
            for line in file:
                if '********** START SB KERNEL ***********' in line:
                    start_new_simulation_read=1                    
                    tmp_time=math.nan
                    tmp_giga_points=math.nan
                    tmp_giga_flops=math.nan
                    tmp_point_updates=math.nan
                    tmp_flops=math.nan
                    tmp_grid=math.nan
                elif '**** START TB KERNEL FIRST ORDER *****' in line:
                    start_new_simulation_read=1
                    tmp_time=math.nan
                    tmp_giga_points=math.nan
                    tmp_giga_flops=math.nan
                    tmp_point_updates=math.nan
                    tmp_flops=math.nan
                    tmp_grid=math.nan
                elif 'GRID' in line and issubclass(type(tmp_grid), str)==False:
                    if np.isnan(tmp_grid):
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
                elif 'R-' in line and 'R-TB1noABC' in line:
                    tmp_method='tb'
                elif 'R-' in line and 'R-TB1withABC' in line:
                    tmp_method='tb_abc'
                elif '# THREADS' in line:
                    tmp_ncores=(float(line.split()[-1]))
                elif ('********** END SB KERNEL *************' in line) or ('********** END TB KERNEL *************' in line):
                    cores.append(tmp_ncores)
                    filenames.append(filename)
                    sim_info.append(info_line)
                    times.append(tmp_time)
                    giga_points.append(tmp_giga_points)
                    giga_flops.append(tmp_giga_flops)
                    point_updates.append(tmp_point_updates)
                    flops.append(tmp_flops)
                    method.append( tmp_method)
                    grids.append(tmp_grid)
                    # code_type.append(tmp.split('_')[1])
                    # num_th.append( float (tmp.split('_')[2]) )
                    sss=1
                previous_line=line
        # print(filename)

# Créer un DataFrame à partir des listes
data = pd.DataFrame({
    'filenames':filenames,
    'method':method,
    'grids':grids,
    'giga_point_s':giga_points,
    'cores':cores,
    'sim_info':sim_info,
    'times':times,
    'giga_flop_s': giga_flops,
    'point_updates': point_updates,
    'Flops': flops,
    })
print(data.columns)
######################  find available grids
print(data['grids'].unique())
print(data['method'].unique())
print(data['cores'].unique())
print(data['method'].unique())
# grids_=['512 x 512 x 512','1024 x 1024 x 512','2048 x 2048 x 512','2048 x 2048 x 1024','2048 x 2048 x 2048']
grids_=['512 x 512 x 512','1024 x 1024 x 512','2048 x 2048 x 512']
print('results will be plotted for grids=',grids_)
data_list_1st_grid=[]
data_list_2nd_grid=[]
metric='giga_point_s'
print(data.columns)

####################### Create summary for dataframe  #######################
data_summary=data.copy()
data_summary=data_summary.drop(columns=['giga_flop_s','point_updates','Flops','sim_info','times'])
data_summary.to_excel(os.path.join("sb_data.xlsx"))

####################### Find best parameters for TB 1st and 2nd order #######################


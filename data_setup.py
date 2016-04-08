####################################################################################################
###     This setup script prepares the raw data and slices it into smaller data sets.            ###
###     These data sets will be found in the 'data/prepared' subdirectory for further            ###
###     computations.                                                                            ###
###                                                                                              ###
###     Author: Fabian Krautgasser, krautgasser@mpia-hd.mpg.de      JAN/2016                     ###
###     Python Version: >2.*                                                                     ###
####################################################################################################

### DATA SETUP ###
import numpy as np
import os

if os.path.exists('data/raw/PLANET3DP.DATA'):
    with open('data/raw/PLANET3DP.DATA', 'r') as d:
        data_list = d.readlines()

    # data is structured in the following parts: General Information, Grid coordinates I,J,K,
    # additional information, and finally the data set.
    # find the important parts in the file:
    # the next few lines find the index in data_list of the first occurrence of the given string
    grid_i_index = [data_list.index(s) for s in data_list if 'I, R[cm]' in s][0]
    grid_j_index = [data_list.index(s) for s in data_list if 'J, THETA[0-pi]' in s][0]
    grid_k_index = [data_list.index(s) for s in data_list if 'K, PHI[0-2pi]' in s][0]

    add_info_index = [data_list.index(s) for s in data_list if 'I     = radial position' in s][0]
    data_index = [data_list.index(s) for s in data_list if 'Data: I,J,K' in s][0]

    if not os.path.exists('data/raw/general_info.txt'):
        gen_info = data_list[:grid_i_index-1]
        gen_info.append('### Additional Info appended below ###')
        gen_info = gen_info + data_list[add_info_index:data_index]
        with open('data/raw/general_info.txt','wt') as d:
            for i in gen_info:
                d.write(i)

    if not os.path.exists('data/raw/GRID_I_INDEX.txt'):

else: print('ERROR: No data file found!')

####################################################################################################
###     This setup script looks for the data directory and the 'raw' and 'prepared'              ###
###     folders it should contain. if they don't exist, they'll be created.                      ###
###     Also, the plot directory and sequence_plots directory will be created if not available   ###
###                                                                                              ###
###     Author: Fabian Krautgasser, krautgasser@mpia-hd.mpg.de      JAN?2016                     ###
###     Python Version: >2.*                                                                     ###
####################################################################################################

### FOLDER SETUP ###
import os

if not os.path.exists('data'):
    os.mkdir('data')
    os.mkdir('data/raw')
    os.mkdir('data/prepared')
    print("'data' directory and subdirectories 'data/raw' and 'data/prepared' were created,")

if os.path.exists('data'):

    if not os.path.exists('data/raw'):
        os.mkdir('data/raw')
        print("subdirectory 'data/raw' was created.")

    if not os.path.exists('data/prepared'):
        os.mkdir('data/prepared')
        print("subdirectory 'data/prepared' was created.")

if not os.path.exists('plots'):
    os.mkdir('plots')
    os.mkdir('plots/sequence_plots')
    print("'plots' directory and subdirectory 'plots/sequence_plots' were created.")

    if not os.path.exists('plots/sequence_plots'):
        os.mkdir('plots/sequence_plots')
        print("subdirectory 'plots/sequence_plots' was created.")

#!/usr/bin/env python

import sys
import os
import sys
import numpy as np
from prospectTools import *

# INPUTS 1
# defining galaxies to be fit in this slurm run
Galaxies_To_Fit = ['12dam_SDSSPS1','12dam_all']
# Galaxies_To_Fit = ['PS19qp']
Galaxies_To_Fit = ['12dam']


# Command Line TaskID
try:
    commandlineargs = sys.argv
    SLURM_ARRAY_TASK_ID  = int(commandlineargs[1])
except:
    print('No TaskID defined in prospect_and_post.py')

Galaxy = Galaxies_To_Fit[SLURM_ARRAY_TASK_ID]

# INPUTS 2
# defining prospector directory inputs.
# All major adjustments should be done inside new_sne_params.py
APPS = os.environ["APPS"]
objname = Galaxy
outfiledir = os.path.join(APPS,'prospector_alpha/results/slsn/')
outfile = os.path.join(outfiledir,Galaxy)
# datfile = os.path.join(APPS,'prospector_alpha/data/SDSSPS1_vs_all_12dam.json')
# datfile = os.path.join(APPS,'prospector_alpha/data/newtargets.json')
# datfile = os.path.join(APPS,'prospector_alpha/data/J-ApJ-830-13.json')
datfile = os.path.join(APPS,'prospector_alpha/data/J-ApJ-830-13_below1micron.json')


dynasty_cmd = os.path.join(APPS,'prospector_dynesty.py')
param_file = os.path.join(APPS,'new_sne_params.py')

# running prospector with ObjStr
command = 'python ' + dynasty_cmd + ' --param_file='+ param_file +' --objname=' + objname + ' --outfile=' + outfile + ' --datfile=' + datfile
print(command)

os.system(command)


###### PART 2 : Pulling correct objectid and output file from input to code

# pulling files, objectids and dates from output file names
files = os.listdir(outfiledir)
dir_objs_dates = np.array([file.rsplit('_',2)[:2] for file in files])
objs, dates = dir_objs_dates.T

# finding specific object files
object_files_condition = np.where(objs == objname)
object_files = np.array(files)[object_files_condition]

print('files:',files)
print('object files:',object_files)

if len(object_files)==0:
    raise 'no object output files with given objectid'

# finding the newest date to find the most recent model
object_files_dates = np.array(dates)[object_files_condition]
newest_date_arg = np.argmax(object_files_dates)
newest_outfile = object_files[newest_date_arg]

print('newest outfile:',newest_outfile)

outfile = newest_outfile

# # deciding which object output file to use
# if len(object_files)>1:
#     try:
#         # determining if specific version is given to command line
#         outfile_inputed = commandlineargs[2]
#         outfile = outfile_inputed
#
#     except:
#         # otherwise, using newest
#         outfile = newest_outfile
# else:
#     # if there's only 1
#     outfile = object_files[0]

# # chosen output file
# print(outfile, "chosen out of", object_files)







###### PART 3 : Running post-process

# running post-process
# prospector_work_dir = '/Users/lockepatton/Desktop/Research/Berger/Odyssey/RemoteCopy/prospector_home_dir/'
# print(prospector_work_dir)


print('Beginning post-processing. Attempting to create Prospie Object.')

Prospie = postProspect(
    objname=objname,
    objstr=outfile,
    datfile=datfile,
    results_dir=outfiledir,
    verbose=True,
)

Prospie.confirmPrint()
Prospie.postProcessSfrMass()
# Prospie.plotTrace()
# Prospie.plotCorner()
Prospie.postProcessComplete(save=True)
# Prospie.postProcessInit()
Prospie.loadPostProcess()
Prospie.postPostProcess()

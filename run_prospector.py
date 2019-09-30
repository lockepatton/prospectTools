import sys
# from prospectTools import *

# INPUTS 1
# defining galaxies to be fit in this slurm run
Galaxies_To_Fit = ['12dam_SDSSPS1','12dam_all']
# Galaxies_To_Fit = ['PS19qp']


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
outfile = os.path.join(APPS,'prospector_alpha/results/slsn/' + Galaxy)
datfile = os.path.join(APPS,'prospector_alpha/data/SDSSPS1_vs_all_12dam.json')
# datfile = os.path.join(APPS,'prospector_alpha/data/newtargets.json')
# datfile = os.path.join(APPS,'prospector_alpha/data/J-ApJ-830-13.json')

dynasty_cmd = os.path.join(APPS,'prospector_dynesty.py')
param_file = os.path.join(APPS,'new_sne_params.py')

# running prospector with ObjStr
command = 'python ' + dynasty_cmd + ' --param_file='+ param_file +' --objname=' + objname + ' --outfile=' + outfile + ' --datfile=' + datfile
print(command)

os.system(command)

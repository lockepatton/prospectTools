import sys
import os
# from prospectTools import *

# INPUTS 1
# defining galaxies to be fit in this slurm run
# Galaxies_To_Fit = ['12dam_SDSSPS1','12dam_all']
# Galaxies_To_Fit = ['PS19qp']
# Galaxies_To_Fit = ['PTF10hgi', 'SN 2010kd', 'PTF12dam', 'SN 2007bi', 'SN 2011ke', 'SN 2012il', 'PTF11rks', 'SN 2010gx', 'SN 2011kf', 'PTF09cnd', 'SN 2005ap', 'MLS121104:021643+204009', 'PTF09cwl', 'SN 2006oz', 'PTF09atu', 'PS1-12bqf', 'PS1-11ap', 'PS1-10bzj', 'PS1-12zn', 'PS1-11bdn', 'PS1-13gt', 'PS1-10awh', 'PS1-10ky', 'PS1-11aib', 'SCP 06F6', 'PS1-10pm', 'PS1-11tt', 'PS1-10afx', 'PS1-11afv', 'PS1-11bam', 'PS1-12bmy']
Galaxies_To_Fit = ['09as  ', '09uy  ', '09atu ', '09cnd ', '09cwl ', '10bfz ', '10bjp ', '10cwr ', '10fel ', '10heh ', '10hgi ', '10jwd ', '10nmn ', '10qaf ', '10qwu ', '10scc ', '10tpz ', '10uhf ', '10vqv ', '10vwg ', '10yyc ', '10aagc', '11dij ', '11dsf ', '11hrq ', '11rks ', '12dam ', '12epg ', '12gwu ', '12mue ', '12mkp ', '12mxx ']


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
datfile = os.path.join(APPS,'prospector_alpha/data/2020_04_12_J-ApJ-830-13.json')

dynasty_cmd = os.path.join(APPS,'prospector_dynesty.py')
param_file = os.path.join(APPS,'new_sne_params.py')

# running prospector with ObjStr
command = 'python ' + dynasty_cmd + ' --param_file='+ param_file +' --objname=' + objname + ' --outfile=' + outfile + ' --datfile=' + datfile
print(command)

os.system(command)

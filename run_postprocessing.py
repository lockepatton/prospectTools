import os
import sys
import numpy as np

###### PART 1 : Pulling correct objectid and output file from input to code

# PULLING OBJECTNAME FROM COMMANDLINE
commandlineargs = sys.argv
objname  = commandlineargs[1]

# IMPUT DATA DIRECTORY
APPS = os.environ["APPS"]
outfiledir = os.path.join(APPS,'prospector_alpha/results/slsn/')

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

# deciding which object output file to use
if len(object_files)>1:
    try:
        # determining if specific version is given to command line
        outfile_inputed = commandlineargs[2]
        outfile = outfile_inputed

    except:
        # otherwise, using newest
        outfile = newest_outfile
else:
    # if there's only 1
    outfile = object_files[0]

# chosen output file
print(outfile, "chosen out of", object_files)


###### PART 2 : Running post-process


# running post-process
# prospector_work_dir = '/Users/lockepatton/Desktop/Research/Berger/Odyssey/RemoteCopy/prospector_home_dir/'
# print(prospector_work_dir)

prospector_work_dir = APPS

datfile = prospector_work_dir + 'prospector_alpha/data/SDSSPS1_vs_all_12dam.json'
results_dir = prospector_work_dir + 'prospector_alpha/results/slsn/'


Prospie = postProspect(
    objname=objname,
    objstr=outfile,
    datfile=datfile,
    results_dir=results_dir,
    verbose=True,
)

Prospie.confirmPrint()
Prospie.postProcessSfrMass()
# Prospie.plotTrace()
# Prospie.plotCorner()
# Prospie.postProcessInit()
Prospie.postProcessComplete(save=True)
Prospie.loadPostProcess()
Prospie.postPostProcess()
# Prospie.plotSFR(saveplots=True)
# Prospie.plotFilters()
# Prospie.plotSED(plot_SFR=True, saveplots=True)

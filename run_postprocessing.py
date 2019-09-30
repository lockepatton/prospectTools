from myutils import *

import os

# PULLING OBJECTNAME FROM COMMANDLINE
commandlineargs = sys.argv
objname  = commandlineargs[1]

# IMPUT DATA DIRECTORY
APPS = os.environ["APPS"]
outfiledir = os.path.join(APPS,'prospector_alpha/results/slsn/')

# objname = '12dam_SDSSPS1'
# test_path = './prospector_home_dir/prospector_alpha/results/slsn/'

files = os.listdir(outfiledir)
dir_objs_dates = np.array([file.rsplit('_',2)[:2] for file in files])
objs, dates = dir_objs_dates.T

print(files)

objname_outfiles_found = np.where(objs == objname)
newest_outfile_arg = np.argmax(dates[objname_outfiles_found])
newest_outfile = files[newest_outfile_arg]


# USE Newest Prospector Output unless given inside command line 2nd location
try:
    outfile = commandlineargs[2]
else
    outfile = newest_outfile

print(outfile,'newest:', newest_outfile)

#
# # running post-process
# prospector_work_dir = '/Users/lockepatton/Desktop/Research/Berger/Odyssey/RemoteCopy/prospector_home_dir/'
# datfile = prospector_work_dir + 'prospector_alpha/data/SDSSPS1_vs_all_12dam.json'
# results_dir = prospector_work_dir + 'prospector_alpha/results/slsn/'
#
# Prospie = postProspect(
#     objname='12dam_all',
#     objstr='12dam_all_1569002029_mcmc.h5',
#     datfile=datfile,
#     results_dir=results_dir,
#     verbose=True,
# )
#
# Prospie.confirmPrint()
# Prospie.postProcessSfrMass()
# Prospie.plotTrace()
# Prospie.plotCorner()
# Prospie.postProcessInit()
# # Prospie.postProcessComplete(save=True)
# Prospie.loadPostProcess()
# Prospie.postPostProcess()
# Prospie.plotSFR(saveplots=True)
# Prospie.plotFilters()
# Prospie.plotSED(plot_SFR=True, saveplots=True)

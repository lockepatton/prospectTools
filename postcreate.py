from myutils import * 
from prospectTools import *
import sys

# running post-process
prospector_work_dir = '/Users/lockepatton/Desktop/Research/Berger/Odyssey/RemoteCopy/prospector_home_dir/'
datfile = prospector_work_dir + 'prospector_alpha/data/J-ApJ-830-13_below1micron.json'
results_dir = prospector_work_dir + 'prospector_alpha/results/slsn/'

objects = {
    '10fel' : '10fel_1570123003_mcmc.h5',
    '10fel_below1micron' : '10fel_below1micron_1570123002_mcmc.h5',
    '10jwd' : '10jwd_1570123002_mcmc.h5',
    '10jwd_below1micron' : '10jwd_below1micron_1570123002_mcmc.h5',
    '10uhf' : '10uhf_1570123003_mcmc.h5',
    '10uhf_below1micron' : '10uhf_below1micron_1570123002_mcmc.h5',
    '12dam' : '12dam_all_1569002029_mcmc.h5',
    '12dam_below1micron' : '12dam_below1micron_1570123002_mcmc.h5',
#     'AT2018hyz' : 'AT2018hyz_1570205265_mcmc.h5',
}

commandlineargs = sys.argv
objname = commandlineargs[1]
objstr = objects[objname]

Prospie = postProspect(
    objname=objname,
    objstr=objstr,
    datfile=datfile,
    results_dir=results_dir,
    verbose=True,
)

Prospie.confirmPrint()
Prospie.postProcessSfrMass()

Prospie.postProcessInit()
Prospie.postProcessComplete(save=True)

Prospie.loadPostProcess()
print('Loaded post-process')

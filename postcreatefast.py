# from myutils import *
from prospectTools import *
import sys
import os

# Inputs
APPS = os.environ["APPS"]
datfile = os.path.join(APPS,'prospector_alpha/data/2020_04_12_J-ApJ-830-13.json')
results_dir = os.path.join(APPS,'prospector_alpha/results/slsn/')

objstr_all = ['09as_1586824933_mcmc.h5',
               '10fel_1586824941_mcmc.h5',
               '10uhf_1586824941_mcmc.h5',
               '12dam_1586824941_mcmc.h5',
               '09atu_1586824933_mcmc.h5',
               '10heh_1586824941_mcmc.h5',
               '10vqv_1586824941_mcmc.h5',
               '12gwu_1586824934_mcmc.h5',
               '09cnd_1586824938_mcmc.h5',
               '10hgi_1586824941_mcmc.h5',
               '10vwg_1586824941_mcmc.h5',
               '12mkp_1586824934_mcmc.h5',
               '09uy_1586824933_mcmc.h5',
               '10nmn_1586824941_mcmc.h5',
               '10yyc_1586824941_mcmc.h5',
               '12mue_1586824934_mcmc.h5',
               '10aagc_1586824941_mcmc.h5',
               '10qaf_1586824941_mcmc.h5',
               '11dij_1586824941_mcmc.h5',
               '12mxx_1586824934_mcmc.h5',
               '10bfz_1586824933_mcmc.h5',
               '10qwu_1586824941_mcmc.h5',
               '11dsf_1586824941_mcmc.h5',
               '10bjp_1586824933_mcmc.h5',
               '10scc_1586824941_mcmc.h5',
               '11hrq_1586824941_mcmc.h5',
               '10cwr_1586824941_mcmc.h5',
               '10tpz_1586824941_mcmc.h5',
               '11rks_1586824941_mcmc.h5']

get_objname_func = lambda x : x.split('_')[0]

commandlineargs = sys.argv
index = int(commandlineargs[1])

objstr = objstr_all[index]
objname = get_objname_func(objstr)

# Running post-process
print("Running postProspect() for index {}, for object {} and objstr {}.".format(index, objname, objstr))
print("datfile: {}".format(datfile.split('/')[-1]))

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


print('Attempting to plot:')
Prospie.plotTrace(saveplots=True)
Prospie.plotCorner(saveplots=True)
Prospie.plotSED_witherr_nuFnu(saveplots=True, plotnuFnu=False)
Prospie.plotSED_witherr_nuFnu(saveplots=True, plotnuFnu=True)

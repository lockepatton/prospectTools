from myutils import * 
from prospectTools import *
import sys

# running post-process
prospector_work_dir = '/Users/lockepatton/Desktop/Research/Berger/Odyssey/RemoteCopy/prospector_home_dir/'
datfile = prospector_work_dir + 'prospector_alpha/data/J-ApJ-830-13_below1micron.json'
results_dir = prospector_work_dir + 'prospector_alpha/results/slsn/'

objects = {
 '09as': '09as_1586824933_mcmc.h5',
 '10fel': '10fel_1586824941_mcmc.h5',
 '10uhf': '10uhf_1586824941_mcmc.h5',
 '12dam': '12dam_1586824941_mcmc.h5',
 '09atu': '09atu_1586824933_mcmc.h5',
 '10heh': '10heh_1586824941_mcmc.h5',
 '10vqv': '10vqv_1586824941_mcmc.h5',
 '12gwu': '12gwu_1586824934_mcmc.h5',
 '09cnd': '09cnd_1586824938_mcmc.h5',
 '10hgi': '10hgi_1586824941_mcmc.h5',
 '10vwg': '10vwg_1586824941_mcmc.h5',
 '12mkp': '12mkp_1586824934_mcmc.h5',
 '09uy': '09uy_1586824933_mcmc.h5',
 '10nmn': '10nmn_1586824941_mcmc.h5',
 '10yyc': '10yyc_1586824941_mcmc.h5',
 '12mue': '12mue_1586824934_mcmc.h5',
 '10aagc': '10aagc_1586824941_mcmc.h5',
 '10qaf': '10qaf_1586824941_mcmc.h5',
 '11dij': '11dij_1586824941_mcmc.h5',
 '12mxx': '12mxx_1586824934_mcmc.h5',
 '10bfz': '10bfz_1586824933_mcmc.h5',
 '10qwu': '10qwu_1586824941_mcmc.h5',
 '11dsf': '11dsf_1586824941_mcmc.h5',
 '10bjp': '10bjp_1586824933_mcmc.h5',
 '10scc': '10scc_1586824941_mcmc.h5',
 '11hrq': '11hrq_1586824941_mcmc.h5',
 '10cwr': '10cwr_1586824941_mcmc.h5',
 '10tpz': '10tpz_1586824941_mcmc.h5',
 '11rks': '11rks_1586824941_mcmc.h5',
#    '10fel' : '10fel_1570123003_mcmc.h5',
#    '10fel_below1micron' : '10fel_below1micron_1570123002_mcmc.h5',
#    '10jwd' : '10jwd_1570123002_mcmc.h5',
#    '10jwd_below1micron' : '10jwd_below1micron_1570123002_mcmc.h5',
#    '10uhf' : '10uhf_1570123003_mcmc.h5',
#    '10uhf_below1micron' : '10uhf_below1micron_1570123002_mcmc.h5',
#    '12dam' : '12dam_all_1569002029_mcmc.h5',
#    '12dam_below1micron' : '12dam_below1micron_1570123002_mcmc.h5',
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

# prospectTools

A directory to help with post processing and running of prospector. Not designed to be snazzy.

## Basic Setup:

```
from prospectTools import *

prospector_work_dir = './prospector_home_dir/'
datfile = prospector_work_dir + 'prospector_alpha/data/SDSSPS1_vs_all_12dam.json'
results_dir = prospector_work_dir + 'prospector_alpha/results/slsn/'

Prospie = postProspect(
    objname='12dam_SDSSPS1',
    objstr='12dam_SDSSPS1_1569002029_mcmc.h5',
    datfile=datfile,
    results_dir=results_dir,
    verbose=True,
)
```


## Additional Tools:

```
Prospie.confirmPrint()

fig, ax = Prospie.plotTrace()
fig, ax = Prospie.plotCorner()
```


## Post Processing:


#### The first time:

```
# Star formation and Mass
Prospie.postProcessSfrMass()

# Complete SED and Spectra Models
Prospie.postProcessInit()               # Will take ~3 minutes the first call to load fsps libraries into RAM.
Prospie.postProcessComplete(save=True)  # Will take ~40 minutes to model all spectra and create output file

# 1 sigma above and below SED calc
Prospie.postPostProcess()
```

#### Every time after that:

```
Prospie.postProcessSfrMass()
Prospie.postProcessInit()               # Will take ~3 minutes the first call to load fsps libraries into RAM.
Prospie.loadPostProcess()
Prospie.postPostProcess()
```

## Post Plotting:


```
fig, ax = Prospie.plotSFR(plotRecent=True);
fig, ax = Prospie.plotFilters();       # Note necessary to call this in order to know reasonable limits for SFH plot
fig, ax, ins = Prospie.plotSED(plot_SFR=True);
fig, ax = Prospie.plotSED(plot_SFR=False);
```

import prospect.io.read_results as pread
import prospect.io.read_results as rr
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import os
import json
import warnings
import dynesty.plotting
import json

from astropy import constants as const
from astropy import units as u
from scipy import stats
from sedpy import observate
# from mpl_toolkits.axes_grid.inset_locator import inset_axes
from new_sne_params import *
# from new_sne_params_seb_extraSFRbin import *
from prospect.io.read_results import results_from, get_sps
from mpl_toolkits.axes_grid.inset_locator import inset_axes

class postProspect(object):

    def __init__(self,
                 objname='',
                 objstr='',
                 datfile='',
                 results_dir='',
                 out_dir='',
                 plots_dir='./plots/',
                 saveplots=False,
                 verbose=True):

        # output bools
        self.verbose = verbose
        self.saveplots = saveplots

        if self.verbose == True:
            print(objname, objstr, datfile, results_dir, verbose)

        # object specific definitions
        self.objname = objname
        self.objstr = objstr
        self.datfile = datfile
        self.results_dir = results_dir

        # output files
        self.out_dir = out_dir
        self.plots_dir = plots_dir

        # pulling prospector output
        self.res, self.obs, self.mod = rr.results_from(
            self.results_dir + self.objstr)

        # loading model
        self.model = load_model(datfile=self.datfile, objname=self.objname, agelims=[
                                0.0, 7.4772, 8.0, 8.5, 9.0, 9.5, 9.8, 10.0])
        self.model.params['nebemlineinspec'] = True

        # reading input datfile
        with open(self.datfile, 'r') as f:
            self.data = json.load(f)

        # loading sps
        self.sps = load_sps()

        self.wavelength_unitconv = 0.0001 # anstroms to micrometers
        self.flux_unitconv = 3631e-23 # maggies to erg/s/cm^2/Hz
        #1 maggie is the flux density in Janskys divided by 3631
        # u.maggie = (u.Jy/3631).to(u.erg/u.s/u.cm**2/u.Hz).value

        if verbose:
            print('Object Built')

    def confirmPrint(self):

        print('res keys:',self.res.keys(),'\n')
        print('weights shape:',self.res['weights'].shape,'\n')
        # print(self.res['samples_id'])
        print('chain shape:',self.res['chain'].shape,'\n')
        # print(self.res['chain'].T.shape,'\n')
        print('theta labels:',self.res['theta_labels'],'\n')

        # print(self.res)
        # print(self.mod)
        print('obs:', self.obs)

    def postProcessSfrMass(self, p_all = [.50-.341, .50, .50+.341]):

        # SFR

        self.agebins = self.model.params['agebins']

        logsfr_all = [ x for x in self.res['theta_labels'] if "logsfr_ratios" in x ]
        self.sfr_theta = np.array([self.res['chain'].T[self.res['theta_labels'].index(i)] for i in logsfr_all])
        self.massmet_0 = self.res['chain'].T[self.res['theta_labels'].index('massmet_1')]
        self.massmet_1 = self.res['chain'].T[self.res['theta_labels'].index('massmet_2')]

        self.Masses = np.array([logmass_to_masses(massmet=[massmet_0_it, massmet_1_it],
                                                  logsfr_ratios=sfr_it,
                                                  agebins=self.agebins)
                                for sfr_it, massmet_0_it, massmet_1_it in zip(self.sfr_theta.T,
                                                                              self.massmet_0,
                                                                              self.massmet_1)])


        # STELLAR MASS

        self.p_all = p_all
        self.dt = (10**self.agebins[:,1]-10**self.agebins[:,0]) # agebins in log(yr)

        self.stellarmass = dynesty.plotting._quantile(10**self.massmet_0, self.p_all, weights=self.res['weights'])

        if self.verbose:
            print('stellar mass:')
            print('median, dlow, dhi')
            print ('{:.2e}, {:.2e}, {:.2e}'.format(self.stellarmass[1],
                                                   -self.stellarmass[1]+self.stellarmass[0],
                                                   -self.stellarmass[1]+self.stellarmass[2]))
            print('low, middle, hi')
            print('{:.2e}, {:.2e}, {:.2e}'.format(self.stellarmass[0],
                                                  self.stellarmass[1],
                                                  self.stellarmass[2]))


        # METALLICITY

        self.metallicity = dynesty.plotting._quantile(10**self.massmet_1, self.p_all, weights=self.res['weights'])

        if self.verbose:
            print('metallicity:')
            print('median, dlow, dhi')
            print ('{:.2e}, {:.2e}, {:.2e}'.format(self.metallicity[1],
                                                   -self.metallicity[1]+self.metallicity[0],
                                                   -self.metallicity[1]+self.metallicity[2]))
            print('low, middle, hi')
            print('{:.2e}, {:.2e}, {:.2e}'.format(self.metallicity[0],
                                                  self.metallicity[1],
                                                  self.metallicity[2]))
            print('median, dlow, dhi')
            print(self.metallicity[1], -self.metallicity[1]+self.metallicity[0], -self.metallicity[1]+self.metallicity[2])
            print('low, middle, hi')
            print(self.metallicity[0], self.metallicity[1], self.metallicity[2])



        # LOCAL SFR

        def calcLocalSFR(mass):
            t0, sfr0 = self.agebins[:,1][0],(mass[0]/self.dt)[0]
            t1, sfr1 = self.agebins[:,1][1],(mass[1]/self.dt)[1]
            timeScaledAverage = ((10**t0-0) * sfr0 + (10**t1-10**t0) * sfr1)/10**t1
            return timeScaledAverage

        self.LocalSFRs = np.array([calcLocalSFR(mass) for mass in self.Masses])
        percentLocalSFR = lambda LocalSFRs, p_all, weights : np.array([dynesty.plotting._quantile(LocalSFRs, p_all, weights=weights)])

        self.LocalSFRlo, self.LocalSFRmed, self.LocalSFRhi = percentLocalSFR(self.LocalSFRs, self.p_all, self.res['weights'])[0]

        if self.verbose:
            print('local SFR:')
            print('median, dlow, dhi')
            print ('{:.2e}, {:.2e}, {:.2e}'.format(self.LocalSFRmed,
                                                   -1*(self.LocalSFRmed-self.LocalSFRlo),
                                                   -1*(self.LocalSFRmed-self.LocalSFRhi)))
            print('low, middle, hi')
            print ('{:.2e}, {:.2e}, {:.2e}'.format(self.LocalSFRlo, self.LocalSFRmed, self.LocalSFRhi))


        # SFR History / Masses

        percentMass = lambda Masses, p_all, weights : np.array([dynesty.plotting._quantile(mass, p_all, weights=weights) for mass in Masses.T])
        self.percentMasses =  percentMass(self.Masses, self.p_all, self.res['weights']) # Masses Msolar at -1, 0, 1 sigmas

        self.t_uni = max(self.agebins[:,1])

        if self.verbose:
            print('log10(age of universe):',self.t_uni)

        # Maximum Probability Case

        imax = np.argmax(self.res['lnprobability'])
        i = np.unravel_index(imax, self.res['lnprobability'].shape)
        self.theta_max = self.res['chain'][i, :].copy()[0]


    def plotTrace(self, saveplots=False):

        tracefig = rr.traceplot(self.res, linewidth=.05)
        plt.suptitle(self.objname)

        if saveplots:
            plt.savefig(os.path.join(self.plots_dir, self.objstr + "_trace.pdf"))

        plt.close()

        return tracefig


    def plotCorner(self, saveplots=False):

        p_plottingrange = [.01,.99] #[.25,.75]
        range = [dynesty.plotting._quantile(theta_it, p_plottingrange, weights=self.res['weights']) for theta_it in self.res['chain'].T]

        cornerfig = rr.subcorner(self.res, start=0, thin=5, truths=self.theta_max, range=range)
        # , fig=plt.subplots(5,5,figsize=(27,27))[0]

        if self.verbose:
            print('MAP value: {}'.format(self.theta_max))

        if saveplots:
            plt.savefig(os.path.join(self.plots_dir, self.objstr + "_corner.pdf"))

        plt.close()

        return cornerfig


    def postProcessInit(self):
        self.theta = self.model.theta.copy()
        self.mspec, self.mphot, mextra = self.model.mean_model(self.theta, self.obs, sps=self.sps)
        self.mspec_map, self.mphot_map, _ = self.model.mean_model(self.theta_max, self.obs, sps=self.sps)


    def postProcessComplete(self, save=True):
        halpha_line_i = 62
        total_len = len(self.res['chain'][:, 0].copy())

        self.AllModelPull = {}
        self.AllModelPull['mspec'] = []
        self.AllModelPull['mphot'] = []
        self.AllModelPull['nebline_lum'] = []
        self.AllModelPull['nebline_lum_halpha'] = []

        for i in np.arange(total_len):

            theta_max_test = self.res['chain'][(i,), :].copy()[0]
            mspec, mphot, mextra = self.model.mean_model(theta_max_test, self.obs, sps=self.sps)
            self.AllModelPull['mspec'].append(mspec)
            self.AllModelPull['mphot'].append(mphot)
            nebline_lum = np.array(self.sps.get_nebline_luminosity)
            self.AllModelPull['nebline_lum'].append(nebline_lum)
            self.AllModelPull['nebline_lum_halpha'].append(nebline_lum[halpha_line_i])

        if save:
            np.savez_compressed(os.path.join(self.out_dir, self.objstr+'-mspec'),self.AllModelPull['mspec'])
            AllModelPull_NoSpec = self.AllModelPull.copy()
            del AllModelPull_NoSpec['mspec']
            np.savez_compressed(os.path.join(self.out_dir, self.objstr), AllModelPull_NoSpec)

    def loadPostProcess(self, outputstr='.npz'):
        AllModelPull_Loaded = np.load(os.path.join(self.out_dir, self.objstr+outputstr), allow_pickle = True)
        self.AllModelPull = AllModelPull_Loaded['arr_0'].item()

        mspec_Loaded = np.load(os.path.join(self.out_dir, self.objstr+'-mspec'+outputstr), allow_pickle = True)
        self.AllModelPull['mspec'] = mspec_Loaded['arr_0']

    def postPostProcess(self):
        # Test Adding (See if this means I don't have to run postProcessInit)
        self.theta = self.model.theta.copy()

        self.percentMspec = lambda mspecT, p_all, weights : np.array([dynesty.plotting._quantile(i, p_all, weights=weights)
                                                         for i in mspecT])

        self.MspecsT = np.array(self.AllModelPull['mspec']).T
        returns = self.percentMspec(self.MspecsT, self.p_all, self.res['weights'])
        self.siglo, self.med, self.sighi = returns.T


        self.percentMphot = lambda mphotT, p_all, weights: np.array([dynesty.plotting._quantile(i, p_all, weights=weights)
                                                                 for i in mphotT])

        self.MphotT = np.array(self.AllModelPull['mphot']).T
        returns = self.percentMphot(self.MphotT, self.p_all, self.res['weights'])

        self.mphot_siglo, self.mphot_med, self.mphot_sighi = returns.T
        self.mphot_siglo_err = np.abs(self.mphot_med - self.mphot_siglo)
        self.mphot_sighi_err = np.abs(self.mphot_med - self.mphot_sighi)

    def save1SigmaSpecFile(self, ):
        self.a = 1.0 + self.model.params.get('zred', 0.0) # cosmological redshifting
        self.wspec = self.sps.wavelengths.copy()
        self.wspec *= self.a

        ToSave = pd.DataFrame({
            'wavelength[A]': self.wspec,
            '1sig_low[maggies]': self.siglo,
            '1sig_hi[maggies]': self.sighi,
            'median[maggies]': self.med,
            '1sig_low[erg/s/cm2/Hz]': self.siglo*self.flux_unitconv,
            '1sig_hi[erg/s/cm2/Hz]': self.sighi*self.flux_unitconv,
            'median[erg/s/cm2/Hz]': self.med*self.flux_unitconv,
        })
        ToSave.to_csv(os.path.join(self.out_dir, self.objstr+'_lowmedhi.txt'))

    def save1SigmaPhotFile(self, ):

        self.wphot = self.obs['wave_effective'].copy()

        ToSave = pd.DataFrame({
            'wavelength[A]': self.wphot,
            '1sig_low[maggies]': self.mphot_siglo,
            '1sig_hi[maggies]': self.mphot_sighi,
            'median[maggies]': self.mphot_med,
            '1sig_low[erg/s/cm2/Hz]': self.mphot_siglo*self.flux_unitconv,
            '1sig_hi[erg/s/cm2/Hz]': self.mphot_sighi*self.flux_unitconv,
            'median[erg/s/cm2/Hz]': self.mphot_med*self.flux_unitconv,
        })
        ToSave.to_csv(os.path.join(self.out_dir, self.objstr+'_mphot_lowmedhi.txt'))

    def save1SigmaSFRFile(self, ):

        self.SFR_all = {}
        self.SFR_all['time[Gyr]'] = []
        self.SFR_all['agebins_0'] = self.agebins[:,0]
        self.SFR_all['agebins_1'] = self.agebins[:,1]
        self.SFR_all['M_lo[Mdot]'] = []
        self.SFR_all['M_med[Mdot]'] = []
        self.SFR_all['M_hi[Mdot]'] = []
        self.SFR_all['SFR_lo[Mdot/yr]'] = []
        self.SFR_all['SFR_med[Mdot/yr]'] = []
        self.SFR_all['SFR_hi[Mdot/yr]'] = []

        for timestart, timeend, mass in zip(self.agebins[:,0],self.agebins[:,1],self.percentMasses):
            dt_ = 10**timeend - 10**timestart

            self.SFR_all['time[Gyr]'].append(dt_)
            self.SFR_all['M_lo[Mdot]'].append(mass[0])
            self.SFR_all['M_med[Mdot]'].append(mass[1])
            self.SFR_all['M_hi[Mdot]'].append(mass[2])
            self.SFR_all['SFR_lo[Mdot/yr]'].append(mass[0]/dt_)
            self.SFR_all['SFR_med[Mdot/yr]'].append(mass[1]/dt_)
            self.SFR_all['SFR_hi[Mdot/yr]'].append(mass[2]/dt_)

        ToSave = pd.DataFrame(self.SFR_all)
        ToSave.to_csv(os.path.join(self.out_dir, self.objstr+'_SFR_lowmedhi.txt'))

    def plotSFR(self, figax=None, saveplots=True, plotRecent=True,
                fs = 21, fs2 = 17, fs3 = 12,
                c1 = 'lightgrey', c2 = 'darkgrey',
                time_units = 1e-9,
                color_recent = '#38A6A5'):


        if figax is None:
            fig, ax = plt.subplots(1,1, figsize=(6,6))
        else:
            fig, ax = figax

        dt_last = None
        for timestart, timeend, mass in zip(self.agebins[:,0],self.agebins[:,1],self.percentMasses):

            dt_ = 10**timeend - 10**timestart

            if dt_last != None:
                plt.loglog(time_units*10**np.array([timestart, timestart]), [mass_last[1]/dt_last, mass[1]/dt_],
                           lw=0.7, c=c2)
            plt.loglog(time_units*10**np.array([timestart,timeend]), [mass[1],mass[1]]/dt_, lw=0.7, c=c2)
            plt.fill_between(time_units*10**np.array([timestart,timeend]), [mass[0],mass[0]]/dt_, [mass[2],mass[2]]/dt_,
                             lw=0.7, color=c1, alpha=.5)

            mass_last = mass
            dt_last = dt_

        plt.loglog()
        ax.set_xlabel(r't $\rm ( Gyr )$', fontsize=fs2)
        ax.set_ylabel(r'$\rm SFR_{model}$ $\rm (M_\odot \,\, yr^{-1})$   ', fontsize=fs2);
        ax.set_xlim(time_units*10**7,time_units*10**self.t_uni)

        # plt.plot(np.linspace(0,time_units*10**8,100),[timeScaledAverage]*100, 'r--')

        if plotRecent:
            # plt.plot(np.linspace(0,time_units*10**8,100),[self.LocalSFRlo]*100, color=color_recent, linestyle='--', lw=1)
            plt.plot(np.linspace(0,time_units*10**8,100),[self.LocalSFRmed]*100, color=color_recent, lw=1,
                     label='Last 100 Myr')
            # plt.plot(np.linspace(0,time_units*10**8,100),[LocalSFRhi]*100, color=color_recent, linestyle='--', lw=1)

            plt.fill_between(np.linspace(0,time_units*10**8,100), [self.LocalSFRlo]*100, [self.LocalSFRhi]*100,
                             color=color_recent, alpha=.1, label='$\pm1\sigma$ Model')

            ax.legend(loc=4, fontsize=fs3, frameon=True)

        ax.invert_xaxis()
        ax.tick_params(which='major', direction='in', labelsize=fs3)
        ax.tick_params(which='minor', direction='in', labelsize=fs3)
        ax.tick_params(top=True, bottom=True, left=True, right=True, which='both')
        ax.minorticks_on()

        if saveplots:
            plt.savefig(os.path.join(self.plots_dir,'SFR'+self.objstr+'.pdf'))

        return fig, ax

    def plotFilters(self, saveplots=True):
        fig, ax = plt.subplots(1,1, figsize=(18,2.5))
        self.n_obs = len(self.obs['filters'])
        obscolors = cm.jet(np.linspace(0,1,self.n_obs))
        ymin, ymax = 0, 1

        # plot transmission curves
        for f,c in zip(self.obs['filters'],obscolors):
            w, t = f.wavelength.copy(), f.transmission.copy()
            t *= 1/t.max()
        #     while t.max() > 1:
        #         t /= 10.
            t = 0.01*(ymax-ymin)*t + ymin
            ax.semilogx(w, t, lw=1.2, color=c, alpha=0.7)

            filtername = ' '.join(f.name.split('_'))
            ax.text(f.wave_effective, min(t), filtername,
                     rotation=90, color=c, horizontalalignment="left", verticalalignment="bottom")


        self.xmin_filters, self.xmax_filters = ax.get_xlim()

        if saveplots:
            plt.savefig(os.path.join(self.plots_dir,self.objstr+"_filters.pdf"))

        return fig, ax


    def plotSED(self, figax = None, saveplots=True, savedir='./', plot_SFR=True,
                fs = 21, fs2 = 17, fs3 = 12, fs_ticks = 19,
                colors = ['#DAA51B', '#38A6A5', '#2F8AC4', '#5F4690'],
                ins_width=2, ins_height=2,
                ins_bbox_to_anchor=[.25, 0.5, 600+30, 1000/2-215 - 5],
                ins_loc=1, time_units=1e-9, c1='lightgrey', c2='darkgrey'):

        self.wphot = self.obs['wave_effective'].copy()

        self.a = 1.0 + self.model.params.get('zred', 0.0) # cosmological redshifting
        self.wspec = self.sps.wavelengths.copy()
        self.wspec *= self.a

        self.wavelength_unitconv = 0.0001 # anstroms to micrometers
        self.flux_unitconv = 3631e-23 # maggies to erg/s/cm^2


        if figax is None:
            fig, ax = plt.subplots(1,1, figsize=(10,7))
        else:
            fig, ax = figax

        color_model = colors[2]

        # ax.errorbar(self.wphot, self.mphot_map, label='Model photometry (MAP)',
        #             marker='s', markersize=8, alpha=0.8, ls='',
        #             markerfacecolor='none', markeredgecolor=color_model,
        #             markeredgewidth=3)

        # REMOVED THIS ONE
        # ax.plot(self.wphot*self.wavelength_unitconv, self.mphot_map*self.flux_unitconv, alpha=0)

        # ax.scatter(self.wphot*self.wavelength_unitconv, self.mphot_map*self.flux_unitconv, c=color_model, marker='s',
        #            facecolor='none',label='Model photometry | max probability',
        #            zorder=2)

        # print(xmin, xmax)

        self.n_obs = len(self.obs['filters'])
        self.obscolors = cm.rainbow_r(np.linspace(0,1,self.n_obs))
        # obscolors = cm.Set3(np.linspace(0,1,self.n_obs))

        # obscolors = 'blue'
        # obscolors = obs_all
        ax.errorbar(self.wphot*self.wavelength_unitconv, self.obs['maggies']*self.flux_unitconv, yerr=self.obs['maggies_unc']*self.flux_unitconv,
                    ls='',  c='white',
                    markeredgecolor='black', markerfacecolor='none',
                    zorder=3)
        ax.scatter(self.wphot*self.wavelength_unitconv, self.obs['maggies']*self.flux_unitconv, c='white',
                   marker='o', facecolor='none',
                   zorder=3)
        ax.errorbar(self.wphot*self.wavelength_unitconv, self.obs['maggies']*self.flux_unitconv, yerr=self.obs['maggies_unc']*self.flux_unitconv,
                    ls='', color=self.obscolors, alpha=.7,
                    markeredgecolor='black', markerfacecolor='none',
                    zorder=3)
        ax.scatter(self.wphot*self.wavelength_unitconv, self.obs['maggies']*self.flux_unitconv, c=self.obscolors,  alpha=.7,
                   marker='o', facecolor='none', label='Observed Photometry',
                   zorder=3)

        # ax.errorbar(self.wphot, self.obs['maggies'], yerr=self.obs['maggies_unc'],
        #             label='Observed photometry', ecolor=self.obscolors,
        #             marker='o', markersize=8, ls='', lw=3, alpha=0.8,
        #             markerfacecolor='none', markeredgecolor='red',
        #             markeredgewidth=3)

        # ax.loglog(self.wspec*self.wavelength_unitconv, mspec_map*self.flux_unitconv, label='Model spectrum | max probability',
        #           lw=0.7, color=color_model, alpha=0.7, zorder=0)


        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()

        # ax.plot(self.wspec*self.wavelength_unitconv, siglo*self.flux_unitconv, lw=0.7, c='lightgrey')
        ax.plot(self.wspec*self.wavelength_unitconv, self.med*self.flux_unitconv, lw=0.7, c='darkgrey', zorder=1)
        # ax.plot(self.wspec*self.wavelength_unitconv, sighi*self.flux_unitconv, lw=0.7, c='lightgrey')
        ax.fill_between(self.wspec*self.wavelength_unitconv, self.siglo*self.flux_unitconv, self.sighi*self.flux_unitconv,
                        label='$\pm1\sigma$ Model', lw=0.7, color='lightgrey', alpha=.5, zorder=-1)

        asymmetric_error = np.array([self.mphot_siglo_err, self.mphot_sighi_err])
        ax.errorbar(self.wphot*self.wavelength_unitconv, self.mphot_med*self.flux_unitconv,
                    yerr=asymmetric_error*self.flux_unitconv,
                    c='slategrey', ecolor='slategrey', fmt='s', alpha=.7, label='Model Photometry')

        ax.set_xlim(self.xmin_filters*self.wavelength_unitconv, self.xmax_filters*self.wavelength_unitconv)

        if figax == None:
            dy_filterspace = 10**-1.2*ymax-ymin # Changed this
            ax.set_ylim(ymin-dy_filterspace, ymax)
        else:
            ax.set_ylim(ymin, ymax)


        # print(len(obscolors))
        # colors_short = ['#DAA51B', '#5F4690'] #  '#38A6A5',
        # colors_short = ['#DAA51B', '#38A6A5', '#5F4690']


        ax2 = ax.twinx()

        obs_all = []
        # plot transmission curves
        for f,c in zip(self.obs['filters'],self.obscolors):
        # for i,f in enumerate(self.obs['filters']):
        #     i = i % 3
        #     c = colors_short[i]
            w, t = f.wavelength.copy(), f.transmission.copy()
            t *= 1/t.max()
            t *=.06
            ax2.plot(w*self.wavelength_unitconv, t, lw=1, color=c, alpha=.65)
            obs_all.append(c)


        ax2.set_ylim(0,1)
        ax2.axes.get_yaxis().set_visible(False)

        # ax.set_xlim(10**3.3*self.wavelength_unitconv, 10**5*self.wavelength_unitconv)
        ax.set_xlim(self.xmin_filters*self.wavelength_unitconv, self.xmax_filters*self.wavelength_unitconv)

        ax.set_xlabel(r'Wavelength $\rm (\mu m)$', fontsize=fs)
        ax.set_ylabel(r'Flux Density $\rm (erg \,\, s^{-1} \,\, cm^{-2} \,\, Hz^{-1})$', fontsize=fs)
        ax.legend(loc=2, fontsize=fs3, frameon=True)

        ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
        # plt.xticks([9000,10000,11000])
        print(self.obs['filters'])
        # plt.xticks([min(f.wavelength)*self.wavelength_unitconv for f in self.obs['filters']])
        # plt.xticks([.3,1,3])

        ax.tick_params(which='major', direction='in', labelsize=fs_ticks)
        ax.tick_params(which='minor', direction='in', labelsize=fs_ticks)
        ax.tick_params(top=True, bottom=True, left=True, right=True, which='both')
        ax.minorticks_on()

        returns = [fig, ax]

        if plot_SFR:
            ins = inset_axes(ax,
                             width=ins_width,
                             height=ins_height,
                             bbox_to_anchor=ins_bbox_to_anchor,
                             loc=ins_loc)
            # [left, bottom, width, height], or a tuple of [left, bottom]


            c1 = 'lightgrey'
            c2 = 'darkgrey'

            dt_last = None
            for timestart, timeend, mass in zip(self.agebins[:,0],self.agebins[:,1],self.percentMasses):

                dt_ = 10**timeend - 10**timestart

                if dt_last != None:
                    ins.loglog(time_units*10**np.array([timestart, timestart]), [mass_last[1]/dt_last, mass[1]/dt_],
                               lw=0.7, c=c2)
                ins.loglog(time_units*10**np.array([timestart,timeend]), [mass[1],mass[1]]/dt_, lw=0.7, c=c2)
                ins.fill_between(time_units*10**np.array([timestart,timeend]), [mass[0],mass[0]]/dt_, [mass[2],mass[2]]/dt_,
                                 lw=0.7, color=c1, alpha=.5)

                mass_last = mass
                dt_last = dt_

            ins.loglog()
            ins.set_xlabel(r't $\rm ( Gyr )$', fontsize=fs2)
            ins.set_ylabel(r'$\rm SFR_{model}$ $\rm (M_\odot \,\, yr^{-1})$   ', fontsize=fs2);
            ins.set_xlim(time_units*10**7,time_units*10**self.t_uni)

            # plt.plot(np.linspace(0,time_units*10**8,100),[timeScaledAverage]*100, 'r--')

            color_recent = colors[1]

            # plt.plot(np.linspace(0,time_units*10**8,100),[self.LocalSFRlo]*100, color=color_recent, linestyle='--', lw=1)
            plt.plot(np.linspace(0,time_units*10**8,100),[self.LocalSFRmed]*100, color=color_recent, lw=1,
                     label='Last 100 Myr')
            # plt.plot(np.linspace(0,time_units*10**8,100),[LocalSFRhi]*100, color=color_recent, linestyle='--', lw=1)

            ins.fill_between(np.linspace(0,time_units*10**8,100), [self.LocalSFRlo]*100, [self.LocalSFRhi]*100,
                             color=color_recent, alpha=.1)

            ins.legend(loc=4, fontsize=fs3, frameon=True)
            ins.invert_xaxis()
            ins.tick_params(which='major', direction='in', labelsize=fs3)
            ins.tick_params(which='minor', direction='in', labelsize=fs3)
            ins.tick_params(top=True, bottom=True, left=True, right=True, which='both')
            ins.minorticks_on()

            returns = [fig, ax, ins_height]

        if saveplots:
            plt.savefig(os.path.join(self.plots_dir,self.objstr+"_modelspec_phot_newshape_zoom.pdf"))

        return returns



    def plotSED_witherr(self, figax = None, saveplots=True, savedir='./',
                        plot_SFR=True, plotFilterNames=True,
                        fs = 21, fs2 = 17, fs3 = 12, fs_ticks = 12, figsize=(10,7),
                        colors = ['#DAA51B', '#38A6A5', '#2F8AC4', '#5F4690'],
#                         ax_space = (0, .1, 1, 1),
#                         ax2_space = (.7, .25, .35, .35),
#                         ax3_space = (0, 0, 1, .1),
                        ax_space = (.1, .2, .8, .7),
                        ax2_space = (.67, .57, .3, .3),
                        ax3_space = (.1, .1, .8, .1),
                        time_units=1e-9, c1='lightgrey', c2='darkgrey'):

        self.wphot = self.obs['wave_effective'].copy()

        self.a = 1.0 + self.model.params.get('zred', 0.0) # cosmological redshifting
        self.wspec = self.sps.wavelengths.copy()
        self.wspec *= self.a

        # If including SFR subplot:
        if plot_SFR:
            if figax is None:
                fig = plt.figure(figsize=figsize)
                ax=fig.add_axes(ax_space)
                ax2_space = np.array(ax2_space)
                ax2_space[2] = ax2_space[2]*figsize[1]/figsize[0]
                ax2=fig.add_axes(tuple(ax2_space))
                ax3=fig.add_axes(ax3_space)
            else:
                fig, [ax,ax2,ax3] = figax
        # If not including SFR subplot
        else:
            if figax is None:
                fig = plt.figure(figsize=figsize)
                ax=fig.add_axes(ax_space)
                ax3=fig.add_axes(ax3_space)
            else:
                fig, [ax,ax3] = figax

        color_model = colors[2]

        ax.tick_params(which='major', direction='in', labelsize=fs_ticks)
        ax.tick_params(which='minor', direction='in', labelsize=fs_ticks)
        ax.tick_params(top=True, bottom=True, left=True, right=True, which='both')
        ax.minorticks_on()
        ax.get_xaxis().set_visible(False)

        # REMOVED THIS ONE
        # ax.plot(self.wphot*self.wavelength_unitconv, self.mphot_map*self.flux_unitconv, alpha=0)

        self.n_obs = len(self.obs['filters'])
        self.obscolors = cm.rainbow_r(np.linspace(0,1,self.n_obs))
        # self.obscolors = cm.Set3(np.linspace(0,1,self.n_obs))


        if figax == None:
            label='Observed Photometry'
        else:
            label=''

        ax.errorbar(self.wphot*self.wavelength_unitconv, self.obs['maggies']*self.flux_unitconv, yerr=self.obs['maggies_unc']*self.flux_unitconv,
                    ls='',  c='white',
                    markeredgecolor='black', markerfacecolor='none',
                    zorder=3)
        ax.scatter(self.wphot*self.wavelength_unitconv, self.obs['maggies']*self.flux_unitconv, c='white',
                   marker='o', facecolor='none',
                   zorder=3)
        ax.errorbar(self.wphot*self.wavelength_unitconv, self.obs['maggies']*self.flux_unitconv, yerr=self.obs['maggies_unc']*self.flux_unitconv,
                    ls='', color=self.obscolors, alpha=.7,
                    markeredgecolor='black', markerfacecolor='none',
                    zorder=3)
        ax.scatter(self.wphot*self.wavelength_unitconv, self.obs['maggies']*self.flux_unitconv, c=self.obscolors,  alpha=.7,
                   marker='o', facecolor='none', label=label,
                   zorder=3)

        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()


        if figax == None:
            label='Model Photometry'
        else:
            label=''

        asymmetric_error = np.array([self.mphot_siglo_err, self.mphot_sighi_err])
        ax.errorbar(self.wphot*self.wavelength_unitconv, self.mphot_med*self.flux_unitconv,
                    yerr=asymmetric_error*self.flux_unitconv,
                    c='slategrey', ecolor='slategrey', fmt='s', alpha=.7, label=label)

        if figax == None:
            label='$\pm1\sigma$ Model'
        else:
            label=''

        ax.plot(self.wspec*self.wavelength_unitconv, self.med*self.flux_unitconv, lw=0.7, c=c2, zorder=1)
        ax.fill_between(self.wspec*self.wavelength_unitconv, self.siglo*self.flux_unitconv, self.sighi*self.flux_unitconv,
                        label=label, lw=0.7, color=c1, alpha=.5, zorder=-1)

        ax.set_xlim(self.xmin_filters*self.wavelength_unitconv, self.xmax_filters*self.wavelength_unitconv)

        dy_filterspace = 10**-1.2*ymax-ymin # Changed this
        ax.set_ylim(ymin-dy_filterspace, ymax)

        # Plotting filter curves on twin axis
        axtwin = ax.twinx()

        obs_all = []
        # plot transmission curves
        for f,c in zip(self.obs['filters'],self.obscolors):
            w, t = f.wavelength.copy(), f.transmission.copy()
            t *= 1/t.max()
            t *=.06
            axtwin.plot(w*self.wavelength_unitconv, t, lw=1, color=c, alpha=.65)


            filtername = ' '.join(f.name.split('_'))
            if plotFilterNames:
                axtwin.text(f.wave_effective*self.wavelength_unitconv, max(t), filtername,
                            rotation=90, color=c, horizontalalignment="left", verticalalignment="bottom")

            obs_all.append(c)

        axtwin.set_ylim(0,1)
        axtwin.axes.get_yaxis().set_visible(False)

        ax.set_xlim(self.xmin_filters*self.wavelength_unitconv, self.xmax_filters*self.wavelength_unitconv)
        axtwin.set_xlim(self.xmin_filters*self.wavelength_unitconv, self.xmax_filters*self.wavelength_unitconv)

        ax.set_xlabel(r'Wavelength $\rm (\mu m)$', fontsize=fs)
        ax.set_ylabel(r'Flux Density $\rm (erg \,\, s^{-1} \,\, cm^{-2} \,\, Hz^{-1})$', fontsize=fs)
        ax.legend(loc=2, fontsize=fs3, frameon=True)

        ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

        ax.tick_params(which='major', direction='in', labelsize=fs_ticks)
        ax.tick_params(which='minor', direction='in', labelsize=fs_ticks)
        ax.tick_params(top=True, bottom=True, left=True, right=True, which='both')
        ax.minorticks_on()

        # MAGNITUDE DIFFERENCE LOWER PLOT
        ax3.tick_params(which='major', direction='in', labelsize=fs_ticks)
        ax3.tick_params(which='minor', direction='in', labelsize=fs_ticks)
        ax3.tick_params(top=True, bottom=True, left=True, right=True, which='both')
        ax3.minorticks_on()

        # ax3.axhline(y=0, color='k', linestyle='--')

        diff_obs_mod = self.obs['maggies']*self.flux_unitconv - self.mphot_med*self.flux_unitconv
        yerr_obs = self.obs['maggies_unc']*self.flux_unitconv
        yerr_mod = np.abs(self.mphot_siglo - self.mphot_sighi)
        diff_obs_mod_err = np.sqrt(yerr_obs**2 + yerr_mod**2)

        ax3.plot(self.wphot*self.wavelength_unitconv, diff_obs_mod,
                 marker='o', alpha=0, linewidth=0)
        ax3.scatter(self.wphot*self.wavelength_unitconv, diff_obs_mod,
                    c=self.obscolors, alpha=.7, marker='o', facecolor='none',
                    zorder=3)

        ax3_ylims = ax3.get_ylim()
        limit = max(np.abs(ax3_ylims))*1.3
        ax3.set_ylim(-limit,limit)

        ax3.set_xlim(self.xmin_filters*self.wavelength_unitconv,
                     self.xmax_filters*self.wavelength_unitconv)

        ax3.axhline(y=0, color='k', linestyle='--', alpha=0.2)

        ax3.set_xlabel(r'Wavelength $\rm (\mu m)$', fontsize=fs)
        ax3.set_ylabel(r'$\Delta$', fontsize=fs);

        returns = [fig, ax, ax3]

        if plot_SFR:

            ax2.tick_params(which='major', direction='in', labelsize=fs_ticks)
            ax2.tick_params(which='minor', direction='in', labelsize=fs_ticks)
            ax2.tick_params(top=True, bottom=True, left=True, right=True, which='both')
            ax2.minorticks_on()

            dt_last = None
            for timestart, timeend, mass in zip(self.agebins[:,0],self.agebins[:,1],self.percentMasses):

                dt_ = 10**timeend - 10**timestart

                if dt_last != None:
                    ax2.loglog(time_units*10**np.array([timestart, timestart]), [mass_last[1]/dt_last, mass[1]/dt_],
                               lw=0.7, c=c2)
                ax2.loglog(time_units*10**np.array([timestart,timeend]), [mass[1],mass[1]]/dt_, lw=0.7, c=c2)
                ax2.fill_between(time_units*10**np.array([timestart,timeend]), [mass[0],mass[0]]/dt_, [mass[2],mass[2]]/dt_,
                                 lw=0.7, color=c1, alpha=.5)

                mass_last = mass
                dt_last = dt_

            ax2.loglog()
            ax2.set_xlabel(r't $\rm ( Gyr )$', fontsize=fs2)
            ax2.set_ylabel(r'$\rm SFR_{model}$ $\rm (M_\odot \,\, yr^{-1})$   ', fontsize=fs2);
            ax2.set_xlim(time_units*10**7,time_units*10**self.t_uni)

            color_recent = colors[1]

            if figax == None:
                label='Last 100 Myr'
            else:
                label=''

            ax2.plot(np.linspace(0,time_units*10**8,100),[self.LocalSFRmed]*100, color=color_recent, lw=1,
                     label=label)
            ax2.fill_between(np.linspace(0,time_units*10**8,100), [self.LocalSFRlo]*100, [self.LocalSFRhi]*100,
                             color=color_recent, alpha=.1)

            ax2.legend(loc=4, fontsize=fs3, frameon=True)
            ax2.invert_xaxis()
            ax2.tick_params(which='major', direction='in', labelsize=fs3)
            ax2.tick_params(which='minor', direction='in', labelsize=fs3)
            ax2.tick_params(top=True, bottom=True, left=True, right=True, which='both')
            ax2.minorticks_on()

            returns = [fig, ax, ax2, ax3]

        if saveplots:
            plt.savefig(os.path.join(self.plots_dir,self.objstr+"_modelspec_phot_newshape_zoom.pdf"))

        return returns


    def plotSED_witherr_nuFnu(self, figax = None, saveplots=True, savedir='./', plotnuFnu=False,
                              plot_SFR=True, plotFilterNames=True,
                              fs = 21, fs2 = 17, fs3 = 12, fs_ticks = 12, figsize=(10,7),
                              colors = ['#DAA51B', '#38A6A5', '#2F8AC4', '#5F4690'],
                              # ax_space = (0, .1, 1, 1),
                              # ax2_space = (.7, .25, .35, .35),
                              # ax3_space = (0, 0, 1, .1),
                              ax_space = (.1, .2, .8, .7),
                              ax2_space = (.67, .57, .3, .3),
                              ax3_space = (.1, .1, .8, .1),
                              time_units=1e-9, c1='lightgrey', c2='darkgrey'):

        self.wphot = self.obs['wave_effective'].copy()

        self.a = 1.0 + self.model.params.get('zred', 0.0) # cosmological redshifting
        self.wspec = self.sps.wavelengths.copy()
        self.wspec *= self.a

        # If including SFR subplot:
        if plot_SFR:
            if figax is None:
                fig = plt.figure(figsize=figsize)
                ax=fig.add_axes(ax_space)
                ax2_space = np.array(ax2_space)
                ax2_space[2] = ax2_space[2]*figsize[1]/figsize[0]
                ax2=fig.add_axes(tuple(ax2_space))
                ax3=fig.add_axes(ax3_space)
            else:
                fig, [ax,ax2,ax3] = figax
        # If not including SFR subplot
        else:
            if figax is None:
                fig = plt.figure(figsize=figsize)
                ax=fig.add_axes(ax_space)
                ax3=fig.add_axes(ax3_space)
            else:
                fig, [ax,ax3] = figax

        color_model = colors[2]

        ax.tick_params(which='major', direction='in', labelsize=fs_ticks)
        ax.tick_params(which='minor', direction='in', labelsize=fs_ticks)
        ax.tick_params(top=True, bottom=True, left=True, right=True, which='both')
        ax.minorticks_on()
        ax.get_xaxis().set_visible(False)

        # REMOVED THIS ONE
        # ax.plot(self.wphot*self.wavelength_unitconv, self.mphot_map*self.flux_unitconv, alpha=0)

        self.n_obs = len(self.obs['filters'])
        self.obscolors = cm.rainbow_r(np.linspace(0,1,self.n_obs))
        # self.obscolors = cm.Set3(np.linspace(0,1,self.n_obs))


        if figax == None:
            label='Observed Photometry'
        else:
            label=''

        wave = self.wphot * self.wavelength_unitconv # micrometers
        self.nu = (const.c / (wave * u.micrometer)).to(u.Hz).value
        wave_spec = self.wspec*self.wavelength_unitconv
        self.nu_spec = (const.c / (wave_spec * u.micrometer)).to(u.Hz).value


        if plotnuFnu == True:
            self.flux_nu_factor = self.nu
            self.flux_nu_spec_factor = self.nu_spec
        else:
            self.flux_nu_factor = 1.
            self.flux_nu_spec_factor = 1.


        ax.errorbar(self.wphot*self.wavelength_unitconv, self.flux_nu_factor*self.obs['maggies']*self.flux_unitconv,
                    yerr=self.flux_nu_factor*self.obs['maggies_unc']*self.flux_unitconv,
                    ls='',  c='white',
                    markeredgecolor='black', markerfacecolor='none',
                    zorder=3)
        ax.scatter(self.wphot*self.wavelength_unitconv, self.flux_nu_factor*self.obs['maggies']*self.flux_unitconv,
                   c='white', marker='o', facecolor='none',
                   zorder=3)
        ax.errorbar(self.wphot*self.wavelength_unitconv, self.flux_nu_factor*self.obs['maggies']*self.flux_unitconv,
                    yerr=self.flux_nu_factor*self.obs['maggies_unc']*self.flux_unitconv,
                    ls='', color=self.obscolors, alpha=.7,
                    markeredgecolor='black', markerfacecolor='none',
                    zorder=3)
        ax.scatter(self.wphot*self.wavelength_unitconv, self.flux_nu_factor*self.obs['maggies']*self.flux_unitconv,
                   c=self.obscolors,  alpha=.7,
                   marker='o', facecolor='none', label=label,
                   zorder=3)

        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()


        if figax == None:
            label='Model Photometry'
        else:
            label=''

        asymmetric_error = np.array([self.mphot_siglo_err, self.mphot_sighi_err])
        ax.errorbar(self.wphot*self.wavelength_unitconv, self.flux_nu_factor*self.mphot_med*self.flux_unitconv,
                    yerr=self.flux_nu_factor*asymmetric_error*self.flux_unitconv,
                    c='slategrey', ecolor='slategrey', fmt='s', alpha=.7, label=label)

        if figax == None:
            label='$\pm1\sigma$ Model'
        else:
            label=''

        ax.plot(self.wspec*self.wavelength_unitconv, self.flux_nu_spec_factor*self.med*self.flux_unitconv,
                lw=0.7, c=c2, zorder=1)
        ax.fill_between(self.wspec*self.wavelength_unitconv,
                        self.flux_nu_spec_factor*self.siglo*self.flux_unitconv,
                        self.flux_nu_spec_factor*self.sighi*self.flux_unitconv,
                        label=label, lw=0.7, color=c1, alpha=.5, zorder=-1)

        ax.set_xlim(self.xmin_filters*self.wavelength_unitconv, self.xmax_filters*self.wavelength_unitconv)

        dy_filterspace = 10**-1.2*ymax-ymin # Changed this
        ax.set_ylim(ymin-dy_filterspace, ymax)

        # Plotting filter curves on twin axis
        axtwin = ax.twinx()

        obs_all = []
        # plot transmission curves
        for f,c in zip(self.obs['filters'],self.obscolors):
            w, t = f.wavelength.copy(), f.transmission.copy()
            t *= 1/t.max()
            t *=.06
            axtwin.plot(w*self.wavelength_unitconv, t, lw=1, color=c, alpha=.65)


            filtername = ' '.join(f.name.split('_'))
            if plotFilterNames:
                axtwin.text(f.wave_effective*self.wavelength_unitconv, max(t), filtername,
                            rotation=90, color=c, horizontalalignment="left", verticalalignment="bottom")

            obs_all.append(c)

        axtwin.set_ylim(0,1)
        axtwin.axes.get_yaxis().set_visible(False)

        ax.set_xlim(self.xmin_filters*self.wavelength_unitconv, self.xmax_filters*self.wavelength_unitconv)
        axtwin.set_xlim(self.xmin_filters*self.wavelength_unitconv, self.xmax_filters*self.wavelength_unitconv)

        ax.set_xlabel(r'Wavelength $\rm (\mu m)$', fontsize=fs)
        if plotnuFnu == True:
            ylabel = r'$\nu \,\, F_\nu$ $\rm (erg \,\, s^{-1} \,\, cm^{-2})$'
        else:
            ylabel = r'Flux Density $\rm (erg \,\, s^{-1} \,\, cm^{-2} \,\, Hz^{-1})$'
        ax.set_ylabel(ylabel, fontsize=fs)

        ax.legend(loc=2, fontsize=fs3, frameon=True)

        ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

        ax.tick_params(which='major', direction='in', labelsize=fs_ticks)
        ax.tick_params(which='minor', direction='in', labelsize=fs_ticks)
        ax.tick_params(top=True, bottom=True, left=True, right=True, which='both')
        ax.minorticks_on()

        # MAGNITUDE DIFFERENCE LOWER PLOT
        ax3.tick_params(which='major', direction='in', labelsize=fs_ticks)
        ax3.tick_params(which='minor', direction='in', labelsize=fs_ticks)
        ax3.tick_params(top=True, bottom=True, left=True, right=True, which='both')
        ax3.minorticks_on()

        # ax3.axhline(y=0, color='k', linestyle='--')

        self.diff_obs_mod = self.flux_nu_factor * (self.obs['maggies'] - self.mphot_med)  * self.flux_unitconv
        self.yerr_obs = self.flux_nu_factor * self.obs['maggies_unc'] * self.flux_unitconv
        self.yerr_mod = self.flux_nu_factor * np.abs(self.mphot_siglo - self.mphot_sighi) * self.flux_unitconv
        self.diff_obs_mod_err = np.sqrt(self.yerr_obs**2 + self.yerr_mod**2)

        # https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.chisquare.html

        ax3.plot(self.wphot*self.wavelength_unitconv, self.diff_obs_mod,
                 marker='o', alpha=0, linewidth=0)
        ax3.errorbar(self.wphot*self.wavelength_unitconv, self.diff_obs_mod,
                     yerr=self.diff_obs_mod_err,
                     ls='', color=self.obscolors, alpha=.7,
                     markeredgecolor='black', markerfacecolor='none',
                     zorder=3)
        ax3.scatter(self.wphot*self.wavelength_unitconv, self.diff_obs_mod,
                    c=self.obscolors, alpha=.7, marker='o', facecolor='none',
                    zorder=3)

        ax3_ylims = ax3.get_ylim()
        limit = max(np.abs(ax3_ylims))*1.3
        ax3.set_ylim(-limit,limit)

        ax3.set_xlim(self.xmin_filters*self.wavelength_unitconv,
                     self.xmax_filters*self.wavelength_unitconv)

        ax3.axhline(y=0, color='k', linestyle='--', alpha=0.2)

        ax3.set_xlabel(r'Wavelength $\rm (\mu m)$', fontsize=fs)
        ax3.set_ylabel(r'$\Delta$', fontsize=fs);

        returns = [fig, ax, ax3]

        if plot_SFR:

            ax2.tick_params(which='major', direction='in', labelsize=fs_ticks)
            ax2.tick_params(which='minor', direction='in', labelsize=fs_ticks)
            ax2.tick_params(top=True, bottom=True, left=True, right=True, which='both')
            ax2.minorticks_on()

            dt_last = None
            for timestart, timeend, mass in zip(self.agebins[:,0],self.agebins[:,1],self.percentMasses):

                dt_ = 10**timeend - 10**timestart

                if dt_last != None:
                    ax2.loglog(time_units*10**np.array([timestart, timestart]), [mass_last[1]/dt_last, mass[1]/dt_],
                               lw=0.7, c=c2)
                ax2.loglog(time_units*10**np.array([timestart,timeend]), [mass[1],mass[1]]/dt_, lw=0.7, c=c2)
                ax2.fill_between(time_units*10**np.array([timestart,timeend]), [mass[0],mass[0]]/dt_, [mass[2],mass[2]]/dt_,
                                 lw=0.7, color=c1, alpha=.5)

                mass_last = mass
                dt_last = dt_

            ax2.loglog()
            ax2.set_xlabel(r't $\rm ( Gyr )$', fontsize=fs2)
            ax2.set_ylabel(r'$\rm SFR_{model}$ $\rm (M_\odot \,\, yr^{-1})$   ', fontsize=fs2);
            ax2.set_xlim(time_units*10**7,time_units*10**self.t_uni)

            color_recent = colors[1]

            if figax == None:
                label='Last 100 Myr'
            else:
                label=''

            ax2.plot(np.linspace(0,time_units*10**8,100),[self.LocalSFRmed]*100, color=color_recent, lw=1,
                     label=label)
            ax2.fill_between(np.linspace(0,time_units*10**8,100), [self.LocalSFRlo]*100, [self.LocalSFRhi]*100,
                             color=color_recent, alpha=.1)

            ax2.legend(loc=4, fontsize=fs3, frameon=True)
            ax2.invert_xaxis()
            ax2.tick_params(which='major', direction='in', labelsize=fs3)
            ax2.tick_params(which='minor', direction='in', labelsize=fs3)
            ax2.tick_params(top=True, bottom=True, left=True, right=True, which='both')
            ax2.minorticks_on()

            returns = [fig, ax, ax2, ax3]

        if saveplots:
            plt.savefig(os.path.join(self.plots_dir,self.objstr+"_modelspec_phot_newshape_zoom.pdf"))

        return returns

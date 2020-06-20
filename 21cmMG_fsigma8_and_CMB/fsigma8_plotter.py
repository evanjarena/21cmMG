"""
We will use examine the effects that Horndeski models of 
modified gravity/dark energy have on the growth rate of 
structure fsigma8(z) from redshift z=6 to z=0.

Parameterization of background: 'lcdm'

Parameterization of gravity: 'hill_scale'
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

from peturbations import Peturbations


class Fsigma8:

    def __init__(self):
        self.f, self.ax = plt.subplots(2, 1, figsize=(15./1.5, 15./1.5), sharex='col')

    def reldev(self, lcdm, mg):
        """
        Calculate the relative percent deviation between LCDM 
        and non-LCDM (MG) outputs
        """
        return 100. * (mg - lcdm) / lcdm

    def chisqr(self, fs8_model, fs8_meas, fs8_meas_err):
        """
        Compute the chi squared
        Eq. 12 of arXiv:1712.04865
        """
        chisqr = 0
        for i in range(len(fs8_meas)):
            chisqr = chisqr + ((fs8_model[i] - fs8_meas[i])**2.) / (fs8_meas_err[i])**2.
        return chisqr

    def get_fsigma8_data(self, datafile):
        """
        Get fsigma8 data
        We list current measurements between redshifts 0 and 2
        from arXiv:1808.01337v3 Table II.
        Format: Index, z, fs8, fs8_error, #Dataset
        """
        data = np.loadtxt(datafile)
        meas_fsigma8 = data[:,2]
        meas_fsigma8_error = data[:,3]
        meas_zlist = data[:,1]
        index = data[:,0]
        return index, meas_zlist, meas_fsigma8, meas_fsigma8_error

    def make_plot(self):
        """
        Create two-panel plot for fsigma8 vs. z
        and rel. dev. from LCDM vs. z
        """
        self.ax[0].set_ylabel(r'$f\sigma_8$')
        self.ax[1].set_ylabel('$\mathrm{rel. dev. [\%]$}')
        self.ax[1].set_xlabel(r'$z$')

    def plot_errorbar(self, index, z, fsigma8, fsigma8_error):
        """
        Plot errorbar used for fsigma8 data.
        Each marker is a letter corresponding to 
        the name of the data set.  This is determined by
        the index number in arXiv:1808.01337v3 Table II.
        We will only be using a subset of this table
        (see datafile).
        """
        # WiggleZ
        if (index >= 10) and (index <= 13):
            marker = r'${\rm W}$'
        # 6dfGS
        elif (index == 13):
            marker = r'${\rm 6}$'
        # FastSound
        elif (index == 27):
            marker = r'${\rm F}$'
        # BOSS DR12
        elif (index >= 32) and (index <= 34):
            marker = r'${\rm B}$'
        # Vipers
        elif (index >= 41) and (index <= 42):
            marker = r'${\rm V}$'
        # SDSS-IV
        elif (index >= 60) and (index <= 63):
            marker = r'${\rm S}$'
        else:
            marker = '.'

        self.ax[0].errorbar(z, fsigma8, yerr=fsigma8_error, 
                            marker=marker, 
                            mfc='k', 
                            mew=0.5, mec='k', 
                            elinewidth=0.5, ecolor='k', 
                            capsize=5, capthick=0.5)

    def plot_fsigma8_data(self, fsigma8_datafile):
        """
        Plot errorbars for fsigma8 data
        """
        index, meas_zlist, meas_fsigma8, meas_fsigma8_error = self.get_fsigma8_data(fsigma8_datafile)
        for i in range(len(meas_zlist)):
            self.plot_errorbar(index[i], meas_zlist[i], meas_fsigma8[i], meas_fsigma8_error[i])

    def plot_fsigma8_GR(self, zlist, fsigma8_datafile):
        """
        Plot fsigma8 vs. z, rel. dev. from LCDM vs. z,
        and compute chi-squared from measurements
        for General Relativity.
        """
        # Get measured fsigma8 data
        __, meas_zlist, meas_fsigma8, meas_fsigma8_error = self.get_fsigma8_data(fsigma8_datafile)

        # fsigma8 for our redshift range
        lcdm_fsigma8 = Peturbations().fsigma8_GR(zlist)

        # fsigma8 for measured redshift range
        lcdm_fsigma8_meas_z = Peturbations().fsigma8_GR(meas_zlist)

        # Compute chi squared
        lcdm_chisqr = self.chisqr(lcdm_fsigma8_meas_z, meas_fsigma8, meas_fsigma8_error)

        # Plot fsigma8
        self.ax[0].plot(zlist, lcdm_fsigma8, 
                        dashes=[6,4], color='black', 
                        label=r'$\Lambda\mathrm{CDM}$', 
                        zorder=2)
        # Plot rel. dev.
        self.ax[1].plot(zlist, self.reldev(lcdm_fsigma8, lcdm_fsigma8), 
                        dashes=[6,4], color='black', 
                        label=r'$\chi^2={}$'.format(lcdm_chisqr), 
                        zorder=2)
        """
        # Save fsigma8 output to .dat file
        output = np.vstack((zlist, lcdm_fsigma8)).T
        np.savetxt('fsigma8_output/fsigma8_LCDM.dat', 
                   output, 
                   delimiter='    ', fmt='%0.4f', 
                   header=str('Redshift z, fsigma8(z)'))
        """

    def plot_fsigma8_MG(self, zlist, fsigma8_datafile, a_B, a_M, a_trans, r):
        """
        Plot fsigma8 vs. z, rel. dev. from LCDM vs. z,
        and compute chi-squared from measurements
        for Modified Gravity.
        """
        # Get measured fsigma8 data
        __, meas_zlist, meas_fsigma8, meas_fsigma8_error = self.get_fsigma8_data(fsigma8_datafile)
        
        # LCDM fsigma8
        lcdm_fsigma8 = Peturbations().fsigma8_GR(zlist)

        # MG fsigma8 for our redshift range
        mg_fsigma8 = Peturbations().fsigma8_MG(zlist, a_B, a_M, a_trans, r)
        
        # MG fsigma8 for measured redshift range
        mg_fsigma8_meas_z = Peturbations().fsigma8_MG(meas_zlist, a_B, a_M, a_trans, r)
        
        # Compute chi squared
        mg_chisqr = self.chisqr(mg_fsigma8_meas_z, meas_fsigma8, meas_fsigma8_error)

        # Plot fsigma8
        self.ax[0].plot(zlist, mg_fsigma8, 
                        label=r'$\hat{{\alpha}}_\mathrm{{B}}={}, \hat{{\alpha}}_\mathrm{{M}}={}, a_t=1/(1+{}), r={}$'.format(
                         a_B, a_M, (1./a_trans)-1., r), 
                        zorder=1) 
        # Plot rel. dev.
        self.ax[1].plot(zlist, self.reldev(lcdm_fsigma8, mg_fsigma8), 
                        label=r'$\chi^2={}$'.format(mg_chisqr), 
                        zorder=1)
        """
        # Save fsigma8 output to .dat file
        output = np.vstack((zlist, mg_fsigma8)).T
        np.savetxt('fsigma8_output/fsigma8_MG_a_B_'+str(a_B)+'_a_M_'+str(a_M)+'_a_trans_'+str(a_trans)+'_r_'+str(r)+'.dat', 
                   output, 
                   delimiter='    ', fmt='%0.4f', 
                   header=str('Redshift z, fsigma8(z)'))
        """

    def show_plot(self):
        """
        Show  and save two-panel plot for 
        fsigma8 vs. z
        and rel. dev. from LCDM vs. z
        """
        # Tight layout
        plt.tight_layout()
        # Remove whitespace between upper and lower plots
        plt.subplots_adjust(hspace=0, wspace=0.3) 
        # Tick marks on all sides of each plot and show legend
        for j in range(2):
            axes=self.ax[j]
            axes.tick_params(axis='both', which='both', direction='in',
                             top=True, right=True)
            legend=axes.legend(framealpha=0)
        # Save and show
        plt.savefig('fsigma8_LCDM_MG_measurements.pdf', format='pdf')
        plt.show()


















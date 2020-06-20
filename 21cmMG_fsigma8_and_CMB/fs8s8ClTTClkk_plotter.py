"""
We will use examine the effects that Horndeski models of 
modified gravity/dark energy have on fsigma8(z), sigma8(z),
the CMB lensing potential, and the CMB temperature power
spectrum.

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


class Fs8s8ClTTClkk:

    def __init__(self):
        self.f, self.ax = plt.subplots(2, 2, figsize=(10., 5.), sharex='col')

    def make_plot(self):
        """
        Create four-panel plot for 
          fsigma8 vs. z    C_ell^TT vs. ell
           sigma8 vs. z    C_ell^kappakappa vs. ell
        """
        self.ax[0, 0].set_ylabel(r'$f\sigma_8$')
        self.ax[1, 0].set_ylabel(r'$\sigma_8$')
        self.ax[1, 0].set_xlabel(r'$z$')

        self.ax[0, 1].set_ylabel(r'$[\ell(\ell+1)/2{{\pi}}]C_{{\ell}}^{{TT}}$')
        self.ax[1, 1].set_ylabel(r'$C_{{\ell}}^{{\kappa\kappa}}$')
        self.ax[1, 1].set_xlabel(r'$\ell$')

    def plot_fsigma8_GR(self, zlist):
        """
        Plot fsigma8 vs. z for General Relativity
        """
        # fsigma8 for our redshift range
        lcdm_fsigma8 = Peturbations().fsigma8_GR(zlist)

        # Plot fsigma8
        self.ax[0, 0].plot(zlist, lcdm_fsigma8, 
                           dashes=[6,4], color='black', 
                           label=r'$\Lambda\mathrm{CDM}$', 
                           zorder=2)

    def plot_fsigma8_MG(self, zlist, a_B, a_M, a_trans, r):
        """
        Plot fsigma8 vs. z for Modified Gravity.
        """
        # LCDM fsigma8
        #lcdm_fsigma8 = Peturbations().fsigma8_GR(zlist)

        # MG fsigma8 for our redshift range
        mg_fsigma8 = Peturbations().fsigma8_MG(zlist, a_B, a_M, a_trans, r)

        # Plot fsigma8
        #self.ax[0, 0].plot(zlist, mg_fsigma8, 
        #                   label=r'$\hat{{\alpha}}_\mathrm{{B}}={}, \hat{{\alpha}}_\mathrm{{M}}={}, a_t=1/(1+{}), r={}$'.format(
        #                       a_B, a_M, (1./a_trans)-1., r), 
        #                   zorder=1) 
        self.ax[0, 0].plot(zlist, mg_fsigma8, 
                           label=r'$\hat{{\alpha}}_\mathrm{{B}}={}, \hat{{\alpha}}_\mathrm{{M}}={}$'.format(a_B, a_M), 
                           zorder=1) 

    def plot_sigma8_GR(self, zlist):
        """
        Plot sigma8 vs. z for General Relativity
        """
        # sigma8 for our redshift range
        lcdm_sigma8 = Peturbations().sigma8_GR(zlist)

        # Plot sigma8
        self.ax[1, 0].plot(zlist, lcdm_sigma8, 
                           dashes=[6,4], color='black', 
                           zorder=2)

    def plot_sigma8_MG(self, zlist, a_B, a_M, a_trans, r):
        """
        Plot sigma8 vs. z for Modified Gravity.
        """
        # LCDM sigma8
        #lcdm_sigma8 = Peturbations().sigma8_GR(zlist)

        # MG sigma8 for our redshift range
        mg_sigma8 = Peturbations().sigma8_MG(zlist, a_B, a_M, a_trans, r)

        # Plot sigma8
        self.ax[1, 0].plot(zlist, mg_sigma8, 
                           zorder=1) 

    def plot_Cl_TT_GR(self, zlist, llist):
        """
        Plot C_ell^TT vs. ell for General Relativity
        """
        # C_ell^TT
        lcdm_Cl_TT = Peturbations().Cl_TT_GR(zlist, llist)

        # Plot C_ell^TT
        self.ax[0, 1].semilogx(llist, lcdm_Cl_TT, 
                             dashes=[6,4], color='black', 
                             zorder=2)

    def plot_Cl_TT_MG(self, zlist, llist, a_B, a_M, a_trans, r):
        """
        Plot C_ell^TT vs. ell for Modified Gravity
        """
        # LCDM C_ell^kappakappa
        #lcdm_Cl_TT = Peturbations().Cl_TT_GR(zlist, llist)

        # MG C_ell^kappakappa
        mg_Cl_TT = Peturbations().Cl_TT_MG(zlist, llist, a_B, a_M, a_trans, r)
        
        # Plot C_ell^kappakappa
        self.ax[0, 1].semilogx(llist, mg_Cl_TT, 
                             zorder=1) 

    def plot_Cl_kappa_kappa_GR(self, zlist, llist):
        """
        Plot C_ell^kappakappa vs. ell for General Relativity
        """
        # C_ell^kappakappa
        lcdm_Cl_kap_kap = Peturbations().Cl_kappa_kappa_GR(zlist, llist)

        # Plot C_ell^kappakappa
        self.ax[1, 1].semilogx(llist, lcdm_Cl_kap_kap, 
                             dashes=[6,4], color='black', 
                             zorder=2)

    def plot_Cl_kappa_kappa_MG(self, zlist, llist, a_B, a_M, a_trans, r):
        """
        Plot C_ell^kappakappa vs. ell for Modified Gravity
        """
        # LCDM C_ell^kappakappa
        #lcdm_Cl_kap_kap = Peturbations().Cl_kappa_kappa_GR(zlist, llist)

        # MG C_ell^kappakappa
        mg_Cl_kap_kap = Peturbations().Cl_kappa_kappa_MG(zlist, llist, a_B, a_M, a_trans, r)
        
        # Plot C_ell^kappakappa
        self.ax[1, 1].semilogx(llist, mg_Cl_kap_kap, 
                             zorder=1) 

    def show_plot(self):
        """
        Show  and save four-panel plot for 
        fsigma8 vs. z, sigma8 vs. z,
        C_ell^kappakappa, and C_ell^TT vs ell.
        """
        # Tight layout
        plt.tight_layout()
        # Remove whitespace between upper and lower plots
        plt.subplots_adjust(hspace=0, wspace=0.2) 
        # Tick marks on all sides of each plot and show legend
        for i in range(2):
            for j in range(2):
                axes=self.ax[i][j]
                axes.tick_params(axis='both', which='both', direction='in',
                                 top=True, right=True)
                legend=axes.legend(framealpha=0)
        # Save and show
        plt.savefig('fs8s8ClTTClkk_LCDM_MG_with_no_slip.pdf', format='pdf')
        plt.show()


















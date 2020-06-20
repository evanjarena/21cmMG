"""
We will use examine the effects that Horndeski models of 
modified gravity/dark energy have on the CMB lensing potential

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


class CMB_lensing_potential:

    def __init__(self):
        #self.peturb = Peturbations()
        #self.f, self.ax = plt.subplots(2, 1, figsize=(15/1.5, 10/1.5), sharex='col')
        self.f, self.ax = plt.subplots(2, 1, figsize=(15./1.5, 15./1.5), sharex='col')

    def reldev(self, lcdm, mg):
        """
        Calculate the relative percent deviation between LCDM 
        and non-LCDM (MG) outputs
        """
        return 100. * (mg - lcdm) / lcdm

    def make_plot(self):
        """
        Create two-panel plot for C_ell^kappakappa vs. ell
        and rel. dev. from LCDM vs. ell
        """
        self.ax[0].set_ylabel(r'$C_{{\ell}}^{{\kappa\kappa}}$')
        self.ax[1].set_ylabel('$\mathrm{rel. dev. [\%]$}')
        self.ax[1].set_xlabel(r'$\ell$')

    def plot_Cl_kappa_kappa_GR(self, llist):
        """
        s
        """
        # Cl_kappa_kappa
        lcdm_Cl_kap_kap = Peturbations().Cl_kappa_kappa_GR(llist)

        # Plot Cl_kappa_kappa
        self.ax[0].loglog(llist, lcdm_Cl_kap_kap, 
                          dashes=[6,4], color='black', 
                          label=r'$\Lambda\mathrm{CDM}$', 
                          zorder=2)
        # Plot rel. dev.
        self.ax[1].semilogx(llist, self.reldev(lcdm_Cl_kap_kap, lcdm_Cl_kap_kap), 
                            dashes=[6,4], color='black', 
                            zorder=2)

    def plot_Cl_kappa_kappa_MG(self, llist, a_B, a_M, a_trans, r):
        """
        s
        """
        # LCDM Cl_kappa_kappa
        lcdm_Cl_kap_kap = Peturbations().Cl_kappa_kappa_GR(llist)

        # MG Cl_kappa_kappa
        mg_Cl_kap_kap = Peturbations().Cl_kappa_kappa_MG(llist, a_B, a_M, a_trans, r)
        
        # Plot Cl_kappa_kappa
        self.ax[0].loglog(llist, mg_Cl_kap_kap, 
                          label=r'$\hat{{\alpha}}_\mathrm{{B}}={}, \hat{{\alpha}}_\mathrm{{M}}={}, a_t=1/(1+{}), r={}$'.format(
                              a_B, a_M, (1./a_trans)-1., r), 
                          zorder=1) 
        # Plot rel. dev.
        self.ax[1].semilogx(llist, self.reldev(lcdm_Cl_kap_kap, mg_Cl_kap_kap), 
                            zorder=1)

    def show_plot(self):
        """
        Show  and save two-panel plot for 
        C_ell^kappakappa vs. ell
        and rel. dev. from LCDM vs. ell
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
        plt.savefig('CMB_lensing_potential_LCDM_MG.pdf', format='pdf')
        plt.show()


















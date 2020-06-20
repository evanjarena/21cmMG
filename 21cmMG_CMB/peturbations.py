import numpy as np
from classy import Class

class Peturbations:

    def __init__(self):
        self.cosmo = Class()

    def cosmo_Cl_GR(self):
        """
        Compute cosmology for General Relativity 
        with Planck15 parameters.
        """
        #zstr=','.join(map(str,zlist+zlist[-1]+2))
        lcdmpars = {'output': 'tCl, pCl, lCl',
                    'lensing': 'yes',
                    #'P_k_max_h/Mpc':20,
                    #'z_pk': zstr,
                    'background_verbose': 1, #Info
                    'tau_reio': 0.07,
                    'omega_cdm': 0.11987,     
                    'A_s': 2.204e-9, 
                    'h': 0.6715918, #computed from 100*theta_s=1.042143 
                    'N_ur': 3.046-1.,  
                    'N_ncdm': 1.,           
                    'm_ncdm': 0.06,       
                    'omega_b': 0.022252,     
                    'n_s': 0.96475, 
                   }

        self.cosmo.set(lcdmpars)
        self.cosmo.compute()

    def cosmo_Cl_MG(self, alpha_b, alpha_m, a_trans, r):
        """
        Compute cosmology for Modified gravity with a 
        LCDM background expansion (with Planck15
        parameters) and 'hill_scale' parameterization 
        for the property functions.
        """
        #zstr=','.join(map(str,zlist+zlist[-1]+2))
        mgpars =   {'output': 'tCl, pCl, lCl',
                    'lensing': 'yes',
                    'background_verbose': 1, #Info
                    'tau_reio': 0.07,
                    'omega_cdm': 0.11987,     
                    'A_s': 2.204e-9, 
                    'h': 0.6715918, #computed from 100*theta_s=1.042143 
                    'N_ur': 3.046-1.,  
                    'N_ncdm': 1.,           
                    'm_ncdm': 0.06,       
                    'omega_b': 0.022252,     
                    'n_s': 0.96475, 

                    'Omega_Lambda': 0,
                    'Omega_fld': 0,
                    'Omega_smg': -1,
                    'expansion_model': 'lcdm',
                    'gravity_model': 'hill_scale',
                    }
        
        # Planck Mass
        Mstarsq = 1.
        # Property Functions
        alpha_t = 0.
        alpha_k = 0.001

        mgpars['parameters_smg'] = "{}, {}, {}, {}, {}, {}, {}".format(
                                   Mstarsq, 
                                   alpha_k, alpha_b, alpha_m, alpha_t,
                                   a_trans, r)
        self.cosmo.set(mgpars)
        self.cosmo.compute()

    def Cl_kappa_kappa_GR(self, llist):
        """
        Get C_ell^kappakappa for General Relativity
        """
        self.cosmo_Cl_GR() 
        lmax = llist[-1]
        Cl_phi_phi_GR = self.cosmo.lensed_cl(lmax)['pp'][2:]
        Cl_kap_kap_GR = (1./4.) * ((llist * (llist + 1.))**2.) * Cl_phi_phi_GR
        return Cl_kap_kap_GR

    def Cl_kappa_kappa_MG(self, llist, alpha_b, alpha_m, a_trans, r):
        """
        Get C_ell^kappakappa for Modified Gravity with a LCDM background 
        expansion and hill_scale parameterization for the property
        functions.
        """
        self.cosmo_Cl_MG(alpha_b, alpha_m, a_trans, r)
        lmax = llist[-1]
        Cl_phi_phi_MG = self.cosmo.lensed_cl(lmax)['pp'][2:]
        Cl_kap_kap_MG = (1./4.) * ((llist * (llist + 1.))**2.) * Cl_phi_phi_MG
        return Cl_kap_kap_MG

